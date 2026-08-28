#!/usr/bin/env bash
# _addGenesMeth_to_TE.sh
#
# Adds gene annotation (2 upstream and 2 downstream).
# Intersects TE-flank BED files with per-genome methylation data and outputs
# a BED file with methylation sums and cytosine context counts appended.
#
# Usage:
#   _addMeth_to_flanks.sh \
#       --ids       <ids_path>    \   # one genome ID per line
#       --te        <te_path>     \   # dir containing <genome>.bed
#       --gene      <gene_path>   \   # dir containing gene annotations
#       --meth      <meth_path>   \   # dir containing methylation files
#       --out       <out_path>    \   # output root directory
#       --meth_format <fmt>       \   # fast5_freq | bedGraph.gz
#  --gene_format  Gene annotation format: liftoff | ensembl | bed10
#  --gene_suffix  Filename suffix after the genome ID for gene annotation files
#                 e.g. "_LiftOver.TAIR10_genes.bed" or ".genes.gff3"
#       --meth_suffix <suffix>    \   # see _prep_meth.sh for naming conventions
#       [--threads  <N>]              # parallelism (default: 4)
#
# Output per genome:
#   <out_path>/TEperGenome_meth/<genome>.meth.bed
#
# Columns:
#   1-7    : original TE BED columns
#   8-17   : upstream gene annotation   (2 genes × name + strand + start + end + distance)
#   18-27  : downstream gene annotation (2 genes × name + strand + start + end + distance)
#   then, for each of the 9 regions (up.inside up.0_100 up.100_200 up.200_300
#                                     down.inside down.0_100 down.100_200 down.200_300 te):
#     +2 cols : meth_sum (CG),  context_count (CG)
#     +2 cols : meth_sum (CHG), context_count (CHG)
#     +2 cols : meth_sum (CHH), context_count (CHH)
#   last 3 cols : genome-wide average methylation (CG, CHG, CHH)
 


################################################################################################
################################### Argument parsing ###########################################
################################################################################################

usage() {
    cat <<EOF
Usage: $0 --ids <ids_path> --te <te_path> --gene <gene_path> --meth <meth_path> \\
          --out <out_path> --meth_format <fmt> --meth_suffix <suffix> \\
          [--threads <N>]
 
  --ids         File with one genome ID per line
  --te          Directory containing <genome>.bed files
  --gene        Directory containing gene annotation files
  --meth        Directory containing methylation files
  --out         Output root directory
  --meth_format Input methylation format: fast5_freq | bedGraph.gz
  --meth_suffix Filename suffix after the genome ID (use {context} placeholder
                for bedGraph.gz format, e.g.
                  fast5_freq  : "fast5s.C.call_mods.frequency.tsv.gz"
                  bedGraph.gz : "clean_bismark_bt2.deduplicated.sorted.{context}.bedGraph.gz")
  --gene_format  Gene annotation format: liftoff | ensembl | bed10
  --gene_suffix  Filename suffix after the genome ID for gene annotation files
                 e.g. "_LiftOver.TAIR10_genes.bed" or ".genes.gff3"
  --threads     Parallel genome workers (default: 4)
EOF
    exit 1
}
 
ids_path=""
te_path=""
gene_path=""
gene_format=""
gene_suffix=""
meth_path=""
out_path=""
meth_format=""
meth_suffix=""
threads=4
 
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ids)          ids_path="$2";    shift 2 ;;
        --te)           te_path="$2";     shift 2 ;;
        --gene)         gene_path="$2";   shift 2 ;;
        --gene_format)  gene_format="$2"; shift 2 ;;
        --gene_suffix)  gene_suffix="$2"; shift 2 ;;
        --meth)         meth_path="$2";   shift 2 ;;
        --out)          out_path="$2";    shift 2 ;;
        --meth_format)  meth_format="$2"; shift 2 ;;
        --meth_suffix)  meth_suffix="$2"; shift 2 ;;
        --threads)      threads="$2";     shift 2 ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done
 
if [[ -z "$ids_path"  || -z "$te_path"   || -z "$gene_path" || \
      -z "$gene_format" || -z "$gene_suffix" || \
      -z "$meth_path" || -z "$out_path"  || \
      -z "$meth_format" || -z "$meth_suffix" ]]; then
    echo "ERROR: missing required argument(s)." >&2
    usage
fi


################################################################################################
################################### Source helpers #############################################
################################################################################################
 
# Resolve the directory this script lives in so _prep_meth.sh is always found,
# regardless of the caller's working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_prep_meth.sh"
source "$SCRIPT_DIR/_prep_genes.sh"
 
mkdir -p "$out_path/TEperGenome_meth"


################################################################################################
################################### Helpers ####################################################
################################################################################################

# Regions iterated in every context loop (exported so subshells see it)
REGIONS_STR="up.inside up.0_100 up.100_200 up.200_300 down.inside down.0_100 down.100_200 down.200_300 te"
export REGIONS_STR
 

 

# make_windows WDIR GENOME
#
# Derives 8 window BED files from the TE body ($wdir/$genome.bed).
# All windows carry the same 7 columns as the source BED.
make_windows() {
    local wdir="$1"
    local genome="$2"
    local copy="$wdir/$genome.bed"
 
    # Upstream inside: first 100 bp of TE body from the start
    awk 'BEGIN{OFS="\t"} {print $1, $2, $2+100, $4, $5, $6, $7}' \
        "$copy" > "$wdir/$genome.up.inside.bed" 
 
    # Upstream flanks: 0–100, 100–200, 200–300 bp before TE start
    for i in 0 100 200; do
        local j=$(( i + 100 ))
        awk -v i="$i" -v j="$j" \
            'BEGIN{OFS="\t"} {
                start = ($2-j > 0 ? $2-j : 0)
                end   = ($2-i > 0 ? $2-i : 0)
                print $1, start, end, $4, $5, $6, $7
            }' "$copy" > "$wdir/$genome.up.${i}_${j}.bed" 
    done
 
    # Downstream inside: last 100 bp of TE body before the end
    awk 'BEGIN{OFS="\t"} {print $1, ($3-100 > 0 ? $3-100 : 0), $3, $4, $5, $6, $7}' \
        "$copy" > "$wdir/$genome.down.inside.bed" 
 
    # Downstream flanks: 0–100, 100–200, 200–300 bp after TE end
    for i in 0 100 200; do
        local j=$(( i + 100 ))
        awk -v i="$i" -v j="$j" \
            'BEGIN{OFS="\t"} {print $1, $3+i, $3+j, $4, $5, $6, $7}' \
            "$copy" > "$wdir/$genome.down.${i}_${j}.bed" 
    done
}
export -f make_windows

#: '
#make_windows() {
#    local wdir=$1
#    local genome=$2
#    local copy="$wdir/$genome.bed"
#
#
#    # All 8 window files are independent — launch in parallel
#    # Upstream inside (first 100 bp of TE body from the start)
#    awk 'BEGIN{OFS="\t"}{print $1,$2,$2+100,$4,$5,$6,$7}' \
#        "$copy" > "$wdir/$genome.up.inside.bed"
#
#    # Upstream flanks: 0–100, 100–200, 200–300 bp before TE start
#    for i in 0 100 200; do
#        local j=$((i+100))
#        awk -v i=$i -v j=$j 'BEGIN{OFS="\t"}{
#            start = $2-j > 0 ? $2-j : 0
#            end   = $2-i > 0 ? $2-i : 0
#            print $1, start, end, $4,$5,$6,$7
#        }' "$copy" > "$wdir/$genome.up.${i}_${j}.bed" 
#    done
#
#    # Downstream inside (last 100 bp of TE body before the end)
#    awk 'BEGIN{OFS="\t"} {print $1,($3-100>0?$3-100:0),$3,$4,$5,$6,$7}' \
#        "$copy" > "$wdir/$genome.down.inside.bed" 
# 
#    # Downstream flanks: 0–100, 100–200, 200–300 bp after TE end
#    for i in 0 100 200; do
#        local j=$((i+100))
#        awk -v i=$i -v j=$j 'BEGIN{OFS="\t"} {
#            print $1, $3+i, $3+j, $4,$5,$6,$7
#        }' "$copy" > "$wdir/$genome.down.${i}_${j}.bed" 
#    done
#}        
#export -f make_windows
#'

# do_intersects WDIR GENOME METHBASE
#
# Intersects each region BED with each per-context methylation BED.
# -loj : left outer join (unmatched TEs get a "." row so nothing is lost)
# -wa  : print the original region interval
# -wb  : append the matching methylation interval
#
# Output: $wdir/$genome.intersect.<context>.<region>.txt
do_intersects() {
    local wdir="$1"
    local genome="$2"
    local methbase="$3"   # full path prefix, e.g. $wdir/$genome.meth
 
    local -a regions
    read -ra regions <<< "$REGIONS_STR"
 
    local context region input
    for context in CG CHG CHH; do
        for region in "${regions[@]}"; do
            if [[ "$region" == "te" ]]; then
                input="$wdir/$genome.bed"
            else
                input="$wdir/$genome.$region.bed"
            fi
            bedtools intersect -loj -wa -wb \
                -a "$input" \
                -b "${methbase}.${context}.bed" \
                > "$wdir/$genome.intersect.$context.$region.txt"
        done
    done
}
export -f do_intersects

: '
### VERSION TO UPDATE
do_intersects() {
    local wdir=$1
    local genome=$2
    local wdir=$3   # directory holding the per-context .bed files for this genome
    local -a regions
    read -ra regions <<< "$REGIONS_STR"   # rebuild array from exported string

 
    local context region input
    for context in CG CHG CHH; do
        for region in "${regions[@]}"; do
            # TE body uses the copy bed; windows have their own files
            if [[ "$region" == "te" ]]; then
                input="$wdir/$genome.bed"
            else
                input="$wdir/$genome.$region.bed"
            fi
            bedtools intersect -loj -wa -wb \
                -a "$input" \
                -b "$wdir/$genome.fast5s.C.call_mods.frequency.$context.bed" \
                > "$wdir/$genome.intersect.$context.$region.txt" 
        done
    done
}
export -f do_intersects
'

# do_groupby WDIR GENOME MODE OP FREQ_COL
#
# Groups intersect output by the 7 TE/window columns and aggregates the
# methylation frequency column.
#
# MODE     – output label: "meth_sum" | "context"
# OP       – bedtools groupby -o operation: "sum" | "count"
# FREQ_COL – 1-based column index of the frequency value in the intersect file
#            (= number of window BED columns + 5, because the meth BED has
#             chrom/start/end/label/freq in cols 1-5 relative to the appended block)
do_groupby() {
    local wdir="$1"
    local genome="$2"
    local mode="$3"
    local op="$4"
    local freq_col="$5"
 
    local -a regions
    read -ra regions <<< "$REGIONS_STR"
 
    local context region
    for context in CG CHG CHH; do
        for region in "${regions[@]}"; do
            bedtools groupby \
                -i  "$wdir/$genome.intersect.$context.$region.txt" \
                -grp 1-7 \
                -c  "$freq_col" \
                -o  "$op" \
                | cut -f 8 \
                > "$wdir/$genome.$mode.$context.$region.txt" 
        done
    done
}
export -f do_groupby

: '
do_groupby() {
    local wdir=$1
    local genome=$2
    local mode=$3   # "meth_sum" or "context"
    local op=$4     # bedtools -o operation: "sum" or "count"
    local -a regions
    read -ra regions <<< "$REGIONS_STR"   # rebuild array from exported string
 

 
    local context region
    for context in CG CHG CHH; do
        for region in "${regions[@]}"; do
            bedtools groupby \
                -i "$wdir/$genome.intersect.$context.$region.txt" \
                -grp 1-7 \
                -c 12 \
                -o "$op" \
                | cut -f 8 \
                > "$wdir/$genome.$mode.$context.$region.txt" 
        done
    done
}
export -f do_groupby
'

# join_meth WDIR GENOME CONTEXT REGION
#
# Pastes meth_sum and context_count for one context × region pair.
join_meth() {
    local wdir="$1"
    local genome="$2"
    local context="$3"
    local region="$4"
    paste \
        "$wdir/$genome.meth_sum.$context.$region.txt" \
        "$wdir/$genome.context.$context.$region.txt" \
        > "$wdir/$genome.$context.$region.meth.txt"
}
export -f join_meth


# finalize_region WDIR GENOME REGION
#
# Pastes CG + CHG + CHH meth files for one region into a single file.
finalize_region() {
    local wdir="$1"
    local genome="$2"
    local region="$3"
    paste \
        "$wdir/$genome.CG.$region.meth.txt" \
        "$wdir/$genome.CHG.$region.meth.txt" \
        "$wdir/$genome.CHH.$region.meth.txt" \
        > "$wdir/$genome.TEs.meth.$region.txt"
}
export -f finalize_region







################################################################################################
################################### Per-genome worker ##########################################
################################################################################################
process_genome() {
    local GENOME=$1
    local te_path=$2
    local gene_path=$3
    local gene_format="$4"
    local gene_suffix="$5"
    local meth_path="$6"
    local out_path="$7"
    local meth_format="$8"
    local meth_suffix="$9"
 
    local TEdir="$out_path/TEperGenome_meth"
    
    # Per-genome scratch directory avoids filename collisions under parallelism
    local wdir="$TEdir/$GENOME"
    mkdir -p "$wdir"

    # Use wdir for ALL files, including methylation intermediates
    local methbase="$wdir/$GENOME.meth"

    local -a regions
    read -ra regions <<< "$REGIONS_STR"   # rebuild array from exported string
 

    echo "[$$] Processing $GENOME..."
 
    # ── Prepare TE BED ───────────────────────────────────────────────────────
    cp "$te_path/$GENOME.bed" "$wdir/$GENOME.bed"

    # Count TE BED columns so we can find the frequency column after intersect.
    # bedtools intersect -loj -wb appends the full meth BED row; freq is its col 5,
    # so the merged column index = te_cols + 5.
    local te_cols
    te_cols=$(awk 'NR==1{print NF; exit}' "$wdir/$GENOME.bed")
    local freq_col=$(( te_cols + 5 ))
 
    # ── Decompress, filter, split methylation into per-context BEDs ─────────
    echo "[$$] $GENOME — preparing methylation ($meth_format)"
    prep_meth "$GENOME" "$meth_path" "$meth_format" "$meth_suffix" "$methbase" || {
        echo "ERROR [$$]: prep_meth failed for $GENOME" >&2
        return 1
    }


    # ── Genome-wide average methylation (chromosomes Chr1–Chr5 only) ─────────
    echo "[$$] $GENOME — computing genome-wide averages"
    local averageCG averageCHG averageCHH
    averageCG=$(awk  '$1~/^Chr[1-5]$/ {sum+=$5; n++} END{printf "%.6f\n", (n>0?sum/n:0)}' \
                    "$methbase.CG.bed")
    averageCHG=$(awk '$1~/^Chr[1-5]$/ {sum+=$5; n++} END{printf "%.6f\n", (n>0?sum/n:0)}' \
                    "$methbase.CHG.bed")
    averageCHH=$(awk '$1~/^Chr[1-5]$/ {sum+=$5; n++} END{printf "%.6f\n", (n>0?sum/n:0)}' \
                    "$methbase.CHH.bed")
 

    # ── Gene annotation: normalise to 5-column BED ──────────────────────────
    echo "[$$] $GENOME — preparing gene annotation ($gene_format)"
    prep_genes "$GENOME" "$gene_path" "$gene_format" "$gene_suffix" \
               "$wdir/$GENOME.genes.cut.bed" || {
        echo "ERROR [$$]: prep_genes failed for $GENOME" >&2
        return 1
    }
 
    ### FIX GENE SUFFIXES
    # ── Gene annotation: keep only gene features, slim to 5 columns ─────────
    #awk '$8=="gene"' "$gene_path/${GENOME}_LiftOver.TAIR10_genes.bed" \
    #    | awk -F'\t' -v OFS='\t' '{
    #        s=$10; sub(/;.+/,"",s); sub(/ID=/,"",s)
    #        print $1,$2,$3,s,$6
    #      }' \
    #    > "$wdir/$GENOME.genes.cut.bed"


 
    # ── Closest gene up/down: both directions  ────────────────────
    for dir in down up; do
        (
        bedtools closest \
            -a "$wdir/$GENOME.bed" \
            -b "$wdir/$GENOME.genes.cut.bed" \
            -D ref -i"${dir::1}" -t first -k 2 \
            | bedtools groupby \
                -grp 1-7 \
                -c 11,9,10,12,13 \
                -o collapse \
            | awk -F'[,\t]' -v OFS='\t' '{
                print $1,$2,$3,$4,$5,$6,$7,
                      $8,$10,$12,$14,$16,
                      $9,$11,$13,$15,$17
                }' \
            > "$wdir/$GENOME.te_genes.$dir.bed" || true
        ) 
    done


    #local averageCG  averageCHG  averageCHH
    #averageCG=$(< "$wdir/$GENOME.avg.CG.txt")
    #averageCHG=$(< "$wdir/$GENOME.avg.CHG.txt")
    #averageCHH=$(< "$wdir/$GENOME.avg.CHH.txt")

    # ── Merge upstream + downstream gene columns ─────────────────────────────
    paste "$wdir/$GENOME.te_genes.up.bed" "$wdir/$GENOME.te_genes.down.bed" \
        | cut -f 1-17,25- \
        > "$wdir/$GENOME.te_genes.updown.bed"
 
 
     # ── Windows (8 awk jobs) ───────────────────────────
    echo "[$$] $GENOME — creating windows"
    make_windows "$wdir" "$GENOME"
 
    # ── Intersect (27 bedtools jobs) ──────────────────────────────
    echo "[$$] $GENOME — intersecting"
    do_intersects "$wdir" "$GENOME" "$methbase"
 
    # ── Groupby sum + count (54 jobs total) ─────────
    echo "[$$] $GENOME — groupby"
    do_groupby "$wdir" "$GENOME" "meth_sum" "sum" "$freq_col"
    do_groupby "$wdir" "$GENOME" "context"  "count" "$freq_col"
    

    # ── join + finalize: all 9 regions  ────────────────────────────
    for region in "${regions[@]}"; do
        (
            join_meth    "$wdir" "$GENOME" CG  "$region"
            join_meth    "$wdir" "$GENOME" CHG "$region"
            join_meth    "$wdir" "$GENOME" CHH "$region"
            finalize_region "$wdir" "$GENOME" "$region"
        ) 
    done

    # ── Accumulate columns across regions (serial: each paste depends on previous) ──
    echo "[$$] $GENOME — accumulating region columns"
    local tmp="$wdir/$GENOME.te_genes.updown.bed"
    for region in "${regions[@]}"; do
        paste "$tmp" "$wdir/$GENOME.TEs.meth.$region.txt" \
            > "$wdir/$GENOME.TEs.genes.meth.$region.txt"
        tmp="$wdir/$GENOME.TEs.genes.meth.$region.txt"
    done
    
 
    # ── Append genome-wide averages ──────────────────────────────────────────
    awk -v cg="$averageCG" -v chg="$averageCHG" -v chh="$averageCHH" \
        'BEGIN{OFS="\t"} {print $0, cg, chg, chh}' \
        "$wdir/$GENOME.TEs.genes.meth.te.txt" \
        > "$TEdir/$GENOME.meth.bed"
 

 
 
    # ── Cleanup ──────────────────────────────────────────────────────────────
    rm -rf "$wdir"   # entire per-genome scratch dir
 
    echo "[$$] $GENOME — done → $TEdir/$GENOME.meth.bed"
}
export -f process_genome






################################################################################################
################################### Parallel main loop #########################################
################################################################################################

# Export all options so worker subshells can see them
export te_path gene_path meth_path out_path meth_format meth_suffix REGIONS_STR
 
# Each xargs worker sources _prep_meth.sh itself (process_genome is exported,
# but the functions it calls from _prep_meth.sh must also be available in the
# subshell). We pass SCRIPT_DIR explicitly via the environment.
export SCRIPT_DIR
 
xargs -a "$ids_path" -I {} -P "$threads" \
    bash -c '
        source "$SCRIPT_DIR/_prep_meth.sh"
        source "$SCRIPT_DIR/_prep_genes.sh"
        process_genome "$@"
    ' _ {} "$te_path" "$gene_path" "$gene_format" "$gene_suffix" \
           "$meth_path" "$out_path" "$meth_format" "$meth_suffix"
 
echo "All genomes processed."