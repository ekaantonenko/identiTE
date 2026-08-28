#!/usr/bin/env bash
# _addMeth_to_flanks.sh
#
# Intersects TE-flank BED files with per-genome methylation data and outputs
# a BED file with methylation sums and cytosine context counts appended.
#
# Usage:
#   _addMeth_to_flanks.sh \
#       --ids       <ids_path>    \   # one genome ID per line
#       --flanks    <flanks_path> \   # dir containing <genome>.flanks.bed
#       --meth      <meth_path>   \   # dir containing methylation files
#       --out       <out_path>    \   # output root directory
#       --meth_format <fmt>       \   # fast5_freq | bedGraph.gz
#       --meth_suffix <suffix>    \   # see _prep_meth.sh for naming conventions
#       [--threads  <N>]              # parallelism (default: 4)
#
# Output per genome:
#   <out_path>/FlanksPerGenome_meth/<genome>.flanks.meth.bed
#
# Columns:
#   1-6  : original flank BED columns
#   7    : methylation frequency sum  (CG)
#   8    : cytosine count             (CG)
#   9    : methylation frequency sum  (CHG)
#   10   : cytosine count             (CHG)
#   11   : methylation frequency sum  (CHH)
#   12   : cytosine count             (CHH)


################################################################################################
################################### Argument parsing ###########################################
################################################################################################

usage() {
    cat <<EOF
Usage: $0 --ids <ids_path> --flanks <flanks_path> --meth <meth_path> \\
          --out <out_path> --meth_format <fmt> --meth_suffix <suffix> \\
          [--threads <N>]
 
  --ids         File with one genome ID per line
  --flanks      Directory containing <genome>.flanks.bed files
  --meth        Directory containing methylation files
  --out         Output root directory
  --meth_format Input methylation format: fast5_freq | bedGraph.gz
  --meth_suffix Filename suffix after the genome ID (use {context} placeholder
                for bedGraph.gz format, e.g.
                  fast5_freq  : "fast5s.C.call_mods.frequency.tsv.gz"
                  bedGraph.gz : "clean_bismark_bt2.deduplicated.sorted.{context}.bedGraph.gz")
  --threads     Parallel genome workers (default: 4)
EOF
    exit 1
}
 
ids_path=""
flanks_path=""
meth_path=""
out_path=""
meth_format=""
meth_suffix=""
threads=4
 
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ids)          ids_path="$2";    shift 2 ;;
        --flanks)       flanks_path="$2"; shift 2 ;;
        --meth)         meth_path="$2";   shift 2 ;;
        --out)          out_path="$2";    shift 2 ;;
        --meth_format)  meth_format="$2"; shift 2 ;;
        --meth_suffix)  meth_suffix="$2"; shift 2 ;;
        --threads)      threads="$2";     shift 2 ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done
 
if [[ -z "$ids_path"   || -z "$flanks_path" || -z "$meth_path" || \
      -z "$out_path"   || -z "$meth_format" || -z "$meth_suffix" ]]; then
    echo "ERROR: missing required argument(s)." >&2
    usage
fi




################################################################################################
################################### Source helpers #############################################
################################################################################################
 
# Resolve the directory this script lives in so _prep_meth.sh is found
# regardless of the caller's working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_prep_meth.sh"

mkdir -p "$out_path/FlanksPerGenome_meth"

################################################################################################
################################### Helpers ####################################################
################################################################################################

# do_intersects WDIR GENOME METHBASE
#
# Intersects the genome's flank BED with each per-context methylation BED.
# -loj  : left outer join  (unmatched flanks get a "." row so nothing is lost)
# -wa   : print the original flank interval
# -wb   : append the matching methylation interval
#
# Output: $wdir/$genome.intersect.<CTX>.txt
do_intersects() {
    local wdir="$1"
    local genome="$2"
    local methbase="$3"   # full path prefix, e.g. $wdir/$genome
 
    local context
    for context in CG CHG CHH; do
        bedtools intersect -loj -wa -wb \
            -a  "$wdir/${genome}.bed" \
            -b  "${methbase}.${context}.bed" \
            > "$wdir/${genome}.intersect.${context}.txt"
    done
}
export -f do_intersects

# do_groupby WDIR GENOME MODE OP FREQ_COL
#
# Groups the intersect output by the 6 flank columns and aggregates the
# methylation frequency column.
#
# WDIR      – per-genome scratch directory
# GENOME    – genome ID
# MODE      – output label: "meth_sum" | "context"
# OP        – bedtools groupby -o operation: "sum" | "count"
# FREQ_COL  – 1-based column index of the frequency value in the intersect file
#             (= number of flank BED columns + 5, because the meth BED has
#              chrom/start/end/label/freq in cols 1-5 relative to the appended block)
do_groupby() {
    local wdir="$1"
    local genome="$2"
    local mode="$3"
    local op="$4"
    local freq_col="$5"
 
    local context
    for context in CG CHG CHH; do
        bedtools groupby \
            -i  "$wdir/${genome}.intersect.${context}.txt" \
            -grp 1-6 \
            -c  "$freq_col" \
            -o  "$op" \
            | cut -f 7 \
            > "$wdir/${genome}.${mode}.${context}.txt"
    done
}
export -f do_groupby




################################################################################################
################################### Per-genome worker ##########################################
################################################################################################
process_genome() {
    local GENOME="$1"
    local flanks_path="$2"
    local meth_path="$3"
    local out_path="$4"
    local meth_format="$5"
    local meth_suffix="$6"
 
    local flankdir="$out_path/FlanksPerGenome_meth"
 
    # Per-genome scratch directory avoids filename collisions under parallelism
    local wdir="$flankdir/$GENOME"
    mkdir -p "$wdir"
 
    # All methylation intermediates live under the per-genome scratch dir
    local methbase="$wdir/$GENOME.meth"
 
    echo "[$$] Processing $GENOME..."
    # ── Prepare flanks BED ──────────────────────────────────────────────────
    # Append _RagTag to every chromosome name (col 1)
    cp "$flanks_path/$GENOME.flanks.bed" "$wdir/$GENOME.bed"

    
    # Count flank BED columns so we can find the frequency column after intersect.
    # bedtools intersect -loj -wb appends the full meth BED row; freq is its col 5,
    # so the merged column index = flank_cols + 5.
    local flank_cols
    flank_cols=$(awk 'NR==1{print NF; exit}' "$wdir/$GENOME.bed")
    local freq_col=$(( flank_cols + 5 ))


 
    # ── Decompress, filter, split methylation into per-context BEDs ─────────
    echo "[$$] $GENOME — preparing methylation ($meth_format)"
    prep_meth "$GENOME" "$meth_path" "$meth_format" "$meth_suffix" "$methbase" || {
        echo "ERROR [$$]: prep_meth failed for $GENOME" >&2
        return 1
    }


 
    # ── Intersect flanks × methylation ───────────────────────────────────────
    echo "[$$] $GENOME — intersecting TE flanks with methylation"
    do_intersects "$wdir" "$GENOME" "$methbase"


 
    # ── Group-by: sum methylation frequency per flank ────────────────────────
    echo "[$$] $GENOME — summing methylation"
    do_groupby "$wdir" "$GENOME" "meth_sum" "sum" "$freq_col"
 
    # ── Group-by: count covered cytosines per flank ──────────────────────────
    echo "[$$] $GENOME — counting cytosine contexts"
    do_groupby "$wdir" "$GENOME" "context" "count" "$freq_col"
 
    # ── Paste methylation sums and contexts into a single BED ────────────────
    paste "$wdir/$GENOME.bed" \
          "$wdir/$GENOME.meth_sum.CG.txt" \
          "$wdir/$GENOME.context.CG.txt" \
          "$wdir/$GENOME.meth_sum.CHG.txt" \
          "$wdir/$GENOME.context.CHG.txt" \
          "$wdir/$GENOME.meth_sum.CHH.txt" \
          "$wdir/$GENOME.context.CHH.txt" \
        > "$flankdir/$GENOME.flanks.meth.bed"

 
    # ── Clean up per-genome scratch dir (uncomment when satisfied) ───────────
    rm -rf "$wdir"
 
    echo "[$$] $GENOME — done → $flankdir/$GENOME.flanks.meth.bed"
}
export -f process_genome
 
 
################################################################################################
################################### Parallel main loop #########################################
################################################################################################
 
# Export all options so worker subshells can see them
export flanks_path meth_path out_path meth_format meth_suffix
 
# Each xargs worker sources _prep_meth.sh itself (process_genome is exported,
# but the functions it calls from _prep_meth.sh must also be available in the
# subshell).  We pass SCRIPT_DIR explicitly via the environment.
export SCRIPT_DIR
 
xargs -a "$ids_path" -I {} -P "$threads" \
    bash -c '
        source "$SCRIPT_DIR/_prep_meth.sh"
        process_genome "$@"
    ' _ {} "$flanks_path" "$meth_path" "$out_path" "$meth_format" "$meth_suffix"
 
echo "All genomes processed."


