#!/usr/bin/env bash
# _prep_genes.sh
# Gene annotation preprocessing helpers.
# Source this file; do NOT execute it directly.
#
# Public entry point:
#   prep_genes <genome> <gene_path> <gene_format> <gene_suffix> <outfile>
#
#   genome      – genome ID (e.g. "sample_01")
#   gene_path   – directory containing gene annotation files
#   gene_format – one of: liftoff | ensembl | bed10
#   gene_suffix – filename suffix after the genome ID
#                   e.g. "_LiftOver.TAIR10_genes.gff3"  for liftoff / ensembl
#                        "_LiftOver.TAIR10_genes.bed"   for bed10
#                 Full path resolved as: <gene_path>/<genome><gene_suffix>
#   outfile     – path for the normalised output BED
#
# Output BED (5-col, tab-separated):
#   1=chrom  2=start(0-based)  3=end  4=gene_id  5=strand


################################################################################################
################################### Format handlers ############################################
################################################################################################

# ── Liftoff GFF3 ─────────────────────────────────────────────────────────────────────────────
# Header: 3 non-data lines at the top (skipped by NR<=3).
# Feature type in col 3; gene ID in col 9 as the first semicolon-delimited
# field with an "ID=" prefix, e.g. "ID=AT1G01010;Parent=...".
# Coordinates are 1-based closed → converted to 0-based half-open BED.
_prep_genes_liftoff() {
    local infile="$1"
    local outfile="$2"
    awk 'BEGIN{OFS="\t"}
         NR<=3 { next }
         $3=="gene" {
             id=$9; sub(/;.*/,"",id); sub(/^ID=/,"",id)
             print $1, $4-1, $5, id, $7
         }' "$infile" > "$outfile"
}

# ── Ensembl GFF3 ─────────────────────────────────────────────────────────────────────────────
# Header: comment lines starting with "#" (skipped by /^#/).
# Same column layout as Liftoff, but the ID field uses an optional "gene:"
# prefix, e.g. "ID=gene:AT1G01010" → stripped to "AT1G01010".
# Coordinates are 1-based closed → converted to 0-based half-open BED.
_prep_genes_ensembl() {
    local infile="$1"
    local outfile="$2"
    awk 'BEGIN{OFS="\t"}
         /^#/ { next }
         $3=="gene" {
             id=$9; sub(/;.*/,"",id); sub(/^ID=/,"",id); sub(/^gene:/,"",id)
             print $1, $4-1, $5, id, $7
         }' "$infile" > "$outfile"
}

# ── BED10 ─────────────────────────────────────────────────────────────────────────────────────
# Columns: chr(1) start(2) end(3) type_n(4) .(5) strand(6) method(7) type(8) .(9) attributes(10)
# Feature type in col 8 (not col 3); gene ID in col 10 as "ID=<id>[;...]".
# Coordinates are already 0-based half-open — no adjustment needed.
# Chromosome names may carry a _RagTag suffix that is stripped.
_prep_genes_bed10() {
    local infile="$1"
    local outfile="$2"
    awk 'BEGIN{OFS="\t"}
         $8=="gene" {
             chr=$1; sub(/_RagTag$/,"",chr)
             id=$10; sub(/;.*/,"",id); sub(/^ID=/,"",id)
             print chr, $2, $3, id, $6
         }' "$infile" > "$outfile"
}


################################################################################################
################################### Public entry point #########################################
################################################################################################

prep_genes() {
    local genome="$1"
    local gene_path="$2"
    local gene_format="$3"
    local gene_suffix="$4"
    local outfile="$5"

    local infile="$gene_path/${genome}${gene_suffix}"

    if [[ ! -f "$infile" ]]; then
        echo "ERROR [prep_genes]: input file not found: $infile" >&2
        return 1
    fi

    case "$gene_format" in
        liftoff)  _prep_genes_liftoff "$infile" "$outfile" ;;
        ensembl)  _prep_genes_ensembl "$infile" "$outfile" ;;
        bed10)    _prep_genes_bed10   "$infile" "$outfile" ;;
        *)
            echo "ERROR [prep_genes]: unknown gene_format '$gene_format'" >&2
            echo "       Supported formats: liftoff, ensembl, bed10" >&2
            return 1
            ;;
    esac
}

export -f prep_genes _prep_genes_liftoff _prep_genes_ensembl _prep_genes_bed10