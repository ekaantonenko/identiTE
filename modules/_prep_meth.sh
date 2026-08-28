#!/usr/bin/env bash
# _prep_meth.sh
# Methylation preprocessing helpers.
# Source this file; do NOT execute it directly.
#
# Public entry point:
#   prep_meth <genome> <meth_path> <meth_format> <meth_suffix> <methbase>
#
#   genome      – genome ID (e.g. "sample_01")
#   meth_path   – directory containing per-genome methylation data
#   meth_format – one of: fast5_freq | bedGraph.gz
#   meth_suffix – filename suffix after the genome ID, with "{context}" as a
#                 literal placeholder for the context token where applicable.
#
#                 fast5_freq  : no {context} needed, single file per genome
#                   e.g. "fast5s.C.call_mods.frequency.tsv.gz"
#                   → full path: <meth_path>/<genome>.<meth_suffix>
#
#                 bedGraph.gz : must contain the literal string "{context}",
#                   which is replaced at runtime with CG / CHG / CHH
#                   e.g. "clean_bismark_bt2.deduplicated.sorted.{context}.bedGraph.gz"
#                   → full path: <meth_path>/<genome>.<meth_suffix[{context}→CG]>  etc.
#
#   methbase    – full path prefix for output BED files
#                   e.g. "$wdir/$GENOME"
#                 Outputs written:
#                   $methbase.CG.bed
#                   $methbase.CHG.bed
#                   $methbase.CHH.bed
#                 (A unified $methbase.bed is also written by fast5_freq as an
#                  intermediate; it is empty/absent for bedGraph.gz.)
#
# Output BED columns (all formats, 6-col):
#   1=chrom  2=start(0-based)  3=end  4=label(dot)  5=freq(0-1)  6=context_string


export LC_NUMERIC=C

################################################################################################
################################### Format handlers ############################################
################################################################################################

# ── fast5 / ONT frequency TSV ───────────────────────────────────────────────────────────────
# File: <meth_path>/<genome>.<meth_suffix>
# Input columns (1-based, tab-separated):
#   1=chrom  2=pos  3=strand_string  4=pos_on_+_strand  5..8=misc
#   9=coverage  10=freq(0-1)  11=context_string
# Filter: coverage (col 9) >= 5
# Output BED (6-col): chrom, pos, pos+1, strand_string, freq, context_string
_prep_meth_fast5_freq() {
    local genome="$1"
    local meth_path="$2"
    local meth_suffix="$3"
    local methbase="$4"

    local infile="$meth_path/${genome}.${meth_suffix}"

    if [[ ! -f "$infile" ]]; then
        echo "ERROR [_prep_meth_fast5_freq]: input file not found: $infile" >&2
        return 1
    fi

    gunzip -c "$infile" \
        | awk 'BEGIN{OFS="\t"}
               NR==1 && $1=="chrom" { next }   # skip header if present
               $9 >= 5 { print $1, $2, $2+1, $3, $10, $11 }' \
        > "$methbase.bed"
}


# ── Bismark bedGraph (three separate context files) ──────────────────────────────────────────
# Files: <meth_path>/<genome>.<meth_suffix with {context} replaced by CG|CHG|CHH>
# Input columns (1-based, tab-separated):
#   1=chrom  2=start(0-based)  3=end
#   4=methylation%  5=count_methylated  6=count_unmethylated
# Filter: total coverage (col5+col6) >= 5
# Output BED (6-col): chrom, start, end, ".", freq(0-1), context_string
#
# context_string values assigned here are compatible with the grep patterns in
# _split_meth_contexts (they are not actually used for splitting here since we
# write directly to per-context files, but they keep the schema consistent).
_prep_meth_bismark_bedgraph() {
    local genome="$1"
    local meth_path="$2"
    local meth_suffix="$3"   # must contain the literal token {context}
    local methbase="$4"

    local min_cov=5

    # Bismark context label → output file suffix
    declare -A context_out=( [CG]="CG" [CHG]="CHG" [CHH]="CHH" )

    # A synthetic context string that satisfies the CG/CHG/CHH grep patterns in
    # _split_meth_contexts, in case that function is ever called on these files.
    declare -A context_str=(
        [CG]="NNCGN"
        [CHG]="NNCHGN"
        [CHH]="NNCHHN"
    )

    local context
    for context in CG CHG CHH; do
        # Replace the {context} placeholder with the actual context token
        local suffix_resolved="${meth_suffix/\{context\}/$context}"
        local infile="$meth_path/${genome}.${suffix_resolved}"
        local outfile="$methbase.${context_out[$context]}.bed"

        if [[ ! -f "$infile" ]]; then
            echo "ERROR [_prep_meth_bismark_bedgraph]: input file not found: $infile" >&2
            return 1
        fi

        gunzip -c "$infile" \
            | awk -v min_cov="$min_cov" -v cs="${context_str[$context]}" \
                'BEGIN{OFS="\t"}
                 NR==1 && $1=="track" { next }   # skip UCSC track header if present
                 {
                     cov = $5 + $6
                     if (cov >= min_cov) {
                         freq = $4 / 100
                         print $1, $2, $3, ".", freq, cs
                     }
                 }' \
            > "$outfile"
    done

    # Touch a zero-byte unified .bed so downstream code that tests -s does not
    # try to split it.
    : > "$methbase.bed"
}


################################################################################################
################################### Context splitter ###########################################
################################################################################################

# Split a unified .bed (written by fast5_freq) into per-context files.
# Skipped automatically when $methbase.bed is empty (bismark_bedgraph already
# wrote the per-context files directly).
#
# Context-string regexes match column 6 of the 6-col output BED.
# Patterns (strand-aware, matching the ONT context_string field):
#   CG  : + strand  …NNCGN   / - strand  …NGCNN
#   CHG : + strand  …NNCHGN  / - strand  …NGHNN  (H = A/T/C)
#   CHH : + strand  …NNCHHN  / - strand  …NHHNN
_split_meth_contexts() {
    local methbase="$1"

    # Nothing to split if the unified .bed is empty (bismark_bedgraph path)
    if [[ ! -s "$methbase.bed" ]]; then
        return 0
    fi

    # Use process substitution + tee for a single-pass split
    #tee \
    #    >(grep -E "(\+.*[ATCG][ATCG]CG[ATCG])|(\-.*[ATCG]GC[ATCG][ATCG])"       \
    #        > "$methbase.CG.bed")  \
    #    >(grep -E "(\+.*[ATCG][ATCG]C[ATC]G)|(\-.*G[ATC]C[ATCG][ATCG])"          \
    #        > "$methbase.CHG.bed") \
    #    >(grep -E "(\+.*[ATCG][ATCG]C[ATC][ATC])|(\-.*[ATC][ATC]C[ATCG][ATCG])"  \
    #        > "$methbase.CHH.bed") \
    #    < "$methbase.bed" > /dev/null
    # wait is needed because tee's process substitutions are async
    #wait

    tee \
    >(grep -E "[ATCG]{2}CG[ATCG]"        > "$methbase.CG.bed")  \
    >(grep -E "[ATCG]{2}C[ACT]G"         > "$methbase.CHG.bed") \
    >(grep -E "[ATCG]{2}C[ACT][ACT]"     > "$methbase.CHH.bed") \
    < "$methbase.bed" > /dev/null
    wait
}


################################################################################################
################################### Public entry point #########################################
################################################################################################

prep_meth() {
    local genome="$1"
    local meth_path="$2"
    local meth_format="$3"
    local meth_suffix="$4"
    local methbase="$5"

    case "$meth_format" in
        fast5_freq)
            _prep_meth_fast5_freq       "$genome" "$meth_path" "$meth_suffix" "$methbase" || return 1
            ;;
        bedGraph.gz)
            _prep_meth_bismark_bedgraph "$genome" "$meth_path" "$meth_suffix" "$methbase" || return 1
            ;;
        *)
            echo "ERROR [prep_meth]: unknown meth_format '$meth_format'" >&2
            echo "       Supported formats: fast5_freq, bedGraph.gz" >&2
            return 1
            ;;
    esac

    _split_meth_contexts "$methbase"

    # Strip _RagTag suffix from chromosome names in all per-context BEDs
    for context in CG CHG CHH; do
        local f="$methbase.$context.bed"
        [[ -f "$f" ]] && sed -i 's/_RagTag//' "$f"
    done
}

export -f prep_meth \
           _prep_meth_fast5_freq \
           _prep_meth_bismark_bedgraph \
           _split_meth_contexts