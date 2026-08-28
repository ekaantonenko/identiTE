#!/usr/bin/env python3
import argparse
import os

import warnings
warnings.filterwarnings("ignore")

import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess

import datetime

import numpy as np
import pandas as pd
import math


from io import StringIO
#from Bio.Blast.Applications import NcbiblastnCommandline
from pybedtools import BedTool

from _formats import (
    load_library_graffite, get_info_parser,
    load_library_bed,
    read_fasta,
    load_gene_gff,
    load_pericentromeres,
    process_graffite, process_bed
)


# =========================
# FASTA UTILITIES
# =========================
def load_one_genome(idx, seq_path, seq_suffix):
    return idx, read_fasta(f'{seq_path}/{idx}{seq_suffix}')





# =========================
# REMOVE PERICENTROMERES
# =========================
def removePeri(df, PERI, 
                chr_col='CHROM',
                start_col='POS',
                end_col='POS'
                ):
    dfs_noperi = []

    ### Chromosomes 1, 2, 3, 5 ###
    for i in [0,1,2,5]:
        tmp = df[(df[chr_col] == PERI.loc[i,'Chr']) & 
                ((df[end_col] < PERI.loc[i,'start']) | 
                 (df[start_col] > PERI.loc[i,'end']))]
        dfs_noperi.append(tmp)

    ### Chromosome 4 (as pericentromeres consist of two regions) ###
    tmp = df[(df[chr_col] == 'Chr4') & 
            ((df[end_col] < PERI.loc[3,'start']) | # up- from the 1st part
             (df[start_col] > PERI.loc[4,'end']) | # down- from the 2nd part
            ((df[start_col] > PERI.loc[3,'end']) & 
              (df[end_col] < PERI.loc[4,'start'])) #between the 1st and 2nd parts
            )]
    dfs_noperi.append(tmp)
        
    df_noperi = pd.concat(dfs_noperi)
    #df_noperi = df_noperi.sort_values(by=['Chr', start_col, 'end_ref'])
    return df_noperi

# =========================
# REF GENES RETRIEVAL
# =========================
_ref_genes_cache = None
_ref_lock = threading.Lock()

def get_ref_genes():
    global _ref_genes_cache
    with _ref_lock:
        if _ref_genes_cache is None:
            genes = pd.read_csv(gene_ref_path, sep='\t', skiprows=1,
                                header=None, encoding='cp1252')
            genes = genes[genes[2] == 'gene']
            genes = genes[genes[0].isin(['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5'])]
            genes = genes[[0,3,4,6,8]]
            genes.columns = ['Chr', 'start', 'end', 'strand', 'id']
            genes[['start', 'end']] = genes[['start', 'end']].astype(int)
            genes['id'] = genes['id'].str.split(';').str[0].str.lstrip('ID=')
            genes = genes.sort_values(by=['Chr', 'start', 'end'])
            _ref_genes_cache = genes
    return _ref_genes_cache

# =========================
# ADD CLOSEST GENES
# =========================
def addGenesRef(df):
    genes_ref = get_ref_genes()  # free if already loaded
    df = df.dropna(subset=['end_ref'])
    
    df['end_ref'] = df['end_ref'].astype(int)
    original_index = df.index.copy()

    df = df.sort_values(['Chr', 'start_ref', 'end_ref']).reset_index(drop=True)

    df_bed = BedTool.from_dataframe(df) 
    genes_ref_bed = BedTool.from_dataframe(genes_ref)

    te_genes_down = df_bed.closest(genes_ref_bed, D='ref', iu=True, t='first', k=2)

    te_genes_down_collapse = te_genes_down.groupby(g=[i for i in range(1,len(df.columns)+1)], 
                                                   c=[i for i in range(len(df.columns)+1,len(df.columns)+7)], 
                                                   o=['collapse'])

    te_genes_up = df_bed.closest(genes_ref_bed, D='ref', id=True, t='first', k=2) 
    te_genes_up_collapse = te_genes_up.groupby(g=[i for i in range(1,len(df.columns)+1)], 
                                               c=[i for i in range(len(df.columns)+1,len(df.columns)+7)], 
                                               o=['collapse'])

    te_genes_updown = te_genes_down_collapse.to_dataframe(header=None)
    te_genes_updown[[i for i in range(len(df.columns)+6,len(df.columns)+12)]] = te_genes_up_collapse.to_dataframe(header=None)[[i for i in range(len(df.columns),len(df.columns)+6)]]

    for i in range (len(te_genes_updown.columns)-12, len(te_genes_updown.columns)):
        te_genes_updown[i+12] = te_genes_updown[i].astype(str).str.split(',').str[1]
        te_genes_updown[i] = te_genes_updown[i].astype(str).str.split(',').str[0]

    te_genes_updown.columns = [*df.columns,
                            'Gene1.down.Chr', 'Gene1.down.start', 'Gene1.down.end', 'Gene1.down.strand', 'Gene1.down.id', 'Gene1.down.dist',
                            'Gene1.up.Chr', 'Gene1.up.start', 'Gene1.up.end', 'Gene1.up.strand', 'Gene1.up.id', 'Gene1.up.dist',
                            'Gene2.down.Chr', 'Gene2.down.start', 'Gene2.down.end', 'Gene2.down.strand', 'Gene2.down.id', 'Gene2.down.dist',
                            'Gene2.up.Chr', 'Gene2.up.start', 'Gene2.up.end', 'Gene2.up.strand', 'Gene2.up.id', 'Gene2.up.dist']
    te_genes_updown.index = original_index
    return te_genes_updown
##############################



BLAST_COLS = ['query', 'subject', 'pc_identity', 'aln_length', 'mismatches',
              'gaps_opened', 'query_start', 'query_end', 'subject_start',
              'subject_end', 'e_value', 'bitscore']

RESULT_SUFFIXES = ['btwgenes_pc_identity', 'btwgenes_start', 'btwgenes_end',
                   'left_pc_identity', 'left_start', 'left_end',
                   'right_pc_identity', 'right_start', 'right_end',
                   'btwflanks_pc_identity', 'btwflanks_start', 'btwflanks_end']


############################################################################################
############################################################################################
############################################################################################

def process_one_genome(idx, i, row, genes_t, leftflank, rightflank, TE,
                       genome_cache, gene_cache, out_path, worker_id):
    """All BLAST work for one (TE, genome) pair. Returns dict of results."""
    local_results = {}
    

    ### neighboring genes (2 upstream & 2 downstream) ###    
    genes_idx = gene_cache[idx]
    matched = genes_idx[genes_idx['id'].isin(genes_t)]
    if matched.empty:
        return local_results

    left  = matched['start_ref'].min()
    right = matched['end_ref'].max()
    subseq = genome_cache[idx][row['Chr']]
    subseq_between_genes = subseq[left:right+1]
    if not subseq_between_genes:
        return local_results

    path_db       = f'{out_path}/blast/TE_db_{worker_id}_{idx}_{i}.fsa'
    path_combined = f'{out_path}/blast/combined_{worker_id}_{idx}_{i}.fa'
    blastout      = f'{out_path}/blast/combined_blast_{worker_id}_{idx}_{i}.tab'


    with open(path_db, 'w') as f:
        f.write(f'>{idx}_{i}\n{subseq_between_genes}\n')
    with open(path_combined, 'w') as f:
        f.write(f'>TE_{i}\n{TE}\n')
        f.write(f'>leftflank_{i}\n{leftflank}\n')
        f.write(f'>rightflank_{i}\n{rightflank}\n')

    result = subprocess.run([
        "blastn",
        "-query", path_combined,
        "-subject", path_db,
        "-out", blastout,
        "-outfmt", "6"
    ], check=False, capture_output=True, text=True)

    #print("STDOUT:", result.stdout)
    #print("STDERR:", result.stderr)
    #print("Return code:", result.returncode)


    left_start = left_end = right_start = right_end = None

    if os.path.exists(blastout) and os.path.getsize(blastout) > 0:
        results = pd.read_csv(blastout, sep='\t', header=None, names=BLAST_COLS)

        te_hit    = results[results['query'] == f'TE_{i}']
        left_hit  = results[results['query'] == f'leftflank_{i}']
        right_hit = results[results['query'] == f'rightflank_{i}']

        if len(te_hit) > 0 and te_hit.iloc[0]['aln_length'] > 50:
            local_results[f'{idx}_btwgenes_pc_identity'] = te_hit.iloc[0]['pc_identity']
            local_results[f'{idx}_btwgenes_start']       = te_hit.iloc[0]['subject_start'] + left
            local_results[f'{idx}_btwgenes_end']         = te_hit.iloc[0]['subject_end']   + left

        if len(left_hit) > 0 and left_hit.iloc[0]['aln_length'] > 50:
            left_start = left_hit.iloc[0]['subject_start'] + left
            left_end   = left_hit.iloc[0]['subject_end'] + left
            local_results[f'{idx}_left_pc_identity'] = left_hit.iloc[0]['pc_identity']
            local_results[f'{idx}_left_start']        = left_start 
            local_results[f'{idx}_left_end']          = left_end  

        if len(right_hit) > 0 and right_hit.iloc[0]['aln_length'] > 50:
            right_start = right_hit.iloc[0]['subject_start'] + left
            right_end   = right_hit.iloc[0]['subject_end'] + left
            local_results[f'{idx}_right_pc_identity'] = right_hit.iloc[0]['pc_identity']
            local_results[f'{idx}_right_start']        = right_start
            local_results[f'{idx}_right_end']          = right_end

    # Between-flanks BLAST — only if both flanks found
    if all(v is not None for v in [left_start, left_end, right_start, right_end]):
        a = (left_start, left_end)
        b = (right_start, right_end)
        lf, rf = (a, b) if min(a) < min(b) else (b, a)
        i_left_flank  = max(lf)
        i_right_flank = min(rf)

        if i_right_flank - i_left_flank >= 50:
            subseq_between_flanks = subseq[i_left_flank:i_right_flank+1]
            path_te_db  = f'{out_path}/blast/TE_db_flanks_{worker_id}_{idx}_{i}.fsa'
            path_te_query = f'{out_path}/blast/TE_query_{worker_id}_{idx}_{i}.fa'
            blastout_te = f'{out_path}/blast/te_blast_{worker_id}_{idx}_{i}.tab'

            with open(path_te_db, 'w') as f:
                f.write(f'>{idx}_{i}\n{subseq_between_flanks}\n')
            with open(path_te_query, 'w') as f:
                f.write(f'>TE_{i}\n{TE}\n')

            subprocess.run([
                "blastn",
                "-query", path_te_query,
                "-subject", path_te_db,
                "-out", blastout_te,
                "-outfmt", "6"
                ], check=True)

            if os.path.exists(blastout_te) and os.path.getsize(blastout_te) > 0:
                results_te = pd.read_csv(blastout_te, sep='\t', header=None, names=BLAST_COLS)
                te_hit2 = results_te[results_te['query'] == f'TE_{i}']
                if len(te_hit2) > 0 and te_hit2.iloc[0]['aln_length'] > 50:
                    local_results[f'{idx}_btwflanks_pc_identity'] = te_hit2.iloc[0]['pc_identity']
                    local_results[f'{idx}_btwflanks_start']       = te_hit2.iloc[0]['subject_start'] + i_left_flank
                    local_results[f'{idx}_btwflanks_end']         = te_hit2.iloc[0]['subject_end']   + i_left_flank

            for p in [path_te_db, path_te_query, blastout_te]:
                try: os.unlink(p)
                except OSError: pass

    for p in [path_db, path_combined, blastout]:
        try: os.unlink(p)
        except OSError: pass

    return local_results



############################################################################################
############################################################################################
############################################################################################
def BLAST_flanks_and_TE(df,
                        out_path, i_batch,
                        seq_ref, genome_cache, gene_cache,
                        n_workers):

    
    if df.shape[0] > 0:
        df = addGenesRef(df)


    if df.shape[0] == 0:
        #df.to_csv(f'{out_path}/Result_batch_{i_batch}.csv')
        return df
    
    result_cols = {f'{idx}_{suffix}': {} for idx in IDS for suffix in RESULT_SUFFIXES}
    
    records = df.to_dict('index')

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {}
        ### loop over TEs in Ref
        for i, row in records.items():
            genes_t    = {row['Gene1.up.id'], row['Gene2.up.id'],
                          row['Gene1.down.id'], row['Gene2.down.id']} - {'.', np.nan}
            leftflank  = seq_ref[row['Chr']][row['start_ref']-301:row['start_ref']-1]
            rightflank = seq_ref[row['Chr']][row['end_ref']:row['end_ref']+300]
            TE         = row['SEQ']

            for idx in IDS:
                worker_id = f'{i}_{idx}'  # globally unique
                future = executor.submit(
                    process_one_genome,
                    idx, i, row, genes_t, leftflank, rightflank, TE,
                    genome_cache, gene_cache, out_path, worker_id
                )
                futures[future] = (i, idx)

        for future in as_completed(futures):
            i, idx = futures[future]
            try:
                local_results = future.result()
                for key, val in local_results.items():
                    result_cols[key][i] = val
            except Exception as e:
                print(f'Error TE {i} genome {idx}: {e}')

    for col, vals in result_cols.items():
        df[col] = pd.Series(vals)

    return df
 

############################################################################################
############################################################################################
############################################################################################


def resultsToBinary(df_results, IDS):

    df_binary = pd.DataFrame(columns = IDS, index = df_results['ID'])

    for idx in IDS:
        left_min  = df_results[[f'{idx}_left_start',  f'{idx}_left_end']].min(axis=1)
        left_max  = df_results[[f'{idx}_left_start',  f'{idx}_left_end']].max(axis=1)
        right_min = df_results[[f'{idx}_right_start', f'{idx}_right_end']].min(axis=1)
        right_max = df_results[[f'{idx}_right_start', f'{idx}_right_end']].max(axis=1)

        flank_5p_start = np.where(left_min <= right_min, left_min,  right_min)
        flank_5p_end   = np.where(left_min <= right_min, left_max,  right_max)
        flank_3p_start = np.where(left_min <= right_min, right_min, left_min)
        flank_3p_end   = np.where(left_min <= right_min, right_max, left_max)
 
        span_start = np.minimum(flank_5p_start, flank_3p_start)  # same as flank_5p_start, but safe
        span_end   = np.maximum(flank_5p_end,   flank_3p_end)
        gap = flank_3p_start.astype(float) - flank_5p_end.astype(float)
        
        df_binary[idx] = np.where(gap < 10, 
                                  0, # confirmed absence
                                  np.nan
                                  )

        df_binary[idx] = np.where(df_results[f'{idx}_btwflanks_pc_identity'] > 90.0, 
                                    1, # confirmed presence
                                    df_binary[idx])

    return df_binary



############################################################################################
####################################  MAIN  ################################################
############################################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--library_path", required=True, type=str)
    parser.add_argument("--library_format", default='graffite',
                        choices=['graffite', 'bed'])
    parser.add_argument("--ids_path", required=True, type=str)
    parser.add_argument("--gene_path", required=True, type=str)
    parser.add_argument("--gene_ref_path", required=True, type=str)
    parser.add_argument("--gene_suffix", default=None, type=str)

    parser.add_argument("--seq_path", required=True, type=str)
    parser.add_argument("--seq_ref_path", required=True, type=str)
    parser.add_argument("--seq_suffix", required=True, type=str)
    parser.add_argument("--peri_path", required=True, type=str)
    parser.add_argument("--out_path", required=True, type=str)

    parser.add_argument("--batch_first", required=True, type=int)
    parser.add_argument("--batch_last", required=True, type=int)
    parser.add_argument("--batch_size", required=True, type=int)
    parser.add_argument("--n_threads", default=None, type=int)

    parser.add_argument("--gene_format",  default='liftoff',
                    choices=['liftoff', 'ensembl', 'bed10'],
                    help="Gene annotation file format")
    parser.add_argument("--vcf_version",  default=None,
                        help="GraffiTE VCF version (auto-detected if not specified)")
    parser.add_argument("--peri_format",  default='bed3',
                        choices=['bed3', 'bed4'])

    args = parser.parse_args()

    library_path = args.library_path
    library_format = args.library_format
    ids_path = args.ids_path
    gene_path = args.gene_path
    gene_ref_path = args.gene_ref_path
    gene_suffix = args.gene_suffix
    seq_path = args.seq_path
    seq_ref_path = args.seq_ref_path
    seq_suffix = args.seq_suffix
    peri_path = args.peri_path

    out_path = args.out_path
    
    batch_first = args.batch_first
    batch_last = args.batch_last
    batch_size = args.batch_size
    n_threads = args.n_threads

    gene_format = args.gene_format
    vcf_version = args.vcf_version
    peri_format = args.peri_format

    ###########################################################################

    IDS = list(pd.read_csv(ids_path, header=None)[0])

    print(f"Loading reference genome... : {str(datetime.datetime.now())}")
    seq_ref = read_fasta(seq_ref_path)

    print(f"Loading gene annotations... : {str(datetime.datetime.now())}")
    gene_cache = {}
    for idx in IDS:
        gene_cache[idx] = load_gene_gff(f'{gene_path}/{idx}{gene_suffix}', format=args.gene_format)

    print(f"Loading library... : {str(datetime.datetime.now())}")
    if library_format=='graffite':
        library_full, vcf_version = load_library_graffite(library_path)
        #info_parser = get_info_parser(vcf_version)
    elif library_format=='bed':
        library_full = load_library_bed(library_path)

    print(len(library_full))
    print(library_full.head())
    

    print(f"Loading and removing pericentromeres... : {str(datetime.datetime.now())}")
    PERI = load_pericentromeres(peri_path, format='bed3')
    print('Length of library before:', len(library_full))
    """removing extra (C, M) chromosomes"""
    #library_full = library_full[library_full['CHROM'].isin([f'Chr{i}' for i in range(1,6)])]
    library_full = library_full[~library_full['CHROM'].isin([f'Chr{i}' for i in ['C', 'M']])]
    """removing pericentromeric regions"""
    library_full = removePeri(library_full, PERI) if library_format=='graffite' else removePeri(library_full, PERI,
                                                                                                chr_col='CHROM',
                                                                                                start_col='start_ref',
                                                                                                end_col='end_ref')
    library_full = library_full.sort_values(by=['CHROM'])
    library_full = library_full.reset_index(drop=True)
    print('Length of library after:', len(library_full))

    print(f"Loading genome sequences... : {str(datetime.datetime.now())}")
    genome_cache = {}
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = {executor.submit(load_one_genome, idx, seq_path, seq_suffix): idx for idx in IDS}
        for future in as_completed(futures):
            idx, genome = future.result()
            genome_cache[idx] = genome

    ###########################################################################

    if batch_last == -1:
        batch_last = math.ceil(len(library_full) / batch_size)

    allResults = pd.DataFrame()

    for i_batch in range(batch_first, batch_last + 1):
        print(f'Batch {i_batch - batch_first + 1} ({i_batch}) out of '
              f'{batch_last - batch_first + 1} : {str(datetime.datetime.now())}')
        
        start = (i_batch - 1) * batch_size
        end   = min(i_batch * batch_size, len(library_full))
        
        library = library_full.iloc[start:end, :]
        if library_format=='graffite':
            library = process_graffite(library, vcf_version, IDS)
        elif library_format=='bed':
            library = process_bed(library, IDS)



        result_batch = BLAST_flanks_and_TE(
                            library,
                            out_path, i_batch,
                            seq_ref, genome_cache, gene_cache,
                            n_workers=n_threads  
                        )
        
        result_batch.to_csv(f'{out_path}/blastedSummary/blastedSummary_batch_{i_batch}.csv')
        allResults = pd.concat([allResults, result_batch])
    
    allResults.to_csv(f'{out_path}/blastedSummary/blastedSummary.csv', index=False)

    allResults_binary = resultsToBinary(allResults, IDS)
    allResults_binary.to_csv(f'{out_path}/blastedSummary/binarySummary.csv')


############################################################################################
############################################################################################
############################################################################################