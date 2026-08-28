import argparse
import os

from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import numpy as np

import warnings
warnings.filterwarnings('ignore')


# =========================
# Split TEs per genome
# =========================
def _process_te_single(df, idx, out_path):
    """Process a single genome's TE data. Designed for parallel execution."""

    required = ['Chr', 'ID', 'FAMS', 'SUPERFAMS', 'strand_ref',
                f'{idx}_btwflanks_start', f'{idx}_btwflanks_end',
                f'{idx}_left_start', f'{idx}_left_end',
                f'{idx}_right_start', f'{idx}_right_end']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"[{idx}] Missing columns: {missing}")

    def assign_strand(te_start, te_end, strand_ref):
        """Determine strand orientation from reference strand and TE direction."""
        forward = te_start < te_end
        print(te_start, te_end)
        print(forward)
        if strand_ref == 1:
            return 1 if forward else -1
        if strand_ref == -1:
            return -1 if forward else 1
        return np.nan

    # Normalise each flank to (min, max) regardless of blast direction
    left_min  = df[[f'{idx}_left_start',  f'{idx}_left_end']].min(axis=1)
    left_max  = df[[f'{idx}_left_start',  f'{idx}_left_end']].max(axis=1)
    right_min = df[[f'{idx}_right_start', f'{idx}_right_end']].min(axis=1)
    right_max = df[[f'{idx}_right_start', f'{idx}_right_end']].max(axis=1)

    # 5p flank = whichever has the smaller left coordinate
    # 3p flank = whichever has the bigger left coordinate
    #flank_5p_start = np.where(left_min <= right_min, left_min,  right_min)
    flank_5p_end   = np.where(left_min <= right_min, left_max,  right_max)
    flank_3p_start = np.where(left_min <= right_min, right_min, left_min)
    #flank_3p_end   = np.where(left_min <= right_min, right_max, left_max)

    dfg = pd.DataFrame({
            'Chr':                  df['Chr'].values,
            'start':                df[f'{idx}_btwflanks_start'].values,
            'end':                  df[f'{idx}_btwflanks_end'].values,
            'strand':               np.nan,
            'id.tair10':            df['ID'].values,
            'TE.family':            df['FAMS'].str.strip("[] '").str.split(',').str[0].str.strip("'"),
            'TE.superfamily':       df['SUPERFAMS'].str.strip("[] '").str.split(',').str[0].str.strip("'"),
            'flank_5p_end':         flank_5p_end,
            'flank_3p_start':       flank_3p_start,
        }, index=df.index)

    # --- Strand assignment ---
    strand_ref = df['strand_ref']
    # Convert only numeric values (-1, 1), keep '.' as is
    strand_num = pd.to_numeric(strand_ref, errors='coerce')
    # Apply inversion only to numeric values
    dfg['strand'] = np.where(
        strand_ref == '.',  # keep '.' unchanged
        '.',
        np.where(dfg['start'] < dfg['end'], strand_num, -strand_num)
    )

    # --- Normalise start < end ---
    swap = dfg['end'] < dfg['start']
    dfg.loc[swap, ['start', 'end']] = dfg.loc[swap, ['end', 'start']].values

    # --- Length filters ---
    dfg = dfg[abs(dfg['end']   - dfg['start'])          > 200]
    dfg = dfg[abs(dfg['start'] - dfg['flank_5p_end']) < 50]
    dfg = dfg[abs(dfg['end']   - dfg['flank_3p_start']) < 50]

    dfg[['start', 'end']] = dfg[['start', 'end']].astype(int)
    dfg = dfg.drop(columns=['flank_5p_end', 'flank_3p_start'])
    dfg = dfg.sort_values(by=['Chr', 'start', 'end'])

    out_file = f'{out_path}/TEperGenome/{idx}.bed'
    dfg.to_csv(out_file, header=None, index=None, sep='\t')
    return idx, len(dfg)   # return metadata useful for logging


def splitByGenome_TE(df, IDS, out_path, n_threads=None):
    os.makedirs(f'{out_path}/TEperGenome', exist_ok=True)
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = {
            executor.submit(_process_te_single, df, idx, out_path): idx
            for idx in IDS
        }
        for future in as_completed(futures):
            idx = futures[future]
            try:
                _, n_rows = future.result()
                print(f'[TE]    {idx}: {n_rows} TEs written')
            except Exception as e:
                print(f'[TE]    {idx}: FAILED — {e}')

                                    
# =========================
# Split Flanks per genome
# =========================
def _process_flanks_single(df, idx, out_path):
    """Process a single genome's flank data. Designed for parallel execution."""

    required = [f'{idx}_left_start', f'{idx}_left_end',
                f'{idx}_right_start', f'{idx}_right_end']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"[{idx}] Missing columns: {missing}")

    # Normalise each flank to (min, max) regardless of blast direction
    left_min  = df[[f'{idx}_left_start',  f'{idx}_left_end']].min(axis=1)
    left_max  = df[[f'{idx}_left_start',  f'{idx}_left_end']].max(axis=1)
    right_min = df[[f'{idx}_right_start', f'{idx}_right_end']].min(axis=1)
    right_max = df[[f'{idx}_right_start', f'{idx}_right_end']].max(axis=1)

    # 5p flank = whichever has the smaller left coordinate
    # 3p flank = whichever has the bigger left coordinate
    flank_5p_start = np.where(left_min <= right_min, left_min,  right_min)
    flank_5p_end   = np.where(left_min <= right_min, left_max,  right_max)
    flank_3p_start = np.where(left_min <= right_min, right_min, left_min)
    flank_3p_end   = np.where(left_min <= right_min, right_max, left_max)

    # Span from 5p flank start to 3p flank end (may include overlap)
    span_start = np.minimum(flank_5p_start, flank_3p_start)  # same as flank_5p_start, but safe
    span_end   = np.maximum(flank_5p_end,   flank_3p_end)

    # Gap between flanks: negative = overlap, positive = gap
    gap = flank_3p_start.astype(float) - flank_5p_end.astype(float)

    dfg = pd.DataFrame({
        'Chr':          df['Chr'].values,
        'start':        span_start,
        'end':          span_end,
        'id.tair10':    df['ID'].values,
        'len':          span_end - span_start,
        'gap':          gap,       # negative = overlap, 0 = adjacent, positive = gap
    }, index=df.index)

    # Proximity filter: flanks must nearly touch (gap close to 0), allow overlap
    dfg = dfg[dfg['gap'] < 10]      # remove cases where flanks are far apart
    dfg = dfg[dfg['len'] > 500]

    dfg[['start', 'end', 'len']] = dfg[['start', 'end', 'len']].astype(int)
    dfg = dfg.sort_values(by=['Chr', 'start', 'end'])

    out_file = f'{out_path}/FlanksPerGenome/{idx}.flanks.bed'
    dfg.to_csv(out_file, header=None, index=None, sep='\t')
    return idx, len(dfg)


def splitByGenome_flanks(df, IDS, out_path, n_threads=None):
    os.makedirs(f'{out_path}/FlanksPerGenome', exist_ok=True)
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = {
            executor.submit(_process_flanks_single, df, idx, out_path): idx
            for idx in IDS
        }
        for future in as_completed(futures):
            idx = futures[future]
            try:
                _, n_rows = future.result()
                print(f'[Flanks] {idx}: {n_rows} flanks written')
            except Exception as e:
                print(f'[Flanks] {idx}: FAILED — {e}')



                                    
############################################################################################
####################################  MAIN  ################################################
############################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--blasted_path", required=True, type=str)
    parser.add_argument("--ids_path", required=True, type=str)
    parser.add_argument("--out_path", required=True, type=str)
    parser.add_argument("--returnTE", default=True, type=bool)
    parser.add_argument("--returnFlanks", default=True, type=bool)
    parser.add_argument("--n_threads", default=None, type=int)
    args = parser.parse_args()
    
    IDS = list(pd.read_csv(args.ids_path, header=None)[0])
    #IDS = ['C24']
    df_blasted = pd.read_csv(args.blasted_path, sep=',')#.set_index('Unnamed: 0')

    if args.returnTE:
        splitByGenome_TE(df_blasted, IDS, args.out_path, args.n_threads)
    if args.returnFlanks:
        splitByGenome_flanks(df_blasted, IDS, args.out_path, args.n_threads)
