import pandas as pd
import numpy as np
from io import StringIO


def load_library_graffite(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')
    # Strip the leading '#' from the first column name
    df.columns = [df.columns[0].lstrip('#')] + list(df.columns[1:])
    return df

def graffiteToBed(graffite):
    """
    Converts Graffite to BED
    Returns raw DataFrame with columns ['Chr', 'start_ref', 'end_ref', 'strand_ref', 'ID', 'SEQ', 'FAMS'].
    """
    df = pd.DataFrame(index=graffite.index)
    df["Chr"] = graffite["CHROM"]
    df["start_ref"] = graffite["POS"]
    df["end_ref"] = np.nan
    df["strand_ref"] = np.nan
    df["ID"] = graffite["ID"]

    df["SEQ"] = np.where(graffite["REF"].str.len() > graffite["ALT"].str.len(),
                            graffite["REF"],
                            graffite["ALT"])
    

    for i in graffite.index:
        # Parse INFO field
        for field in graffite.loc[i, "INFO"].split(";"):
            key, *vals = field.split("=")
            if key == "repeat_ids":
                df.at[i, "FAMS"] = vals[0] 
            elif key == "END":
                df.at[i, "end_ref"] = int(vals[0])
            elif key == 'RM_hit_strands':
                hit_strands = set(vals[0].split(','))
                if len(hit_strands) == 1 and '+' in hit_strands:
                    df.at[i, "strand_ref"] = 1
                elif len(hit_strands) == 1 and 'C' in hit_strands:
                    df.at[i, "strand_ref"] = -1

    df = df.dropna(subset=["end_ref"])
    df["end_ref"] = df["end_ref"].astype(int)
    return df


graffite_fixed = load_library_graffite('/cluster/CBIO/data1/eantonenko/STEVE/89GENOMES/GraffiTE.merged.genotypes.fixed.nodup.vcf')
bed_fixed = graffiteToBed(graffite_fixed)
bed_fixed.to_csv("/cluster/CBIO/data1/eantonenko/STEVE/89GENOMES/graffite_fixed.bed", 
                sep="\t", header=False, index=False)


graffite_poly = load_library_graffite('/cluster/CBIO/data1/eantonenko/STEVE/89GENOMES/GraffiTE.merged.genotypes.poly.nodup.vcf')
bed_poly = graffiteToBed(graffite_poly)
bed_poly.to_csv("/cluster/CBIO/data1/eantonenko/STEVE/89GENOMES/graffite_poly.bed", 
                sep="\t", header=False, index=False)