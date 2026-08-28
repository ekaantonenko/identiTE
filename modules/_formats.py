"""
Format adapters for different input file versions.
Edit this file to support new input formats without touching the main pipeline.
"""

import pandas as pd
import numpy as np

from io import StringIO


# ==============================
# LIBRARY GRAFFITE (VCF) PARSERS
# ==============================

def detect_vcf_version(path):
    """Detect GraffiTE VCF version from header comments."""
    with open(path) as f:
        for line in f:
            if line.startswith('##fileformat='):
                version = line.strip().split('=')[1]
                return 'v4.2' if version >= 'VCFv4.2' else 'v4.1'
            if not line.startswith('#'):
                break
    return 'v4.1'  # default


def load_library_graffite(path, version=None):
    """
    Load GraffiTE VCF, auto-detecting version if not specified.
    Returns raw DataFrame with standard column names.
    """
    if version is None:
        version = detect_vcf_version(path)

    # Skip lines starting with '##', use the '#CHROM' line as header
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')
    # Strip the leading '#' from the first column name
    df.columns = [df.columns[0].lstrip('#')] + list(df.columns[1:])
    return df, version


def read_vcf(filepath):
    
    
    
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')
    
    # Strip the leading '#' from the first column name
    df.columns = [df.columns[0].lstrip('#')] + list(df.columns[1:])
    
    return df


def parse_info_field_v41(info_str):
    """Parse INFO field for GraffiTE v4.1. Returns dict of key->value."""
    result = {}
    for field in info_str.split(';'):
        key, *vals = field.split('=')
        if vals:
            result[key] = vals[0]
    return result


def parse_info_field_v42(info_str):
    """Parse INFO field for GraffiTE v4.2 (adjust if format changed)."""
    return parse_info_field_v41(info_str)  # update if v4.2 differs


def get_info_parser(version):
    """Return the appropriate INFO field parser for a given version."""
    return {
        'v4.1': parse_info_field_v41,
        'v4.2': parse_info_field_v42,
    }.get(version, parse_info_field_v41)


# =========================
# LIBRARY BED PARSERS
# =========================

def load_library_bed(path):
    """
    Load BED
    Returns raw DataFrame with standard column names.
    """

    with open(path, 'r') as f:
        lines = [l for l in f]
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')
    # Strip the leading '#' from the first column name
    df.columns = ['CHROM', 'start_ref', 'end_ref', 'strand_ref', 'ID', 'SEQ', 'FAMS']
    return df

# =========================
# FASTA PARSER
# =========================

def read_fasta(path):
    """
    Read FASTA into {header: sequence} dict.
    Handles both single-line and multi-line (wrapped) sequences.
    Strips '_RagTag' suffix from chromosome names.
    Skips sequences whose names start with 'contig'.
    """
    seqs, header, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None and header is not False:
                    seqs[header] = ''.join(chunks)
                name = line[1:].split()[0].removesuffix('_RagTag')
                header = False if name.startswith('contig') else name
                chunks = []
            else:
                if header is not False:
                    chunks.append(line.upper())
        if header is not None and header is not False:
            seqs[header] = ''.join(chunks)
    return seqs




# =========================
# GENE ANNOTATION PARSERS
# =========================

def load_gene_gff(path, format='liftoff'):
    """
    Load gene annotations from GFF-like file.
    Supported formats: 'liftoff', 'ensembl', 'bed6'
    Returns DataFrame with columns: Chr, start_ref, end_ref, strand, id
    """
    loaders = {
        'liftoff':  _load_gene_gff_liftoff,
        'ensembl':  _load_gene_gff_ensembl,
        'bed10':     _load_gene_bed10,
    }
    if format not in loaders:
        raise ValueError(f"Unknown gene annotation format: '{format}'. "
                         f"Supported: {list(loaders.keys())}")
    return loaders[format](path)


def _load_gene_gff_liftoff(path):
    gff = pd.read_csv(path, sep='\t', skiprows=3, header=None, encoding='cp1252')
    gff = gff[gff[2] == 'gene'][[0, 3, 4, 6, 8]]
    gff.columns = ['Chr', 'start_ref', 'end_ref', 'strand', 'id']
    gff['id'] = gff['id'].str.split(';').map(lambda d: d[0]).str.lstrip('ID=')
    gff[['start_ref', 'end_ref']] = gff[['start_ref', 'end_ref']].astype(int)
    return gff.reset_index(drop=True)


def _load_gene_gff_ensembl(path):
    """Ensembl GFF3 — ID field uses gene: prefix."""
    gff = pd.read_csv(path, sep='\t', comment='#', header=None)
    gff = gff[gff[2] == 'gene'][[0, 3, 4, 6, 8]]
    gff.columns = ['Chr', 'start_ref', 'end_ref', 'strand', 'id']
    gff['id'] = (gff['id'].str.extract(r'ID=([^;]+)')[0]
                           .str.removeprefix('gene:'))
    gff[['start_ref', 'end_ref']] = gff[['start_ref', 'end_ref']].astype(int)
    return gff.reset_index(drop=True)


def _load_gene_bed10(path):
    """BED10 format: chr, start, end, type_n, *, strand, method, type, **, id."""
    bed = pd.read_csv(path, sep='\t', header=None,
                      names=['Chr', 'start_ref', 'end_ref', 'type_n', '*', 'strand', 'method', 'type', '**', 'id'])
    bed['Chr'] = bed['Chr'].str.rstrip('_RagTag')
    bed = bed[bed['type'] == 'gene']
    bed = bed[['Chr', 'start_ref', 'end_ref', 'strand', 'id']]
    bed['id'] = bed['id'].str.split(';').map(lambda d: d[0]).str.lstrip('ID=')
    bed[['start_ref', 'end_ref']] = bed[['start_ref', 'end_ref']].astype(int)
    return bed.reset_index(drop=True)


# =========================
# PERICENTROMERE PARSERS
# =========================

def load_pericentromeres(path, format='bed3'):
    """
    Load pericentromere regions.
    Supported formats: 'bed3' (chr, start, end), 'bed4' (chr, start, end, name)
    """
    if format == 'bed3':
        df = pd.read_csv(path, sep='\t', header=None,
                         names=['Chr', 'start', 'end'])
    elif format == 'bed4':
        df = pd.read_csv(path, sep='\t', header=None,
                         names=['Chr', 'start', 'end', 'name'])
    else:
        raise ValueError(f"Unknown pericentromere format: '{format}'")
    return df


# =========================
# TE FAMILY MAPPING
# =========================
def convert_family_to_superfamily(fam, cache={}):
    """Map family -> superfamily using lookup file, with caching."""
    if not cache:
        df = pd.read_csv(f'./sources/athaliana_fam_superfam.txt',
                        sep="\t", header=None)
        cache.update(dict(zip(df[0], df[1])))
    return cache.get(fam, "Unknown")


# =========================
# GRAFFITE PARSING
# =========================

def process_graffite(graffite, version, IDS):
    """Parse and filter GraffiTE VCF entries."""
    info = pd.DataFrame(index=graffite.index)
    info["Chr"] = graffite["CHROM"]
    info["start_ref"] = graffite["POS"]
    info["end_ref"] = np.nan
    info["strand_ref"] = np.nan
    info["ID"] = graffite["ID"]
    info["SEQ"] = np.where(graffite["REF"].str.len() > graffite["ALT"].str.len(),
                            graffite["REF"],
                            graffite["ALT"])
    info["FAMS"] = 'Unknown'
    info["SUPERFAMS"] = 'Unknown'

    for i in graffite.index:
        reftype = None   # reset each row
        supp_vec = None
        length = None
        # Parse INFO field
        for field in graffite.loc[i, "INFO"].split(";"):
            key, *vals = field.split("=")
            if key == "repeat_ids":
                info.at[i, "FAMS"] = vals[0]
            elif key == "SVLEN":
                info.at[i, "LEN"] = length = abs(int(vals[0]))
            elif key == "END":
                info.at[i, "end_ref"] = int(vals[0])
            elif key == 'SUPP_VEC':
                supp_vec = vals[0]
            elif key == 'SVTYPE':
                reftype = vals[0]
            elif key == 'RM_hit_strands':
                hit_strands = set(vals[0].split(','))
                if len(hit_strands) == 1 and '+' in hit_strands:
                    info.at[i, "strand_ref"] = 1
                elif len(hit_strands) == 1 and 'C' in hit_strands:
                    info.at[i, "strand_ref"] = -1
        
        if length is None:
            info.at[i, "LEN"] = len(info.at[i, "SEQ"])
        #if strand_ref is None:
        #    info.at[i, "strand_ref"] = 1 if info.at[i, "end_ref"] >= info.at[i, "start_ref"] else -1


        # Map to superfamilies
        info.at[i, "SUPERFAMS"] =  ",".join(dict.fromkeys(
                convert_family_to_superfamily(fam.strip())
                for fam in info.at[i, "FAMS"].split(",")
            ))
        

        # Genotyping 
        if version=='v4.2':
            if reftype is None:
                reftype = 'DEL'
            for g in IDS:   
                gt = graffite.at[i, g].split(";")[0].split(":")[0]
                info.at[i, g] = genotype_call(gt, reftype)


        elif version=='v4.1':
            if reftype is not None and supp_vec is not None:
                for i_g, g in enumerate(IDS):
                    if reftype=='INS':
                        #print(int(supp_vec[i_g]))
                        info.at[i,g] = int(supp_vec[i_g])
                    elif reftype=='DEL':
                        info.at[i,g] = 1-int(supp_vec[i_g])


    info = info[~graffite["REF"].str.contains(",")]
    info = info[~graffite["ALT"].str.contains(",")]
    info = info.replace({"di": np.nan, "?": np.nan}).dropna(subset=IDS)

    return info

# =========================
# BED PARSING
# =========================

def process_bed(library, IDS):
    """ Process .bed entries."""
    info = library.copy()
    info["SUPERFAMS"] = info["FAMS"].map(
        lambda fams: ",".join(dict.fromkeys(
            convert_family_to_superfamily(fam.strip())
            for fam in fams.split(",")
        ))
    )
    info["LEN"] = library["SEQ"].str.len()
    info = info.rename(columns={"CHROM": "Chr"})

    for g in IDS:
        info[g] = np.nan

    return info


# =========================
# GENOTYPE CALL (for genotyped graffite files)
# =========================
def genotype_call(gt, reftype):
    """Interpret genotype string into categorical calls."""
    if gt == "0/0" and reftype == "INS":
        return "0"
    if gt == "0/0" and reftype == "DEL":
        return "1"
    if gt == "1/1" and reftype == "INS":
        return "1"
    if gt == "1/1" and reftype == "DEL":
        return "0"
    if gt in {"0/1", "1/0"}:
        return "di"
    return "?"