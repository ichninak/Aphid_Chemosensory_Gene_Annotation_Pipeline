#!/usr/bin/env python3
"""
get_conserved_GRs.py

Get a GFF with only conserved GRs and rename features for merging with Genewise.

Usage:
    python get_conserved_GRs.py ABCENTH_GRs.txt ABCENTH_clean_GR_renamed_all_nofragment.gff3 \
        --genome genome.fasta

Outputs:
    {prefix}_conserved.gff3
    {prefix}_peptides.fasta
"""
import argparse
import os
import subprocess
import gffutils

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_table", help="Table listing conserved GR IDs (e.g. ABCENTH_GRs.txt)")
    parser.add_argument("in_gff", help="Input GFF3 file (e.g. ABCENTH_clean_GR_renamed_all_nofragment.gff3)")
    parser.add_argument("--genome", required=True, help="Genome FASTA for sequence extraction")
    parser.add_argument("--prefix", help="Custom prefix for outputs (defaults to in_table basename)")
    return parser.parse_args()

def main():
    args = parse_args()
    # Determine prefix from in_table if not provided
    prefix = args.prefix or os.path.splitext(os.path.basename(args.in_table))[0]

    # Read conserved IDs from table
    conserved_ids = set()
    with open(args.in_table) as tbl:
        for line in tbl:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            conserved_ids.add(line.split()[0])

    # Build an in-memory GFF database
    db = gffutils.create_db(
        data=args.in_gff,
        dbfn=':memory:',
        force=True,
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True
    )

    # Filter and rename features
    out_gff = f"{prefix}_conserved.gff3"
    with open(out_gff, 'w') as out:
        out.write("##gff-version 3\n")
        for feature in db.all_features(order_by=('seqid', 'start')):
            feat_id = feature.attributes.get('ID', [None])[0]
            parents = feature.attributes.get('Parent', [])
            if feat_id in conserved_ids or any(p in conserved_ids for p in parents):
                # Rename ID and Parent attributes by adding prefix
                attrs = feature.attributes.copy()
                if 'ID' in attrs:
                    attrs['ID'] = [f"{prefix}_{attrs['ID'][0]}"]
                if 'Parent' in attrs:
                    attrs['Parent'] = [f"{prefix}_{p}" for p in attrs['Parent']]
                # Build attribute string
                attr_str = ';'.join(f"{k}={','.join(v)}" for k, v in attrs.items())
                # Write GFF line
                out.write(
                    '\t'.join([
                        feature.seqid,
                        feature.source,
                        feature.featuretype,
                        str(feature.start),
                        str(feature.end),
                        feature.score or '.',
                        feature.strand,
                        feature.frame or '.',
                        attr_str
                    ]) + '\n'
                )

    # Extract peptide sequences using external gff2fasta_v3.pl script
    prot_fasta = f"{prefix}_peptides.fasta"
    cmd = [
        'perl',
        'gff2fasta_v3.pl',
        args.genome,
        out_gff,
        prefix
    ]
    subprocess.run(cmd, check=True)

    # Clean up trailing X residues (like sed 's/X*$//')
    tmp = prot_fasta + '.tmp'
    with open(prot_fasta) as inp, open(tmp, 'w') as outp:
        for line in inp:
            if line.startswith('>'):
                outp.write(line)
            else:
                # Remove trailing X
                seq = line.rstrip('\n').rstrip('X')
                outp.write(seq + '\n')
    os.replace(tmp, prot_fasta)

if __name__ == '__main__':
    main()
