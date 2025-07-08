#!/usr/bin/env python3
"""
get_conserved_GRs.py

Get a GFF with only conserved GRs and rename features for merging with Genewise.

Usage:
    python get_conserved_GRs.py ABCENTH_GRs.txt ABCENTH_clean_GR_renamed_all_nofragment.gff3 

Outputs:
    {prefix}_conserved.gff3
    {prefix}_peptides.fasta
"""
import sys
import re

def get_prefix(intable):
    m = re.match(r'(\S+)_GRs.txt', intable)
    if m:
        return m.group(1)
    else:
        sys.exit(f"Can't find prefix in {intable} file")

def get_conserved_genes(intable):
    conserved = set()
    with open(intable) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) > 9 and fields[9] != 'None':
                conserved.add(fields[4])
    return conserved

def process_gff(gff, conserved, outgff, outgffall):
    with open(gff) as fin, \
         open(outgff, 'w') as fout, \
         open(outgffall, 'w') as foutall:
        fout.write('##gff-version 3\n')
        foutall.write('##gff-version 3\n')
        for line in fin:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            feature = fields[2]
            attrs = fields[8]
            if feature == 'CDS':
                m = re.search(r'Parent=([^;]+)(\S+)', attrs)
                if not m:
                    sys.exit(f"ERROR: It fails detecting Parent ID in {line}")
                genename = m.group(1)
                rest = m.group(2)
                nnamef = f"Abcenth_{genename}"
                new_attrs = f"Parent={nnamef}{rest};"
                out_line = '\t'.join(fields[:8] + [new_attrs]) + '\n'
                foutall.write(out_line)
                if genename in conserved:
                    fout.write(out_line)
            elif feature == 'mRNA':
                m = re.search(r'ID=([^;]+);Parent=[^;]+(;\S+)', attrs)
                if not m:
                    print(f"ERROR: It fails detecting ID in {line}")
                    continue
                genename = m.group(1)
                rest = m.group(2)
                nnamef = f"Abcenth_{genename}"
                gene_line = '\t'.join(fields[:2] + ['gene'] + fields[3:8] + [f"ID=g{nnamef}{rest};"]) + '\n'
                mrna_line = '\t'.join(fields[:8] + [f"ID={nnamef};Parent=g{nnamef}{rest};"]) + '\n'
                foutall.write(gene_line)
                foutall.write(mrna_line)
                if genename in conserved:
                    fout.write(gene_line)
                    fout.write(mrna_line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <ABCENTH.txt> <input.gff3>")
        sys.exit(1)
    intable, gff = sys.argv[1:3]
    prefix = get_prefix(intable)
    outgffall = f"{prefix}_GR_all_renamed.gff3"
    outgff = f"{prefix}_GR_Aglyconserved_renamed.gff3"
    conserved = get_conserved_genes(intable)
    process_gff(gff, conserved, outgff, outgffall)

