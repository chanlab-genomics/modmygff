#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import argparse
import pandas as pd


def open_annotation(anno_path: str):

    anno_df = pd.read_csv(anno_path, engine='python',
                          sep='\t', header=None, usecols=[0, 1], names=["ID", "ref"], index_col="ID").astype(str)

    # Extract the first two columns
    print(anno_df)

    return anno_df


"""
test_anno_path = r".\data\Pgla_CCMP1383\CCMP1383_scaffolds.fa.GenePred_v7.evm.filtered.pep.fa.blastp_UniProt_Combined.evalue-5.TopHit.outfmt6"

parser = argparse.ArgumentParser(description="Test.")
parser.add_argument('--annotation', action='append', nargs=3, required=True,
                    metavar=('path', 'index_col', 'ref_col'))
args = parser.parse_args()

print(args.annotation)
"""

test_anno_path = r".\data\Pgla_CCMP1383\CCMP1383_scaffolds.txt"

anno_df = pd.read_csv(test_anno_path, engine='python',
                      sep="\t", header=None).astype(str)
print(anno_df.shape)
