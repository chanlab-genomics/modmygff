#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import argparse
import os
import sys
import warnings
from functools import lru_cache
from pprint import pprint
from typing import (Any, Callable, Dict, Iterable, List, Optional, Tuple, Type,
                    Union)

import pandas as pd

from gff3 import Gff3

"""
Example:
    (Windows)
    python -m modmygff --gff_path .\data\Slin_CCMP2456\S.linucheae_CCMP2456_eg1.gff --anno_path .\data\Slin_CCMP2456\S.linucheae_CCMP2456_uniprot_annotated.tsv --output_path .\data\Slin_CCMP2456\S.linucheae_CCMP2456_eg1_ext.gff
"""


def test_stuff():
    eg_gff_path = os.path.join(os.getcwd(), 'data', 'Slin_CCMP2456',
                               'S.linucheae_CCMP2456_eg1.gff')

    eg_gff: Gff3 = Gff3(gff_file=eg_gff_path)
    pprint(eg_gff.lines[0].keys())
    pprint(eg_gff.lines[0]['attributes'])
    pprint(eg_gff.lines[0]['attributes']['ID'])
    eg_gff.lines[0]['attributes']['Dbxref'] = 'UniProtKB/Swiss-Prot:P12345'
    eg_gff.lines[0]['attributes']['geneID'] = 'amtB'

    # eg_gff_ext_path = os.path.join(os.getcwd(), 'data', 'Slin_CCMP2456',
    #                                'S.linucheae_CCMP2456_eg1_ext.gff')

    # with open(eg_gff_ext_path, 'w') as eg_gff_ext_file:
    #     eg_gff.write(eg_gff_ext_file)


class Modifier:
    """
    A class that extracts important information from a annotation file to
    modify a gff file.
    """

    def __init__(self, anno_path: str, anno_delim: str = '\t'):
        """
        Creates a new instance of a gff modifier.

        Parameters:
            anno_path:
                A path to the annotation file.

            anno_delimiter:
                The delimiter for the annotation file. The default is a
                single tab.
        """

        if anno_delim == 'None':
            anno_delim = None

        else:
            for old, new in [('\\n', '\n'), ('\\t', '\t'), ('\\r', '\r')]:
                anno_delim = anno_delim.replace(old, new)

        self._anno_path: str = anno_path
        self._anno_delim: str = anno_delim

        # Open the annotations path as a pandas dataframe.
        # Set sep=None so that the Python parsing engine can automatically
        # detect the separator.
        self.__anno_df = self.open_anno_file()

    @lru_cache()
    def __getitem__(self, index: str) -> str:
        """
        Returns the corresponding extended information for the given index.
        """

        if 'mRNA1'.lower() not in index.lower():
            index += '.mRNA1'
        elif index.lower().endswith('.mrna1'):
            pass
        else:
            index = '.'.join(index.split('.')[:3])

        try:
            # Get the corresponding row in the annotation file
            anno_row = self.__anno_df.loc[index]
        except KeyError:
            return {}

        return_dict = {}

        gene_name = self.extract_value(anno_row, "gene_name")

        if gene_name is not None:
            return_dict["geneID"] = gene_name

        accession = self.extract_value(anno_row, "accession")

        if accession is not None:
            accession = self.process_accession(accession)
            return_dict["Dbxref"] = accession

        return return_dict

    def modify_gff(self, gff: Gff3):
        """
        Modifies an existing Gff3 object by adding contents from the input
        annotation file.
        """

        for line in gff.lines:

            # Get the ID to use to index on this modifier
            gene_ID = line['attributes']['ID']
            line['attributes'].update(self[gene_ID])

        return

    def open_anno_file(self):
        """
        Opens the annotations file.
        """

        # Try opening the annotations file using the given delimiter using the
        # given delimiter using the c-engine. We will only use the python engine
        # if specified the delimiter is None.

        if self._anno_delim is None:
            return pd.read_csv(
                self._anno_path, engine='python', sep=None,  index_col='sequence').astype(str)

        anno_df = pd.read_csv(
            self._anno_path, engine='c', sep=self._anno_delim,  index_col='sequence')

        _, df_cols = anno_df.shape

        if df_cols != 0:
            return anno_df.astype(str)

        # If we only have one column then the delimiter we are using probably is
        # not correct. We should have at least two columns.

        warning_message = 'Could not separate annotation dataframe' + \
            'with --delimiter={0!r}. Switching to python engine parser.'
        warning_message = warning_message.format(self._anno_delim)

        warnings.warn(warning_message, RuntimeWarning)

        return pd.read_csv(
            self._anno_path, engine='python', sep=None,  index_col='sequence').astype(str)

    def extract_value(self, anno_row: pd.Series, value: str) -> str:
        """
        Extracts a certain value from a certain row of the annotation file.

        Parameters:
            anno_row
                The row of the annotation file to extract the value.

            value:
                The value to be extracted from the row.

            add_to_dict:
                If True, the extracted value is added to self.__kwargs without
                any further processing.

        Returns:
            Returns the extracted value as a string. None otherwise.
        """

        try:
            extracted_value = anno_row[value]
        except KeyError:
            return None

        # Don't add the gene name if it has na
        if extracted_value.lower() in ['na', 'nan']:
            return None

        return extracted_value

    def process_accession(self, accession: str) -> str:
        """
        Processes the accession value from the annotation file.
        """

        db_xref_prefix, db_xref_id, *_ = accession.split(sep='|')

        db_xref_prefix, db_xref_id = db_xref_prefix.strip(), db_xref_id.strip()

        if db_xref_prefix == "sp":
            return 'db_xref="UniProtKB/Swiss-Prot:{0}"'.format(db_xref_id)
        elif db_xref_prefix == "tr":
            return 'db_xref="UniProtKB/TrEMBL:{0}"'.format(db_xref_id)
        else:
            raise NotImplementedError(
                "Not equipped to handle db xref (" + db_xref_id + ") with prefix: " + repr(db_xref_prefix))


def main():
    eg_anno_path = os.path.join(os.getcwd(), 'data', 'Slin_CCMP2456',
                                'S.linucheae_CCMP2456_uniprot_annotated.tsv')

    test_mod = Modifier(eg_anno_path)

    eg_gff_path = os.path.join(os.getcwd(), 'data', 'Slin_CCMP2456',
                               'S.linucheae_CCMP2456_eg1.gff')

    eg_gff: Gff3 = Gff3(gff_file=eg_gff_path)
    print(eg_gff.lines[0]["attributes"])
    test_mod.modify_gff(eg_gff)
    print(eg_gff.lines[0]["attributes"])


def run_modifier(args):

    modifier = Modifier(args.anno_path, anno_delim=args.anno_delim)
    gff: Gff3 = Gff3(gff_file=args.gff_path)

    # Modify the gff file using the Modifier class
    modifier.modify_gff(gff)

    # Write the modified gff to the output path
    if args.output_path is None:
        gff.write(sys.stdout)

    else:
        with open(args.output_path, "w") as file_out:
            gff.write(file_out)


def main():

    parser = argparse.ArgumentParser(description="Creates a flat file from a "
                                     "given gff file and annotations file.")

    parser.add_argument('--gff_path', type=str, required=True,
                        help='A file path to the gff file.')
    parser.add_argument('--anno_path', type=str, required=True,
                        help="A file path to the annotation file. "
                        "Every annotation MUST have a column named 'sequence' "
                        "which is a unique identifier. Other optional columns "
                        "to include are: 'accession' 'function' 'gene_name' "
                        "'gene_synonym' and 'EC_number'. NOTE: these column "
                        "names are case sensitive.")

    parser.add_argument('--output_path', type=str, required=False, default=None,
                        help='A file path to output the contents of the flatfile. '
                        'Default output file is stdout.')
    parser.add_argument('--anno_delim', type=str,
                        required=False, default='\t',
                        help='A delimiter value for the annotation file.')

    args = parser.parse_args()
    run_modifier(args)

    exit(0)


if __name__ == '__main__':
    main()
