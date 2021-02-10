#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import argparse
import sys
from functools import lru_cache

import pandas as pd

from gff3 import Gff3
from tqdm import tqdm


class Modifier:
    """
    A class that extracts important information from a annotation file to
    modify a gff file.
    """

    def __init__(self, anno_file: list):
        """
        Creates a new instance of a gff modifier.

        Parameters:
            anno_path:
                A path to the annotation file.

            anno_delimiter:
                The delimiter for the annotation file. The default is a
                single tab.
        """

        self.__anno_dfs = []

        self._num_df_dots = None

        for path, ID_index, ref_index in anno_file:

            ID_index = int(ID_index)
            ref_index = int(ref_index)

            open_anno_kwargs = {"anno_path": path,
                                "ID_index": ID_index, "ref_index": ref_index}

            anno_df = self.open_anno_file(**open_anno_kwargs)

            anno_df = anno_df[~anno_df.index.duplicated(keep='first')]

            self.__anno_dfs.append(anno_df)

            self._num_df_dots = len(
                list(filter(lambda c: c == '.', anno_df.index[0]))) + 1

    @lru_cache()
    def __getitem__(self, index: str) -> str:
        """
        Returns the corresponding extended information for the given index.
        """

        return_list = []

        for df in self.__anno_dfs:

            try:
                anno_row = df.loc[index]

            except KeyError:

                if index.startswith('cds.'):
                    index = '.'.join(index.split('.')[1:self._num_df_dots + 1])
                else:
                    index = '.'.join(index.split('.')[0:self._num_df_dots])

                try:
                    # Get the corresponding row in the annotation file
                    anno_row = df.loc[index]
                except KeyError:
                    continue

            accession = self.extract_value(anno_row, "ref")

            if accession is not None:
                accession = self.process_accession(accession)
                return_list.append(("Dbxref", accession))

        return return_list

    def list_to_dict(self, input_list: list):

        new_dict = dict(zip((k for k, _ in input_list), [
                        list() for _ in range(len(input_list))]))

        for k, v in input_list:

            new_dict[k].append(v)

        return new_dict

    def modify_gff(self, gff: Gff3):
        """
        Modifies an existing Gff3 object by adding contents from the input
        annotation file.
        """

        for line in tqdm(iterable=gff.lines, desc='Modify Compilation', ascii=True):

            # Get the ID to use to index on this modifier
            gene_ID = line['attributes']['ID']
            update_list = self[gene_ID]

            line['attributes'].update(self.list_to_dict(update_list))

        return

    def open_anno_file(self, anno_path: str = None, ID_index: int = 0, ref_index: int = 1):
        """
        Opens the annotations file.

        Parameters:
            anno_path:
                The path to the annotation file.

            ID_index:
                The ID index to the corresponding annotation file.

            ref_index:
                The reference index to the corresponding annotation file.

        Return:
            Returns a 2-column pandas dataframe with the first column being
            the gene ID and the second column being the reference index.
        """

        anno_df = pd.read_csv(anno_path, engine='python',
                              sep='\t', header=None, usecols=[ID_index, ref_index], names=["ID", "ref"], index_col="ID").astype(str)

        return anno_df

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
            extracted_value: str = anno_row[value]
        except KeyError:
            return None

        # Don't add the gene name if it has na
        if extracted_value.lower() in ['na', 'nan']:
            return None

        # Check for error by searching for a semi-colon
        if ";" in extracted_value:
            extracted_value, _ = extracted_value.split(";", maxsplit=1)
            extracted_value = extracted_value.strip()

        return extracted_value

    def process_accession(self, accession: str) -> str:
        """
        Processes the accession value from the annotation file.
        """

        if accession.startswith("tr|") or accession.startswith("sp|"):

            db_xref_prefix, db_xref_id, *_ = accession.split(sep='|')

            db_xref_prefix, db_xref_id = db_xref_prefix.strip(), db_xref_id.strip()

            if db_xref_prefix == "sp":
                return 'UniProtKB/Swiss-Prot:{0}'.format(db_xref_id)
            elif db_xref_prefix == "tr":
                return 'UniProtKB/TrEMBL:{0}'.format(db_xref_id)

        elif accession.startswith("PF"):
            return 'PFAM:{0}'.format(accession)

        raise NotImplementedError(
            "Not equipped to handle db xref (" + accession + ")")


def run_modifier(args):

    modifier = Modifier(args.annotation)
    print("Reading gff file")
    gff: Gff3 = Gff3(gff_file=args.gff_path)

    # Modify the gff file using the Modifier class
    modifier.modify_gff(gff)

    print("Writing modified gff file")
    # Write the modified gff to the output path
    if args.output_path is None:
        gff.write(sys.stdout)

    else:
        with open(args.output_path, "w") as file_out:
            gff.write(file_out)


def main():

    # Pgla_CCMP1383 usage: (TODO: update paths)
    #   python .\modmygff.py --gff_path ".\data\Pgla_CCMP1383\Polarella_glacialis_CCMP1383_PredGenes_v1.gff3" --annotation ".\data\Pgla_CCMP1383\CCMP1383_UniProt.tsv" 0 1  --annotation ".\data\Pgla_CCMP1383\CCMP1383_scaffolds_PFAM.tsv" 0 5 --output_path ".\data\Polarella_glacialis_CCMP1383_PredGenes_v1_ext.gff3"

    # Pgla_CCMP2088 usage:
    #   python .\modmygff.py --gff_path ".\data\Pgla_CCMP2088\Polarella_glacialis_CCMP2088.gff3" --annotation ".\data\Pgla_CCMP2088\CCMP2088_UniProt.tsv" 0 1  --annotation ".\data\Pgla_CCMP2088\CCMP2088_pfam.tsv" 0 5 --output_path ".\data\Pgla_CCMP2088\Polarella_glacialis_CCMP2088_ext.gff3"

    parser = argparse.ArgumentParser(description="Creates a flat file from a "
                                     "given gff file and annotations file.")

    parser.add_argument('--gff_path', type=str, required=True,
                        help='A file path to the gff file.')
    parser.add_argument('--annotation', action='append', nargs=3, required=True,
                        metavar=('path', 'index_col', 'ref_col'),
                        help="A file path to the annotation file.")

    parser.add_argument('--output_path', type=str, required=False, default=None,
                        help='A file path to output the contents of the flatfile. '
                        'Default output file is stdout.')

    args = parser.parse_args()
    run_modifier(args)

    exit(0)


if __name__ == '__main__':
    main()
