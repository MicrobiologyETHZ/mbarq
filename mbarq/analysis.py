import logging
import collections
import os
from typing import Optional, List, Union
from mbarq.core import Barcode, BarSeqData
from pathlib import Path
import pandas as pd
import sys
from mbarq.mbarq_logger import get_logger

from pandera import Column, DataFrameSchema, Check, Index, SchemaError


class CountDataSet:

    """
    Merge count files into one or if given one count file validate proper format
    """
    def __init__(self, count_files: Union[str, List[str]], name: str,
                 gene_column_name: str = '', sep=',', output_dir: str = ''):
        self.count_files = count_files
        self.gene_name = gene_column_name
        self.sep = sep
        self.countTable, self.sampleIDs = self.createCountTable()
        self.output_dir = Path(output_dir)
        self.name = name
        self.logger = get_logger('merging-log', self.output_dir / f"{self.name}_CountDataSet.log")

    def _merge_count_files(self) -> pd.DataFrame:
        df_list = []
        sampleIDs = []
        for count_file in self.count_files:
            sampleID = Path(count_file).stem.split("_mbarq")[0]
            df = pd.read_table(Path(count_file), sep=self.sep, index_col=0)  # todo don't write down index in count step
            if not self.gene_name:
                df = df.iloc[:, [0, 1]]
                df.columns = ['barcode', 'barcode_count']
                index_cols = ['barcode']
            elif self.gene_name and self.gene_name in df.columns:
                gene_col_index = df.columns.get_loc(self.gene_name)
                df = df.iloc[:, [0, 1, gene_col_index]]
                df.columns = ['barcode', 'barcode_count', self.gene_name]
                index_cols = ['barcode', self.gene_name]
            else:
                self.logger.error(f'Gene name column provided is not found in the {Path(count_file).name}')
                sys.exit(1)
            df = df.assign(sampleID=sampleID).drop_duplicates()  # todo fix bug in counting that produces duplicates
            df_list.append(df)
            sampleIDs.append(sampleID)
        fdf = pd.concat(df_list)
        fdf = fdf.pivot(index=index_cols, columns='sampleID', values=fdf.columns[1]).fillna(0).reset_index()
        if self.gene_name: # todo temp solution, should be addressed at mapping stage
            fdf[self.gene_name] = fdf[self.gene_name.fillna(fdf['barcode'])]
        return fdf

    def _validate_count_table(self, df_to_validate: pd.DataFrame) -> (pd.DataFrame, List[str]):
        # first needs to be a string, second also string, nulls should not be allowed
        # rest should be floats or ints
        num_of_columns = len(df_to_validate.columns)
        num_of_samples = num_of_columns-2 if self.gene_name else num_of_columns-1
        annotation_validations = [Column(str), Column(str)]
        sample_validations = [Column(float)]*num_of_samples
        col_validations = annotation_validations + sample_validations
        col_validation_dict = {col_name: col_validation for col_name, col_validation in
                               zip(df_to_validate.columns, col_validations)}
        schema = DataFrameSchema(
            col_validation_dict,
            strict=True,
            coerce=False,
        )
        try:
            schema.validate(df_to_validate)
            return df_to_validate, df_to_validate.columns[-num_of_samples:]
        except SchemaError:
            logging.error('Invalid count data frame')
            sys.exit(1)

    def createCountTable(self):
        if type(self.count_files) == str:
            df_to_validate = pd.read_table(self.count_files, sep=self.sep)
            return self._validate_count_table(df_to_validate)
        else:
            countTable = self._merge_count_files()
            return self._validate_count_table(countTable)


class DesignMatrix:
    def __init__(self, count_files: List[str], design_table: str = '',
                 control: str = '', treatment: str = '',
                 batch_name: str = '', treatment_name: str = '',
                 cntrl_barcode_file: str = '',
                 ):
        self.count_files = [Path(count_file) for count_file in count_files]
        self.design_table = design_table
        self.cntrl_barcode_file = cntrl_barcode_file

    def processDesignTable(self):
        pass


class Experiment:
    def __init__(self, countDataSet: CountDataSet, designMatrix: DesignMatrix,
                 batch_name: str = '', treatment_name: str = '',
                 control: str = '', treatment: str = '',
                 cntrl_barcode_file: str = '',
                 ):
        pass

    def runBatchCorrection(self):
        pass

    def checkForBottlenecks(self):
        pass