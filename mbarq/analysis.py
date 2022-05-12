import logging
import collections
import os
import pandera.errors
from typing import Optional, List, Union
from mbarq.core import Barcode, BarSeqData
import numpy as np
from pathlib import Path
import pandas as pd
import sys
from mbarq.mbarq_logger import get_logger
from pandera import Column, DataFrameSchema, Check, Index


class CountDataSet:

    """
    Merge count files into one or if given one count file validate proper format

    """
    def __init__(self, count_files: Union[str, List[str]], name: str,
                 gene_column_name: str = '', sep=',', output_dir: str = ''):
        self.output_dir = Path(output_dir)
        self.name = name
        self.logger = get_logger('merging-log', self.output_dir / f"{self.name}_CountDataSet.log")
        self.count_files = count_files
        self.gene_name: str = gene_column_name
        self.sep: str = sep
        self.count_table: pd.DataFrame = pd.DataFrame()
        self.cpms: pd.DataFrame = pd.DataFrame()
        self.sampleIDs: List[str] = []

    def _merge_count_files(self) -> pd.DataFrame:
        df_list = []
        sampleIDs = []
        index_cols = []
        for count_file in self.count_files:
            sampleID = Path(count_file).stem.split("_mbarq")[0]
            df = pd.read_table(Path(count_file), sep=self.sep)
            if not self.gene_name:
                df = df[['barcode', 'barcode_count']]
                index_cols = ['barcode']
            elif self.gene_name and self.gene_name in df.columns:
                df = df[['barcode', 'barcode_count', self.gene_name]]
                index_cols = ['barcode', self.gene_name]
            else:
                self.logger.warning(f'Gene name column provided is not found in the {Path(count_file).name}. Skipping file')
                continue
            df = df.assign(sampleID=sampleID).drop_duplicates()  # todo fix bug in counting that produces duplicates
            df_list.append(df)
            sampleIDs.append(sampleID)
        if len(df_list) == 0 | len(index_cols) == 0:
            self.logger.error(f'Could not process any of the files. Check that {self.gene_name} is correct')
            sys.exit(1)
        fdf = pd.concat(df_list)
        fdf = (fdf.pivot(index=index_cols, columns='sampleID', values=fdf.columns[1])
               .fillna(0)
               .reset_index())
        return fdf

    def _validate_count_table(self, df_to_validate: pd.DataFrame) -> (pd.DataFrame, List[str]):
        # first needs to be a string, second also string if self.gene_name is provided, nulls should not be allowed
        # rest should be floats or ints
        sample_cols = df_to_validate.columns[2:] if self.gene_name else df_to_validate.columns[1:]
        barcode_val = {df_to_validate.columns[0]: Column(str)}
        gene_val = {df_to_validate.columns[1]: Column(str, nullable=True)} if self.gene_name else {}
        sample_val = {sample_col: Column(float) for sample_col in sample_cols}
        validation_dict = {**barcode_val, **gene_val, **sample_val}
        schema = DataFrameSchema(
            validation_dict,
            strict=True,
            coerce=False,
        )
        try:
            schema.validate(df_to_validate)
            return df_to_validate, sample_cols
        except pandera.errors.SchemaError:
            self.logger.error('Invalid count data frame')
            sys.exit(1)

    def create_count_table(self):
        if type(self.count_files) == list:
            df_to_validate = self._merge_count_files()
        else:
            df_to_validate = pd.read_table(self.count_files, sep=self.sep)
        self.count_table, self.sampleIDs = self._validate_count_table(df_to_validate)
        self.count_table.to_csv(self.output_dir/f"{self.name}_mbarq_merged_counts.csv", index=False)

    def calculate_cpms(self):
        self.cpms = self.count_table[self.count_table.sum(axis=1, numeric_only=True) > 10]
        # Normalized for library depth and log transform
        self.cpms = self.cpms.apply(lambda x: np.log2(x/ x.sum() * 1000000 + 0.5) if np.issubdtype(x.dtype, np.number) else x)



class ControlBarcodes:
    """
            Control file: No header

            [barcode],[conc]*,[phenotype]*
            * second and third columns are optional
            1. Check that the first column contains barcodes
            2. Check that second column is numeric -> should contain concentrations
            3. If there is a thrid column it should be string, and at least one should be ['wt', 'WT', 'wildtype']

            """

    def __init__(self, barcode_csv: Union[Path, str] = '', output_dir: Union[Path, str]='.'):
        self.barcode_csv = barcode_csv
        self.wt_barcodes: List = []
        self.wt_bc_conc: pd.DataFrame = pd.DataFrame()
        self.control_bc_conc: pd.DataFrame = pd.DataFrame()
        self.output_dir: Path = Path(output_dir)
        self.logger = get_logger('control-barcode-log', self.output_dir / f"ControlBarcode.log")

    def read_control_file(self):
        if self.barcode_csv:
            cntrl_df = pd.read_csv(self.barcode_csv, header=None)
            num_cols = cntrl_df.shape[1]
            # Add column validation code here

            col_names = ['barcode', 'concentration', 'genotype']

            schema_dict = {"barcode": Column(str),
                          "concentration": Column(float),
                           "genotype": Column(str, [
                        Check(lambda s: any(s.isin(['wt', 'WT', 'wildtype']))),
                    ])}
            cntrl_df.columns = col_names[0:num_cols]
            control_schema = DataFrameSchema({c:schema_dict[c] for c in list(cntrl_df.columns)},
                                             strict=True,
                                             coerce=False,
                                            )
            try:
                control_schema.validate(cntrl_df)
                if 'genotype' in cntrl_df.columns:
                    self.wt_bc_conc = cntrl_df[cntrl_df.genotype.isin(['wt', 'WT', 'wildtype'])].copy()
                    self.control_bc_conc = cntrl_df.copy()
                else:
                    self.wt_bc_conc = cntrl_df.copy()
                self.wt_barcodes = self.wt_bc_conc.barcode.values

            except pandera.errors.SchemaError:
                self.logger.error('Invalid control barcodes file format')
                sys.exit(1)


# class DesignMatrix:
#     def __init__(self, count_files: List[str], design_table: str = '',
#                  control: str = '', treatment: str = '',
#                  batch_name: str = '', treatment_name: str = '',
#                  cntrl_barcode_file: str = '',
#                  ):
#         self.count_files = [Path(count_file) for count_file in count_files]
#         self.design_table = design_table
#         self.cntrl_barcode_file = cntrl_barcode_file
#
#     def processDesignTable(self):
#         pass
#
#     def read_in_sample_data(sample_data_file, sampleIDs, treatment_col="", batch_col=""):
#         """
#         add data validation code
#
#         """
#         return pd.read_csv(sample_data_file)
#
class Experiment:
    def __init__(self, countDataSet: CountDataSet, controlBarcodes: ControlBarcodes,
                 batch_name: str = '', treatment_name: str = '',
                 control: str = '', treatment: str = '',
                 cntrl_barcode_file: str = '',
                 ):
        self.countDataSet = countDataSet
        self.controlBarcodes = controlBarcodes
        self.wt_bc_conc_counts = pd.DataFrame()
        self.control_bc_counts = pd.DataFrame()

    def _get_wt_bc_counts(self):
        """
        """
        if not self.controlBarcodes.wt_bc_conc.empty:
            self.wt_bc_conc_counts = (self.controlBarcodes.wt_bc_conc.merge(self.countDataSet.cpms,
                                                                            how='left', on='barcode')
                                      .fillna(0))

#     def runBatchCorrection(self):
#         pass
#
#     def checkForBottlenecks(self):
#         pass
#

#
#     def calculate_correlation(self, control_df: pd.DataFrame, sampleIDs: List, cutoff: float = 0.8):
#         """
#         Given a data frame with a 'concentration' column and sampleID (normalised) counts + list of sampleIDs,
#         calculate correlation between concentration
#         return a list of 'good samples', i.e. passing the cutoff
#
#         Assert concentration column is present
#         Assert sampleIDs are in control_df columns
#         """
#         concentrations = np.log2(control_df.concentration)
#         samples = control_df[sampleIDs]
#         corr_df = pd.DataFrame(samples.corrwith(concentrations), columns=['R'])
#         corr_df["R2"] = corr_df.R ** 2
#         good_samples = list(corr_df[corr_df.R2 > cutoff].index)
#         return corr_df, good_samples
#
#     def prepare_mageck_dataset(counts_df, sampleData, control_barcodes, annotation_cols, good_samples, name,
#                                batch_col, treatment_col, outDir):
#
#         """
#
#         Assume the first column of sampleData contains sampleIDs.
#         Assume second has geneName
#         The rest are raw counts for samples in sampleIDs.
#
#         """
#
#         batch_file = outDir / f"{name}_batch.txt"
#         count_file = outDir / f"{name}_count.txt"
#         sampleID_col = sampleData.columns[0]
#
#         batch_df = (sampleData[sampleData[sampleID_col].isin(good_samples)]
#                     [[sampleID_col, batch_col, treatment_col]]
#                     .sort_values([treatment_col, batch_col]))
#
#         batch_df.to_csv(batch_file, index=False, sep='\t')
#         magDf = counts_df[annotation_cols + good_samples].copy()
#         magDf.loc[magDf[annotation_cols[0]].isin(control_barcodes), annotation_cols[1]] = 'control'
#         magDf = magDf.dropna(subset=annotation_cols).fillna(0)
#         magDf.to_csv(count_file, index=False, sep='\t')
#         return batch_file, count_file
#
#
#     def batch_correct(outDir, name, r_path="../snippets/batchCorrect.R"):
#         """
#         Given count df only with good samples
#         sample data df (read in and validated somewhere else) with information about batches etc.
#         batch column name
#
#         """
#         count_path = outDir / f"{name}_count.txt"
#         batch_path = outDir / f"{name}_batch.txt"
#         cmd = f'Rscript {r_path} {count_path} {batch_path} {name} {outDir}'
#         print(cmd)
#         r = run_command(cmd.split())
#
#     def get_contrast_samples(sampleData, good_samples, treat_col='day',
#                              treatment='d1', control='d0', sampleID='sampleID'):
#         sDf = sampleData[sampleData[sampleID].isin(good_samples)]
#         controls = ",".join(sDf[sDf[treat_col] == control][sampleID].unique())
#         treats = ",".join(sDf[sDf[treat_col] == treatment][sampleID].unique())
#         return controls, treats
#
#     def run_mageck(count_file, treated, controls, out_prefix, control_barcode_file):
#         """
#         count file could be produced before or after batchcorrection
#
#         """
#         cmd = (f"mageck test -k {count_file} -t {treated} "
#                f"-c {controls}  -n {out_prefix} "
#                f"--control-sgrna {control_barcode_file}  --normcounts-to-file")
#         print(cmd)
#         r = run_command(cmd.split())
#
#     def write_control_barcodes_to_file(wt_barcodes, name, outDir):
#         fname = outDir / f"{name}_wt_barcodes.txt"
#         with open(fname, "w") as fo:
#             for bc in wt_barcodes:
#                 fo.write(f"{bc}\n")
#         return fname
#
#
#
#     def clean_samples(self, merged_count_file, control_file):
#         # Read in merged_count table and get sampleIDs
#         # Calculate cpms
#         # Read in control file
#         # Get control counts
#         #wt_df, full_control_df = get_control_counts(cntrl_df, wt_barcodes, cpms, annotation_cols)
#         # Figure out good samples
#         # corr_df, good_samples = calculate_correlation(wt_df, sampleIDs)
#         pass
#
#     def correlation_plots(self):
#         pass
#
#     def run_analysis(self, counts, sample_data_file, good_samples, wt_df,
#                      annotation_cols, name, batch_col='experiment', treatment_col='day',
#                      contrasts=['d1', 'd2', 'd3', 'd4'], baseline='d0',
#                      sampleID='sampleID',
#                      outDir=scratchDir
#                      ):
#
#         # Read in sample data
#         sampleData = read_in_sample_data(sample_data_file, good_samples)
#         wt_barcodes = wt_df.barcode.values
#         # subset df on only good samples and write out to file
#         batch_file, count_file = prepare_mageck_dataset(counts, sampleData, wt_barcodes,
#                                                         annotation_cols, good_samples, name,
#                                                         batch_col, treatment_col, outDir)
#
#         # run batch correction
#         batch_correct(scratchDir, name, r_path="../snippets/batchCorrect.R")
#
#         fname = write_control_barcodes_to_file(wt_barcodes, name, scratchDir)
#
#         # run MAGeCK RRA for day 1
#         contrasts_ran = []
#         for contrast in contrasts:
#             controls, treat = get_contrast_samples(sampleData, good_samples, treatment_col,
#                                                    contrast, baseline, sampleID)
#             count_file2 = count_file.with_suffix('.batchcorrected.txt')
#             if len(treat) > 0 and len(controls) > 0:
#                 run_mageck(count_file2, treat, controls, scratchDir / f"{name}-{contrast}", fname)
#                 contrasts_ran.append(contrast)
#             else:
#                 continue
#         res = pd.concat([pd.read_table(scratchDir / f"{name}-{i}.gene_summary.txt").assign(treat=i)
#                          for i in contrasts_ran])
#         fres = res[['id', 'num', 'neg|lfc', 'neg|fdr', 'pos|fdr', 'treat']]
#         fres.columns = [annotation_cols[1], 'number_of_barcodes', 'LFC', 'neg_selection_fdr', 'pos_selection_fdr',
#                         'contrast']
#         fres.to_csv(scratchDir / f'{name}_rra_results.csv')
#         return fres
