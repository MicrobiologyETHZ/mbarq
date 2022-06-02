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
import subprocess
import shlex

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


class SampleData:
    def __init__(self, sample_data_csv:Union[Path, str], name, treatment_name,
                 baseline: str, batch_name: str = '', output_dir='.'):
        self.sample_data_csv = sample_data_csv
        self.batch_name = batch_name
        self.treatment_name = treatment_name
        self.name = name
        self.baseline = baseline
        self.sampleData = pd.DataFrame()
        self.output_dir = Path(output_dir)
        self.logger = get_logger('sample-data-log', self.output_dir / f"{self.name}_SampleData.log")
        self.contrasts = []
        self.sampleIDs = []

    def read_sample_data_csv(self):
        self.sampleData = pd.read_csv(self.sample_data_csv)
        sampleID = self.sampleData.columns[0]
        schema_dict = {sampleID: Column(str),
                       self.treatment_name: Column(str, [Check(lambda s: any(s.isin([self.baseline])))])}
        if self.batch_name:
            schema_dict[self.batch_name] = Column(str)
        schema = DataFrameSchema(schema_dict, strict=False, coerce=False)
        try:
            schema.validate(self.sampleData)
            self.contrasts = [c for c in self.sampleData[self.treatment_name].unique() if c != self.baseline]
            self.sampleIDs = list(self.sampleData[sampleID].unique())
        except pandera.errors.SchemaError:
            self.logger.error('Invalid sample data frame')
            sys.exit(1)


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


class Experiment:
    def __init__(self, count_file, sample_data_file, control_file, name,
                  gene_column, treatment_column, baseline, batch_column,
                  cutoff=0.8, output_dir='.'):
        self.output_dir = Path(output_dir)
        self.name = name
        self.sampleIDs = []  # don't need this and self.good_samples refactor
        self.logger = get_logger('experiment-log', self.output_dir / f"{self.name}_Experiment.log")
        self.cds = CountDataSet(count_file, name=name, gene_column_name=gene_column, output_dir=output_dir)
        self.cbars = ControlBarcodes(control_file, output_dir=output_dir)
        self.sd = SampleData(sample_data_file, name, treatment_column, baseline, batch_column, output_dir=output_dir)
        self._initiate_experiment()

        self.wt_bc_conc_counts = pd.DataFrame()
        self.control_bc_counts = pd.DataFrame()
        self.good_samples = []
        self.corr_df = pd.DataFrame()
        self.cutoff = cutoff
        self.batch_file = self.output_dir / f"{self.name}_batch.txt"
        self.count_file = self.output_dir / f"{self.name}_count.txt"
        self.mageck_bc_file = self.output_dir/ f"{self.name}_wt_barcodes.txt"
        self.sampleID_col = self.sd.sampleData.columns[0]


    def _initiate_experiment(self):
        self.cds.create_count_table()
        self.cds.calculate_cpms()
        self.cbars.read_control_file()
        self.sd.read_sample_data_csv()
        self.sampleIDs = set(self.cds.sampleIDs).intersection(set(self.sd.sampleIDs))
        if len(self.sampleIDs) == 0:
            self.logger.error('No matching sampleIDs between sample file and count file found')
            sys.exit(1)

    def _get_good_samples(self):

        """
        If there is information about wt bc and concentrations,
        calculate correlation between counts and concentrations,
        only keep samples with R2 > self.cutoff

        Else all samples are assumed to be good samples

        Correlation is calculated between log2(CPMS + 0.5) and log2(concentration values)

        """
        if not self.cbars.wt_bc_conc.empty:
            self.wt_bc_conc_counts = (self.cbars.wt_bc_conc.merge(self.cds.cpms,
                                                                            how='left', on='barcode')
                                      .fillna(-1)) # cpms are in log, min val is log2(0.5)
            concentrations = np.log2(self.wt_bc_conc_counts.concentration)
            samples = self.wt_bc_conc_counts[self.cds.sampleIDs]
            self.corr_df = pd.DataFrame(samples.corrwith(concentrations), columns=['R'])
            self.corr_df["R2"] = self.corr_df.R ** 2
            self.good_samples = [s for s in list(self.corr_df[self.corr_df.R2 > self.cutoff].index) if s in self.sampleIDs]
        else:
            self.good_samples = self.sampleIDs


    def prepare_mageck_dataset(self):
        """
        Assume the first column of sd.sampleData contains sampleIDs.
        Assume first column of cds.count_table has barcodes and second has geneName
        The rest of cds.count_table are raw counts for samples in sampleIDs.
        """
        annotation_cols = list(self.cds.count_table.columns[0:2])
        sample_cols_to_keep = [self.sampleID_col, self.sd.treatment_name]
        if self.sd.batch_name:
            sample_cols_to_keep += [self.sd.batch_name]
        batch_df = (self.sd.sampleData[self.sd.sampleData[self.sampleID_col].isin(self.good_samples)][sample_cols_to_keep]
                    .sort_values([self.sd.treatment_name]))
        batch_df.to_csv(self.batch_file, index=False, sep='\t')
        magDf = self.cds.count_table[annotation_cols + self.good_samples].copy()
        magDf.loc[magDf[annotation_cols[0]].isin(self.cbars.wt_barcodes), annotation_cols[1]] = 'control'
        magDf = magDf.dropna(subset=annotation_cols).fillna(0)
        magDf.to_csv(self.count_file, index=False, sep='\t')


    def batch_correct(self, r_path=Path(__file__).parent.resolve()/"batchCorrect.R"):
        """
        Given count df only with good samples
        sample data df (read in and validated somewhere else) with information about batches etc.
        batch column name

        """
        cmd = f'Rscript {r_path} {self.count_file} {self.batch_file} {self.name} {self.output_dir}'
        self.logger.info(cmd)
        return_code = subprocess.check_call(shlex.split(cmd))
        if return_code != 0:
            self.logger.error(f"Batch Correction exited with return code {return_code}")
            sys.exit(1)
        self.count_file = self.output_dir / f"{self.name}_count.batchcorrected.txt"
        self.logger.info(f"Batch Correction Complete")

    def get_contrast_samples(self, treatment='d1'):
        sDf = self.sd.sampleData[self.sd.sampleData[self.sampleID_col].isin(self.good_samples)]
        controls = ",".join(sDf[sDf[self.sd.treatment_name] == self.sd.baseline][self.sampleID_col].unique())
        treats = ",".join(sDf[sDf[self.sd.treatment_name] == treatment][self.sampleID_col].unique())
        return controls, treats

    def write_control_barcodes_to_file(self):
        with open(self.mageck_bc_file, "w") as fo:
            for bc in self.cbars.wt_barcodes:
                fo.write(f"{bc}\n")


    def run_mageck(self, treated, controls, run_name):
        """
        count file could be produced before or after batchcorrection
        """
        prefix = self.output_dir/run_name
        if len(self.cbars.wt_barcodes) > 0:
            cmd = (f"mageck test -k {self.count_file} -t {treated} "
                   f"-c {controls}  -n {prefix} --norm-method control "
                   f"--control-sgrna {self.mageck_bc_file}  --normcounts-to-file")
        else:
            cmd = (f"mageck test -k {self.count_file} -t {treated} "
                   f"-c {controls}  -n {prefix} --norm-method median "
                   f" --normcounts-to-file")

        self.logger.info(cmd)
        r = subprocess.check_call(cmd.split())

    def run_all_contrasts(self):
        for contrast in self.sd.contrasts:
            self.logger.info(f"Running MAGeCK to compare {contrast} to {self.sd.baseline}")
            controls, treats = self.get_contrast_samples(treatment=contrast)
            run_name = f"{self.name}_{contrast}_vs_{self.sd.baseline}"
            self.run_mageck(treats, controls, run_name)

    def process_results(self):
        res = pd.concat([pd.read_table(self.output_dir/ f"{self.name}_{i}_vs_{self.sd.baseline}.gene_summary.txt")
                        .assign(contrast=i) for i in self.sd.contrasts])
        fres = res[['id', 'num', 'neg|lfc', 'neg|fdr', 'pos|fdr', 'contrast']]
        fres.columns = ['Name', 'number_of_barcodes', 'LFC', 'neg_selection_fdr', 'pos_selection_fdr',
                        'contrast']
        fres.to_csv(self.output_dir / f'{self.name}_rra_results.csv')

    def run_experiment(self):
        self._get_good_samples()
        self.prepare_mageck_dataset()
        if self.sd.batch_name:
            self.batch_correct()
        if len(self.cbars.wt_barcodes) > 0:
            self.write_control_barcodes_to_file()
        self.run_all_contrasts()
        self.process_results()

