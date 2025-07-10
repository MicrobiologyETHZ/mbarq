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
import os

class CountDataSet:
    """
    Merge count files into one or if given one count file validate proper format

    """

    def __init__(self, count_files: Union[str, List[str]], name: str,
                 gene_column_name: str = '', sep=',', output_dir: str = '',
                 logger=""):
        self.output_dir = Path(output_dir)
        self.name = name
        self.logger = logger if logger else get_logger('count-dataset-log',
                                                       self.output_dir / f"{self.name}_CountDataSet.log")
        self.count_files = count_files
        self.gene_name: str = gene_column_name
        self.sep: str = sep
        self.count_table: pd.DataFrame = pd.DataFrame()
        self.cpms: pd.DataFrame = pd.DataFrame()
        self.sampleIDs: List[str] = []

    def _merge_count_files(self) -> pd.DataFrame:
        self.logger.info("Merging mbarq count files")
        df_list = []
        sampleIDs = []
        index_cols = []
        for count_file in self.count_files:
            self.logger.info(f"Processing {count_file}")
            sampleID = Path(count_file).stem.split("_mbarq")[0]
            df = pd.read_table(Path(count_file), sep=self.sep)
            if not self.gene_name:
                df = df[['barcode', 'barcode_count']]
                index_cols = ['barcode']
            elif self.gene_name and self.gene_name in df.columns:
                df = df[['barcode', 'barcode_count', self.gene_name]]
                index_cols = ['barcode', self.gene_name]
            else:
                self.logger.warning(
                    f'Gene name column provided is not found in the {Path(count_file).name}. Skipping file')
                continue
            df = df.assign(sampleID=sampleID).drop_duplicates()  # todo fix bug in counting that produces duplicates
            df_list.append(df)
            sampleIDs.append(sampleID)
        if len(df_list) == 0 | len(index_cols) == 0:
            self.logger.error(f'Could not process any of the files. Check that {self.gene_name} is correct')
            sys.exit(1)
        self.logger.info(f"Completed processing count files. Merging.")
        fdf = pd.concat(df_list)
        fdf = (fdf.pivot(index=index_cols, columns='sampleID', values=fdf.columns[1])
               .fillna(0)
               .reset_index())
        self.logger.info('Merging count files complete')
        return fdf

    def _validate_count_table(self, df_to_validate: pd.DataFrame) -> (pd.DataFrame, List[str]):
        self.logger.info('Validating count dataset format')
        # first needs to be a string, second also string if self.gene_name is provided, nulls should not be allowed
        # rest should be floats or ints
        sample_cols = df_to_validate.columns[2:] if self.gene_name else df_to_validate.columns[1:]
        barcode_val = {df_to_validate.columns[0]: Column(str)}
        gene_val = {df_to_validate.columns[1]: Column(str, nullable=True)} if self.gene_name else {}
        sample_val = {sample_col: Column(float, coerce=True) for sample_col in sample_cols}
        validation_dict = {**barcode_val, **gene_val, **sample_val}
        schema = DataFrameSchema(
            validation_dict,
            strict=True,
            coerce=False,
        )
        try:
            schema.validate(df_to_validate)
            self.logger.info("Validation complete.")
            return df_to_validate, sample_cols
        except pandera.errors.SchemaError as e:
            self.logger.error('Validation failed. Invalid count dataset.')
            self.logger.error(e)
            sys.exit(1)

    def create_count_table(self, annotated_only=False):
        if type(self.count_files) == list:
            self.count_table, self.sampleIDs = self._validate_count_table(self._merge_count_files())
            self.logger.info(f'Writing merged count dataset')
            if annotated_only:
                to_file = self.count_table[~self.count_table.iloc[:, 1].isna()]
                to_file.to_csv(self.output_dir / f"{self.name}_mbarq_merged_counts.annotated.csv", index=False)
            else:
                sort_by = self.count_table.columns[1]
                self.count_table.sort_values(sort_by).to_csv(self.output_dir / f"{self.name}_mbarq_merged_counts.csv", index=False)
        else:
            self.logger.info(f'Reading count data from {self.count_files}')
            self.count_table, self.sampleIDs = self._validate_count_table(pd.read_table(self.count_files, sep=self.sep))

    def calculate_cpms(self):
        self.cpms = self.count_table[self.count_table.sum(axis=1, numeric_only=True) > 10]
        # Normalized for library depth and log transform
        self.cpms = self.cpms.apply(
            lambda x: np.log2(x / x.sum() * 1000000 + 0.5) if np.issubdtype(x.dtype, np.number) else x)


class SampleData:
    """
    Process, validate and store relevant info from the metadata file

    """
    def __init__(self, sample_data_csv: Union[Path, str], name, treatment_name,
                 baseline: str, batch_name: str = '', output_dir='.', logger=""):
        self.sample_data_csv = sample_data_csv
        self.batch_name = batch_name
        self.treatment_name = treatment_name
        self.name = name
        self.baseline = baseline
        self.sampleData = pd.DataFrame()
        self.output_dir = Path(output_dir)
        self.logger = logger if logger else get_logger('sample-data-log',
                                                       self.output_dir / f"{self.name}_SampleData.log")
        self.contrasts = []
        self.sampleIDs = []

    def read_sample_data_csv(self):
        self.logger.info(f'Processing sample data file: {self.sample_data_csv}')
        self.sampleData = pd.read_csv(self.sample_data_csv)
        sampleID = self.sampleData.columns[0]
        schema_dict = {sampleID: Column(str),
                       self.treatment_name: Column(str, [Check(lambda s: any(s.isin([self.baseline])))])}
        if self.batch_name:
            schema_dict[self.batch_name] = Column(str)
        schema = DataFrameSchema(schema_dict, strict=False, coerce=False)
        self.logger.info("Validating sample data format.")
        try:
            schema.validate(self.sampleData)
            self.contrasts = [c for c in self.sampleData[self.treatment_name].unique() if c != self.baseline]
            self.sampleIDs = list(self.sampleData[sampleID].unique())
            self.logger.info("Validation complete.")
        except pandera.errors.SchemaError:
            self.logger.error('Validation failed: '
                              'Invalid sample data file. '
                              'Please ensure `--treatment_columns` '
                              'and/or --`batch_column` are found in the file.')
            sys.exit(1)


class ControlBarcodes:
    """
            Control file: No header

            [barcode],[conc]*,[phenotype]*
            * second and third columns are optional
            1. Check that the first column contains barcodes
            2. Check that there are at least 4 barcodes
            3. Check that second column is numeric -> should contain concentrations
            4. If there is a third column it should be string, and at least one should be ['wt', 'WT', 'wildtype', 'wild-type']

            """

    def __init__(self, barcode_csv: Union[Path, str] = '', output_dir: Union[Path, str] = '.',
                 logger=""):
        self.barcode_csv = barcode_csv
        self.wt_barcodes: List = []
        self.wt_bc_conc: pd.DataFrame = pd.DataFrame()
        self.control_bc_conc: pd.DataFrame = pd.DataFrame()
        self.output_dir: Path = Path(output_dir)
        self.logger = logger if logger else get_logger('control-barcode-log', self.output_dir / f"ControlBarcode.log")

    def read_control_file(self):
        if self.barcode_csv:
            self.logger.info('Processing and validating control barcodes file.')
            cntrl_df = pd.read_csv(self.barcode_csv, header=None)
            num_cols = cntrl_df.shape[1]
            col_names = ['barcode', 'concentration', 'genotype']
            schema_dict = {"barcode": Column(str, [Check(lambda s: s.nunique() >= 4)]),
                           "concentration": Column(float),
                           "genotype": Column(str, [
                               Check(lambda s: any(s.isin(['wt', 'WT', 'wildtype', 'wild-type']))),
                           ])}
            cntrl_df.columns = col_names[0:num_cols]
            control_schema = DataFrameSchema({c: schema_dict[c] for c in list(cntrl_df.columns)},
                                             strict=True,
                                             coerce=False,
                                             )
            try:
                control_schema.validate(cntrl_df)
                self.logger.info('Validation complete.')
                if 'genotype' in cntrl_df.columns:
                    self.wt_bc_conc = cntrl_df[cntrl_df.genotype.isin(['wt', 'WT', 'wildtype'])].copy()
                    self.control_bc_conc = cntrl_df.copy()
                else:
                    self.wt_bc_conc = cntrl_df.copy()
                self.wt_barcodes = self.wt_bc_conc.barcode.values
            except pandera.errors.SchemaError:
                self.logger.error('Control file validation failed: invalid control file format or insufficient number of control barcodes (must be at least 4).')
                sys.exit(1)
        else:
            self.logger.info("No information on control barcodes provided.")


class Experiment:
    def __init__(self, count_file, sample_data_file, control_file, name,
                 gene_column, treatment_column, baseline, batch_column='',
                 cutoff=0.8, output_dir='.'):
        self.output_dir = Path(output_dir)
        self.name = name if name else Path(count_file).stem
        self.sampleIDs = []
        self.logger = get_logger('experiment-log', self.output_dir / f"{self.name}_Experiment.log")
        self.gene_column = gene_column
        self.cds = CountDataSet(count_file, name=name, gene_column_name=self.gene_column,
                                output_dir=output_dir, logger=self.logger)
        self.cbars = ControlBarcodes(control_file, output_dir=output_dir, logger=self.logger)
        self.sd = SampleData(sample_data_file, name, treatment_column, baseline,
                             batch_column, output_dir=output_dir, logger=self.logger)
        self._initiate_experiment()
        self.annotation_cols = list(self.cds.count_table.columns[0:2])
        self.wt_bc_conc_counts = pd.DataFrame()
        self.control_bc_counts = pd.DataFrame()
        self.corr_df = pd.DataFrame()
        self.cutoff = cutoff
        self.batch_file = self.output_dir / f"{self.name}_batch.txt"
        self.count_file = self.output_dir / f"{self.name}_count.txt"
        self.experiment_counts = pd.DataFrame()
        self.mageck_bc_file = self.output_dir / f"{self.name}_wt_barcodes.txt"
        self.sampleID_col = self.sd.sampleData.columns[0]

    def _initiate_experiment(self):
        """
        Load count and sample and control barcode data


        """
        self.logger.info('Initiating Experiment.')
        self.logger.info("Loading count data.")
        self.cds.create_count_table()
        self.cds.calculate_cpms()
        self.logger.info("Loading control barcode data.")
        self.cbars.read_control_file()
        self.logger.info("Loading sample data.")
        self.sd.read_sample_data_csv()
        self.sampleIDs = list(set(self.cds.sampleIDs).intersection(set(self.sd.sampleIDs)))
        self.logger.info(f"Loaded data on {len(self.sampleIDs)} samples.")
        if len(self.sampleIDs) == 0:
            self.logger.error('No matching sampleIDs between sample file and count file found')
            sys.exit(1)

    def _get_good_samples(self):

        """
        If there is information about wt bc and concentrations,
        calculate correlation between counts (log of cpms) and  log of concentrations,
        only keep samples with R2 > self.cutoff

        Else all samples are assumed to be good samples

        Correlation is calculated between log2(CPMS + 0.5) and log2(concentration values)

        """
        if not self.cbars.wt_bc_conc.empty:
            self.logger.info(f"Calculating correlations between control concentrations and control counts (CPMs)")
            self.logger.info(f"Filtering out samples with R2 < {self.cutoff}")
            # Merge counts with control barcodes, if any barcodes are missing from counts assume count is 0
            # Since forcing there to be at least 4 barcodes, will always be based at least on 4 points
            self.wt_bc_conc_counts = (self.cbars.wt_bc_conc.merge(self.cds.cpms, how='left', on='barcode')
                                      .fillna(-1))  # cpms are in log2, so zero values -> log2(0 +0.5) = -1
            good_samples = [s for s in self.cds.sampleIDs if (self.wt_bc_conc_counts[s] == -1).sum() / self.wt_bc_conc_counts.shape[0] <= 0.2]
            samples = self.wt_bc_conc_counts[good_samples]
            concentrations = np.log2(self.wt_bc_conc_counts.concentration)
            self.corr_df = pd.DataFrame(samples.corrwith(concentrations), columns=['R'])
            self.corr_df["R2"] = self.corr_df.R ** 2
            self.sampleIDs = [s for s in self.sampleIDs if s in list(self.corr_df[self.corr_df.R2 > self.cutoff].index)]
            self.corr_df.to_csv(self.output_dir/f"{self.name}_correlations.csv")
            self.wt_bc_conc_counts.to_csv(self.output_dir/f"{self.name}_control_counts.csv")
        else:
            self.logger.info('No information on control barcode concentration provided. Keeping all samples')
        self.logger.info(f'Proceeding with analysis of {len(self.sampleIDs)} samples')

    def prepare_mageck_dataset(self):
        """
        Assume the first column of sd.sampleData contains sampleIDs.
        Assume first column of cds.count_table has barcodes and second has geneName
        The rest of cds.count_table are raw counts for samples in sampleIDs.
        """
        sample_cols_to_keep = [self.sampleID_col, self.sd.treatment_name]
        # if self.sd.batch_name:
        #     sample_cols_to_keep += [self.sd.batch_name]
        self.sd.sampleData = (self.sd.sampleData[self.sd.sampleData[self.sampleID_col].isin(self.sampleIDs)][sample_cols_to_keep]
                             .sort_values([self.sd.treatment_name]))
        self.sd.sampleData.to_csv(self.batch_file, index=False, sep='\t')
        magDf = self.cds.count_table[self.annotation_cols + self.sampleIDs].copy()
        magDf.loc[magDf[self.annotation_cols[0]].isin(self.cbars.wt_barcodes), self.annotation_cols[1]] = 'control'
        self.experiment_counts = magDf.dropna(subset=self.annotation_cols).fillna(0)
        #self.experiment_counts.to_csv(self.count_file, index=False, sep='\t') # todo get rid of this

    # def batch_correct(self, r_path=Path(__file__).parent.resolve() / "batchCorrect.R"):
    #     """
    #     Given count df only with good samples
    #     sample data df (read in and validated somewhere else) with information about batches etc.
    #     batch column name
    #
    #     """
    #     self.logger.info('Running batch correction')
    #     cmd = f'Rscript {r_path} {self.count_file} {self.batch_file} {self.name} {self.output_dir}'
    #     self.logger.info(f"R command: {cmd}")
    #     return_code = subprocess.check_call(shlex.split(cmd))
    #     if return_code != 0:
    #         self.logger.error(f"Batch Correction exited with return code {return_code}")
    #         sys.exit(1)
    #     self.count_file = self.output_dir / f"{self.name}_count.batchcorrected.txt"
    #     self.experiment_counts = pd.read_table(self.count_file)
    #     self.logger.info(f"Batch correction complete.")

    def get_contrast_samples(self, treatment='d1', filter_low_counts=0):
        controls = list(self.sd.sampleData[self.sd.sampleData[self.sd.treatment_name] == self.sd.baseline][self.sampleID_col].unique())
        treats = list(self.sd.sampleData[self.sd.sampleData[self.sd.treatment_name] == treatment][self.sampleID_col].unique())
        contrast_table = self.experiment_counts[self.annotation_cols + controls + treats]
        if filter_low_counts:
            filter_low_counts = max(filter_low_counts, len(controls+treats))
            self.logger.info(f"Filtering counts < {filter_low_counts}")
            contrast_table = contrast_table[contrast_table[controls+treats].sum(axis=1) > filter_low_counts]
        return ",".join(controls), ",".join(treats), contrast_table

    def write_control_barcodes_to_file(self):
        with open(self.mageck_bc_file, "w") as fo:
            for bc in self.cbars.wt_barcodes:
                fo.write(f"{bc}\n")

    def run_mageck(self, treated, controls, contrast_table, run_name, normalize_by='', run=True):
        """
        count file could be produced before or after batch correction
        """
        self.logger.info("Running MAGeCK")
        prefix = self.output_dir / run_name
        contrast_file = self.output_dir / f"{run_name}.txt"
        contrast_table.to_csv(contrast_file, index=False, sep='\t')
        if normalize_by and normalize_by not in ('control', 'median', 'total'):
            self.logger.warning('Normalization method not valid. Defaulting to "median"')
            norm_method = 'median'
        elif normalize_by in ['control', 'median', 'total']:
            norm_method = normalize_by
        elif normalize_by == 'control' and len(self.cbars.wt_barcodes) == 0:
            self.logger.warning('No control barcodes found. Defaulting normalization mode to "median"')
            norm_method = 'median'
        elif len(self.cbars.wt_barcodes) > 0:
            norm_method = 'control'
        else:
            norm_method = 'median'
        self.logger.info(f"Normalization method: {norm_method}")
        additional_args = f" --control-sgrna {self.mageck_bc_file} " if norm_method == 'control' else ''
        cmd = (f"mageck test -k {contrast_file} -t {treated} "
               f"-c {controls}  -n {prefix} --norm-method {norm_method} " + additional_args)
        self.logger.info(f"MAGeCK command: {cmd}")
        if run:
            r = subprocess.check_call(cmd.split())
            if r == 0:
                os.remove(contrast_file)
        return cmd

    def run_all_contrasts(self, normalize_by='', filter_low_counts=0):
        contrasts_run = []
        for contrast in self.sd.contrasts:
            self.logger.info(f"Comparing {contrast} to {self.sd.baseline}")
            controls, treats, contrast_table = self.get_contrast_samples(treatment=contrast,
                                                                         filter_low_counts=filter_low_counts)
            if controls and treats:
                run_name = f"{self.name}_{contrast}_vs_{self.sd.baseline}"
                self.run_mageck(treats, controls, contrast_table, run_name, normalize_by)
                contrasts_run.append(contrast)
            else:
                self.logger.warning(f"No sufficient samples found for {contrast}")
        return contrasts_run
    
    @staticmethod
    def pivot_to_wide(long_data, index_cols):
        # Check if dataframe is empty
        if long_data.empty:
            print("Warning: DataFrame is empty")
            return pd.DataFrame()
        
        # Ensure index_cols is a list
        if isinstance(index_cols, str):
            index_cols = [index_cols]
        
        # Check if required columns exist (index_cols + contrast)
        required_cols = index_cols + ['contrast']
        missing_cols = [col for col in required_cols if col not in long_data.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        try:
            # Pivot the data
            wide_data = long_data.pivot(index=index_cols, columns='contrast')
            
            # Handle case where pivot results in empty dataframe
            if wide_data.empty:
                print("Warning: Pivot resulted in empty DataFrame")
                return pd.DataFrame()
            
            # Sort columns by contrast values
            wide_data = wide_data[sorted(wide_data.columns, key=lambda x: x[1])]
            
            # Flatten column names
            wide_data.columns = [f"{c[0]}_{c[1]}" for c in wide_data.columns]
            
            # Reset index
            wide_data = wide_data.reset_index()
            
            return wide_data
            
        except ValueError as e:
            if "Index contains duplicate entries" in str(e):
                print(f"Error: Duplicate combinations found in index columns + contrast")
                print("Consider using pivot_table with an aggregation function instead")
            else:
                print(f"Pivot error: {e}")
            return pd.DataFrame()
        
        except Exception as e:
            print(f"Unexpected error during pivot: {e}")
            return pd.DataFrame()

    def process_results(self, contrasts_run=(), format='long'):
        self.logger.info('Writing out final results.')
        if not contrasts_run:
            contrasts_run = self.sd.contrasts
        contrast_dfs = []
        for con in contrasts_run:
            try:
                contrast_dfs.append(pd.read_table(self.output_dir / f"{self.name}_{con}_vs_{self.sd.baseline}.gene_summary.txt")
                                .assign(contrast=con))
            except:
                print(f"Warning: No such file or directory: {self.output_dir}/{self.name}_{con}_vs_{self.sd.baseline}.gene_summary.txt")
        res = pd.concat(contrast_dfs)
        bc_res = pd.concat([pd.read_table(self.output_dir / f"{self.name}_{i}_vs_{self.sd.baseline}.sgrna_summary.txt")
                        .assign(contrast=i) for i in contrasts_run]).rename({'sgrna': 'barcode', 'Gene': 'Name'}, axis=1)
        fres = res[['id', 'num', 'neg|lfc', 'neg|fdr', 'pos|fdr', 'contrast']]
        fres.columns = [self.gene_column, 'number_of_barcodes', 'LFC', 'neg_selection_fdr', 'pos_selection_fdr',
                        'contrast']
        # Insert pivot to wide if desired
        if format == 'long':
            self.logger.info("Saving output data in long format")
        elif format == 'wide':
            self.logger.info("Saving output data in wide format")
            bc_res = self.pivot_to_wide(bc_res, ['barcode', 'Name'])
            fres = self.pivot_to_wide(fres, [self.gene_column, 'number_of_barcodes'])
        else:
            self.logger.info("Unknown format. Saving output data in long format")
        fres.to_csv(self.output_dir / f'{self.name}_rra_results.csv', index=False)
        bc_res.to_csv(self.output_dir / f'{self.name}_barcodes_results.csv', index=False)


    def run_experiment(self, normalize_by='', filter_low_counts=0, format='long'):
        self.logger.info("Identifying samples")
        self._get_good_samples()
        self.logger.info("Preparing dataset")
        self.prepare_mageck_dataset()
        # if self.sd.batch_name:
        #     self.logger.info(f"Batch column: {self.sd.batch_name}")
        #     self.batch_correct()
        if len(self.cbars.wt_barcodes) > 0:
            self.write_control_barcodes_to_file()
        contrasts_run = self.run_all_contrasts(normalize_by=normalize_by, filter_low_counts=filter_low_counts)
        self.process_results(contrasts_run, format=format)
