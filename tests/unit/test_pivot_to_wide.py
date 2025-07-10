import pytest
import pandas as pd
import subprocess
import shlex
from mbarq.analysis import Experiment


class TestPivotToWide:
    """Test suite for the pivot_to_wide staticmethod"""

    def test_pivot_to_wide_basic_functionality(self):
        """Test basic pivot functionality with valid data"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneA', 'geneB', 'geneB'],
            'number_of_barcodes': [3, 3, 2, 2],
            'LFC': [1.2, 1.5, -0.8, -0.9],
            'neg_selection_fdr': [0.05, 0.03, 0.8, 0.7],
            'pos_selection_fdr': [0.9, 0.8, 0.1, 0.2],
            'contrast': ['d1', 'd2', 'd1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(test_data, ['Name', 'number_of_barcodes'])
        
        # Verify structure
        assert 'Name' in result.columns
        assert 'number_of_barcodes' in result.columns
        assert 'LFC_d1' in result.columns
        assert 'LFC_d2' in result.columns
        assert 'neg_selection_fdr_d1' in result.columns
        assert 'neg_selection_fdr_d2' in result.columns
        assert 'pos_selection_fdr_d1' in result.columns
        assert 'pos_selection_fdr_d2' in result.columns
        assert len(result) == 2  # Two unique genes
        
        # Verify data values
        gene_a_row = result[result['Name'] == 'geneA'].iloc[0]
        assert gene_a_row['LFC_d1'] == 1.2
        assert gene_a_row['LFC_d2'] == 1.5
        assert gene_a_row['number_of_barcodes'] == 3

    def test_pivot_to_wide_single_index_column_string(self):
        """Test with single string index column"""
        test_data = pd.DataFrame({
            'barcode': ['ATCG', 'ATCG', 'GCTA', 'GCTA'],
            'LFC': [0.5, 0.7, -0.3, -0.1],
            'contrast': ['d1', 'd2', 'd1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(test_data, 'barcode')
        
        assert 'barcode' in result.columns
        assert 'LFC_d1' in result.columns
        assert 'LFC_d2' in result.columns
        assert len(result) == 2

    def test_pivot_to_wide_multiple_index_columns(self):
        """Test with multiple index columns as list"""
        test_data = pd.DataFrame({
            'barcode': ['ATCG', 'ATCG', 'GCTA', 'GCTA'],
            'Name': ['geneA', 'geneA', 'geneB', 'geneB'],
            'LFC': [0.5, 0.7, -0.3, -0.1],
            'contrast': ['d1', 'd2', 'd1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(test_data, ['barcode', 'Name'])
        
        assert 'barcode' in result.columns
        assert 'Name' in result.columns
        assert 'LFC_d1' in result.columns
        assert 'LFC_d2' in result.columns
        assert len(result) == 2

    def test_pivot_to_wide_empty_dataframe(self):
        """Test with empty DataFrame"""
        empty_df = pd.DataFrame()
        result = Experiment.pivot_to_wide(empty_df, ['Name'])
        assert result.empty

    def test_pivot_to_wide_no_contrast_column(self):
        """Test when 'contrast' column is missing"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneB'],
            'LFC': [1.2, -0.8]
        })
        
        with pytest.raises(ValueError, match="Missing required columns: \\['contrast'\\]"):
            Experiment.pivot_to_wide(test_data, ['Name'])

    def test_pivot_to_wide_missing_index_columns(self):
        """Test when specified index columns don't exist"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneB'],
            'LFC': [1.2, -0.8],
            'contrast': ['d1', 'd2']
        })
        
        with pytest.raises(ValueError, match="Missing required columns: \\['NonExistentColumn'\\]"):
            Experiment.pivot_to_wide(test_data, ['NonExistentColumn'])

    def test_pivot_to_wide_duplicate_index_contrast_combinations(self):
        """Test duplicate combinations of index + contrast"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneA', 'geneA'],  # Duplicate geneA + d1 combination
            'LFC': [1.2, 1.5, 1.8],
            'contrast': ['d1', 'd1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(test_data, ['Name'])
        assert result.empty  # Should return empty DataFrame due to duplicate combinations

    def test_pivot_to_wide_column_sorting(self):
        """Test that columns are sorted by contrast values"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneA', 'geneA', 'geneA'],
            'LFC': [1.0, 2.0, 3.0, 4.0],
            'contrast': ['d4', 'd1', 'd3', 'd2']  # Unsorted contrasts
        })
        
        result = Experiment.pivot_to_wide(test_data, ['Name'])
        
        # Get column names that start with 'LFC_'
        lfc_cols = [col for col in result.columns if col.startswith('LFC_')]
        expected_order = ['LFC_d1', 'LFC_d2', 'LFC_d3', 'LFC_d4']
        assert lfc_cols == expected_order

    def test_pivot_to_wide_column_naming(self):
        """Test column names are formatted as {metric}_{contrast}"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneA'],
            'LFC': [1.2, 1.5],
            'p_value': [0.05, 0.03],
            'fdr': [0.1, 0.08],
            'contrast': ['d1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(test_data, ['Name'])
        
        expected_columns = ['Name', 'LFC_d1', 'LFC_d2', 'p_value_d1', 'p_value_d2', 'fdr_d1', 'fdr_d2']
        assert set(result.columns) == set(expected_columns)

    def test_pivot_to_wide_index_reset(self):
        """Test that index is properly reset"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneA', 'geneB', 'geneB'],
            'LFC': [1.2, 1.5, -0.8, -0.9],
            'contrast': ['d1', 'd2', 'd1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(test_data, ['Name'])
        
        # Check that index is default integer index
        assert list(result.index) == [0, 1]
        assert result.index.name is None

    def test_pivot_to_wide_with_gene_results_structure(self):
        """Test with gene result data structure similar to process_results"""
        # Simulate data structure from process_results for genes
        gene_data = pd.DataFrame({
            'Name': ['geneA', 'geneA', 'geneB', 'geneB', 'geneC', 'geneC'],
            'number_of_barcodes': [3, 3, 2, 2, 1, 1],
            'LFC': [1.2, 1.5, -0.8, -0.9, 0.1, 0.2],
            'neg_selection_fdr': [0.05, 0.03, 0.8, 0.7, 0.9, 0.95],
            'pos_selection_fdr': [0.9, 0.8, 0.1, 0.2, 0.85, 0.88],
            'contrast': ['d1', 'd2', 'd1', 'd2', 'd1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(gene_data, ['Name', 'number_of_barcodes'])
        
        # Verify structure matches expected output from process_results
        assert len(result) == 3  # Three unique genes
        assert 'Name' in result.columns
        assert 'number_of_barcodes' in result.columns
        
        # Check that all metric columns are present for each contrast
        for metric in ['LFC', 'neg_selection_fdr', 'pos_selection_fdr']:
            assert f'{metric}_d1' in result.columns
            assert f'{metric}_d2' in result.columns

    def test_pivot_to_wide_with_barcode_results_structure(self):
        """Test with barcode result data structure similar to process_results"""
        # Simulate data structure from process_results for barcodes
        barcode_data = pd.DataFrame({
            'barcode': ['ATCG', 'ATCG', 'GCTA', 'GCTA'],
            'Name': ['geneA', 'geneA', 'geneB', 'geneB'],
            'LFC': [1.1, 1.3, -0.7, -0.5],
            'p_value': [0.01, 0.02, 0.05, 0.08],
            'fdr': [0.05, 0.06, 0.1, 0.12],
            'contrast': ['d1', 'd2', 'd1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(barcode_data, ['barcode', 'Name'])
        
        # Verify structure matches expected output from process_results
        assert len(result) == 2  # Two unique barcode-gene combinations
        assert 'barcode' in result.columns
        assert 'Name' in result.columns
        
        # Check that all metric columns are present for each contrast
        for metric in ['LFC', 'p_value', 'fdr']:
            assert f'{metric}_d1' in result.columns
            assert f'{metric}_d2' in result.columns

    def test_pivot_to_wide_with_numeric_contrasts(self):
        """Test with numeric contrast values"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneA', 'geneB', 'geneB'],
            'LFC': [1.2, 1.5, -0.8, -0.9],
            'contrast': [1, 2, 1, 2]  # Numeric contrasts
        })
        
        result = Experiment.pivot_to_wide(test_data, ['Name'])
        
        assert 'LFC_1' in result.columns
        assert 'LFC_2' in result.columns
        assert len(result) == 2

    def test_pivot_to_wide_preserves_data_types(self):
        """Test that data types are preserved during pivot"""
        test_data = pd.DataFrame({
            'Name': ['geneA', 'geneA'],
            'count': [10, 15],  # Integer
            'LFC': [1.2, 1.5],  # Float
            'significant': [True, False],  # Boolean
            'contrast': ['d1', 'd2']
        })
        
        result = Experiment.pivot_to_wide(test_data, ['Name'])
        
        # Check that data types are preserved
        assert result['count_d1'].dtype == 'int64'
        assert result['LFC_d1'].dtype == 'float64'
        assert result['significant_d1'].dtype == 'bool'


class TestPivotToWideIntegration:
    """Integration tests for pivot_to_wide with CLI"""

    def test_cli_analysis_wide_format(self, analysis_test_data, tmpdir):
        """Test CLI with wide format output"""
        _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
        treat_col, bline = 'day', 'd0'
        
        # Run CLI command with wide format
        cmd_str = f'mbarq analyze -i {count_file} -s {sample_file} ' \
                  f'-c {controls} --treatment_column {treat_col} ' \
                  f'--baseline {bline} --format wide ' \
                  f'-o {tmpdir} -g Name -n cli_wide_test'
        
        subprocess.call(shlex.split(cmd_str))
        
        # Check that output files exist
        rra_file = tmpdir.join("cli_wide_test_rra_results.csv")
        bc_file = tmpdir.join("cli_wide_test_barcodes_results.csv")
        
        assert rra_file.check()
        assert bc_file.check()
        
        # Verify wide format structure for gene results
        rra_df = pd.read_csv(str(rra_file))
        
        # Check that index columns are present
        assert 'Name' in rra_df.columns
        assert 'number_of_barcodes' in rra_df.columns
        
        # Check that wide format columns exist for each contrast
        contrasts = ['d1', 'd2', 'd3', 'd4']  # Based on test data
        for contrast in contrasts:
            assert f'LFC_{contrast}' in rra_df.columns
            assert f'neg_selection_fdr_{contrast}' in rra_df.columns
            assert f'pos_selection_fdr_{contrast}' in rra_df.columns
        
        # Verify no 'contrast' column exists (should be pivoted away)
        assert 'contrast' not in rra_df.columns
        
        # Verify barcode results structure
        bc_df = pd.read_csv(str(bc_file))
        
        # Check that index columns are present
        assert 'barcode' in bc_df.columns
        assert 'Name' in bc_df.columns
        
        # Check that wide format columns exist
        for contrast in contrasts:
            # Note: exact column names depend on what MAGeCK outputs
            # This may need adjustment based on actual barcode result structure
            assert any(col.endswith(f'_{contrast}') for col in bc_df.columns)

    def test_cli_analysis_wide_vs_long_format_comparison(self, analysis_test_data, tmpdir):
        """Compare wide vs long format outputs to ensure data integrity"""
        _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
        treat_col, bline = 'day', 'd0'
        
        # Run with long format (default)
        cmd_long = f'mbarq analyze -i {count_file} -s {sample_file} ' \
                   f'-c {controls} --treatment_column {treat_col} ' \
                   f'--baseline {bline} --format long ' \
                   f'-o {tmpdir} -g Name -n cli_long_test'
        
        subprocess.call(shlex.split(cmd_long))
        
        # Run with wide format
        cmd_wide = f'mbarq analyze -i {count_file} -s {sample_file} ' \
                   f'-c {controls} --treatment_column {treat_col} ' \
                   f'--baseline {bline} --format wide ' \
                   f'-o {tmpdir} -g Name -n cli_wide_test'
        
        subprocess.call(shlex.split(cmd_wide))
        
        # Read both results
        long_df = pd.read_csv(tmpdir.join("cli_long_test_rra_results.csv"))
        wide_df = pd.read_csv(tmpdir.join("cli_wide_test_rra_results.csv"))
        
        # Verify that wide format has fewer rows (genes collapsed across contrasts)
        assert len(wide_df) < len(long_df)
        
        # Verify that wide format has more columns (metrics spread across contrasts)
        assert len(wide_df.columns) > len(long_df.columns)
        
        # Check that the same genes are present in both formats
        long_genes = set(long_df['Name'].unique())
        wide_genes = set(wide_df['Name'].unique())
        assert long_genes == wide_genes
        
        # Verify data consistency by checking a specific gene's data
        test_gene = long_df['Name'].iloc[0]
        long_gene_data = long_df[long_df['Name'] == test_gene]
        wide_gene_data = wide_df[wide_df['Name'] == test_gene]
        
        # Should have one row in wide format for this gene
        assert len(wide_gene_data) == 1
        
        # Check that data values match between formats
        for _, row in long_gene_data.iterrows():
            contrast = row['contrast']
            wide_row = wide_gene_data.iloc[0]
            
            # Compare LFC values
            assert abs(row['LFC'] - wide_row[f'LFC_{contrast}']) < 1e-10

    def test_cli_analysis_wide_format_no_control(self, analysis_test_data, tmpdir):
        """Test wide format CLI without control barcodes"""
        _, _, _, _, count_file, sample_file, _, _ = analysis_test_data
        treat_col, bline = 'day', 'd0'
        
        cmd_str = f'mbarq analyze -i {count_file} -s {sample_file} ' \
                  f'--treatment_column {treat_col} ' \
                  f'--baseline {bline} --format wide ' \
                  f'-o {tmpdir} -g Name -n cli_wide_no_control_test'
        
        subprocess.call(shlex.split(cmd_str))
        
        # Check that output files exist
        rra_file = tmpdir.join("cli_wide_no_control_test_rra_results.csv")
        bc_file = tmpdir.join("cli_wide_no_control_test_barcodes_results.csv")
        
        assert rra_file.check()
        assert bc_file.check()
        
        # Verify wide format structure
        rra_df = pd.read_csv(str(rra_file))
        assert 'Name' in rra_df.columns
        assert 'number_of_barcodes' in rra_df.columns
        
        # Should have wide format columns
        contrasts = ['d1', 'd2', 'd3', 'd4']
        for contrast in contrasts:
            assert f'LFC_{contrast}' in rra_df.columns

    def test_cli_analysis_invalid_format_defaults_to_long(self, analysis_test_data, tmpdir):
        """Test that invalid format parameter defaults to long format"""
        _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
        treat_col, bline = 'day', 'd0'
        
        cmd_str = f'mbarq analyze -i {count_file} -s {sample_file} ' \
                  f'-c {controls} --treatment_column {treat_col} ' \
                  f'--baseline {bline} --format invalid_format ' \
                  f'-o {tmpdir} -g Name -n cli_invalid_format_test'
        
        subprocess.call(shlex.split(cmd_str))
        
        # Check that output files exist
        rra_file = tmpdir.join("cli_invalid_format_test_rra_results.csv")
        
        assert rra_file.check()
        
        # Should default to long format (has contrast column)
        rra_df = pd.read_csv(str(rra_file))
        assert 'contrast' in rra_df.columns  # Long format should have contrast column