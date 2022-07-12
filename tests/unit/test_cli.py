import subprocess
import shlex
from pathlib import Path
from test_utils import assert_files_are_same


def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc


def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value

OUTDIR= "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/Projects_NCCR/hardt/nguyenb/tnseq/scratch/tmp"


def test_cli_analysis(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    treat_col, batch_col, bline = 'day', 'experiment', 'd0'
    cmd_str = f'mbarq analyze -i {count_file}  -s {sample_file} ' \
              f'-c {controls} --treatment_column {treat_col} ' \
              f'--batch_column {batch_col} --baseline {bline} ' \
              f' -o {tmpdir} -g Name -n cli_analysis_test1 '
    subprocess.call(shlex.split(cmd_str))
    actual_rra = tmpdir.join("cli_analysis_test1_rra_results.csv")
    expected_rra = dnaid1315_expected_outcomes / "test_process_results_rra_results.csv"
    assert_files_are_same(actual_rra, expected_rra)


def test_cli_analysis_no_batch(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    treat_col, batch_col, bline = 'day', '', 'd0'
    cmd_str = f'mbarq analyze -i {count_file}  -s {sample_file} ' \
              f'-c {controls} --treatment_column {treat_col} ' \
              f'--baseline {bline} ' \
              f' -o {tmpdir} -g Name -n cli_analysis_test2 '
    subprocess.call(shlex.split(cmd_str))
    actual_rra = tmpdir.join("cli_analysis_test2_rra_results.csv")
    expected_rra = dnaid1315_expected_outcomes / "test_run_experiment_rra_results_no_batch.csv"
    assert_files_are_same(actual_rra, expected_rra)


def test_cli_analysis_no_control(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    treat_col, batch_col, bline = 'day', 'experiment', 'd0'
    cmd_str = f'mbarq analyze -i {count_file}  -s {sample_file} ' \
              f'--treatment_column {treat_col} ' \
              f'--batch_column {batch_col} --baseline {bline} ' \
              f' -o {tmpdir} -g Name -n cli_analysis_test3 '
    subprocess.call(shlex.split(cmd_str))
    actual_rra = tmpdir.join("cli_analysis_test3_rra_results.csv")
    expected_rra = dnaid1315_expected_outcomes / "test_run_experiment_rra_results_no_control.csv"
    assert_files_are_same(actual_rra, expected_rra)



def test_cli_analysis_log(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    treat_col, batch_col, bline = 'day', 'experiment', 'd0'
    cmd_str = f'mbarq analyze -i {count_file}  -s {sample_file} ' \
              f'-c {controls} --treatment_column {treat_col} ' \
              f'--batch_column {batch_col} --baseline {bline} ' \
              f' -o {OUTDIR} -g Name -n cli_analysis_test1 '
    subprocess.call(shlex.split(cmd_str))
