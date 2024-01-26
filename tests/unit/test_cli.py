import subprocess
import shlex
from pathlib import Path
import pandas as pd


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


def test_cli_analysis_no_batch(analysis_test_data, tmpdir):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    treat_col, bline = 'day',  'd0'
    cmd_str = f'mbarq analyze -i {count_file}  -s {sample_file} ' \
              f'-c {controls} --treatment_column {treat_col} ' \
              f'--baseline {bline} ' \
              f' -o {tmpdir} -g Name -n cli_analysis_test2 '
    subprocess.call(shlex.split(cmd_str))
    actual_results = pd.read_csv(tmpdir.join("cli_analysis_test2_rra_results.csv")).to_dict()
    sample_results = {'Name': {1283: 'sciX', 1955: 'SL1344_4468', 2115: 'xylA', 1704: 'yieM', 3475: 'basR', 920: 'tatA', 2017: 'dcoA', 945: 'lrhA', 2361: 'SL1344_2691', 213: 'SL1344_1264', 1356: 'SL1344_1477', 1931: 'SL1344_1567', 1742: 'STM3026', 2251: 'SL1344_2738', 3285: 'ampH', 2865: 'prgH', 21: 'rfbX', 993: 'pspF', 3558: 'SL1344_0019', 99: 'ilvI'}, 'number_of_barcodes': {1283: 1, 1955: 5, 2115: 1, 1704: 1, 3475: 1, 920: 2, 2017: 1, 945: 1, 2361: 1, 213: 1, 1356: 1, 1931: 4, 1742: 1, 2251: 2, 3285: 1, 2865: 1, 21: 1, 993: 1, 3558: 3, 99: 2}, 'LFC': {1283: 0.35757, 1955: -0.57837, 2115: -0.16625, 1704: 0.72942, 3475: 0.05004, 920: -6.9198, 2017: -0.38311, 945: -2.4561, 2361: 0.14491, 213: 0.078651, 1356: 0.40869, 1931: 0.38251, 1742: 0.7928, 2251: 0.15916, 3285: -0.38851, 2865: -4.9191, 21: -2.7876, 993: -0.13799, 3558: -0.41019, 99: -0.118}, 'neg_selection_fdr': {1283: 1.0, 1955: 0.924097, 2115: 0.924097, 1704: 1.0, 3475: 0.593285, 920: 6.1e-05, 2017: 0.924097, 945: 6.1e-05, 2361: 0.924097, 213: 1.0, 1356: 1.0, 1931: 0.924097, 1742: 1.0, 2251: 0.924097, 3285: 0.257293, 2865: 1.7e-05, 21: 4.8e-05, 993: 1.0, 3558: 0.521566, 99: 4.8e-05}, 'pos_selection_fdr': {1283: 6e-06, 1955: 1.0, 2115: 0.888145, 1704: 6e-06, 3475: 1.0, 920: 1.0, 2017: 0.888145, 945: 1.0, 2361: 0.888145, 213: 0.43025, 1356: 6e-06, 1931: 0.888145, 1742: 6e-06, 2251: 3.6e-05, 3285: 1.0, 2865: 1.0, 21: 1.0, 993: 0.550381, 3558: 1.0, 99: 0.646139}, 'contrast': {1283: 'd2', 1955: 'd3', 2115: 'd3', 1704: 'd2', 3475: 'd4', 920: 'd2', 2017: 'd3', 945: 'd2', 2361: 'd3', 213: 'd1', 1356: 'd2', 1931: 'd3', 1742: 'd2', 2251: 'd3', 3285: 'd4', 2865: 'd4', 21: 'd1', 993: 'd2', 3558: 'd4', 99: 'd1'}}
    print(actual_results)
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])


def test_cli_analysis_no_control(analysis_test_data, tmpdir):
    _, _, _, _, count_file, sample_file, _, _ = analysis_test_data
    treat_col,  bline = 'day',  'd0'
    cmd_str = f'mbarq analyze -i {count_file}  -s {sample_file} ' \
              f'--treatment_column {treat_col} ' \
              f' --baseline {bline} ' \
              f' -o {tmpdir} -g Name -n cli_analysis_test3 '
    subprocess.call(shlex.split(cmd_str))
    actual_results = pd.read_csv(tmpdir.join("cli_analysis_test3_rra_results.csv"))
    sample_results = {'Name': {107: 'gtgA', 342: 'sseI', 2271: 'aroG', 633: 'yjfM', 1506: 'SL1344_3106', 3406: 'SL1344_0330', 2937: 'SL1344_0832', 1326: 'SL1344_0699', 77: 'sipA', 1209: 'tehB', 2692: 'STnc780', 3074: 'SL1344_3750', 2037: 'adi', 811: 'yibK', 1394: 'yrbD', 1734: 'STnc710', 1431: 'SL1344_0702', 71: 'yjgF', 2450: 'ydjM', 843: 'sopD'}, 'number_of_barcodes': {107: 3, 342: 1, 2271: 1, 633: 2, 1506: 2, 3406: 1, 2937: 3, 1326: 3, 77: 4, 1209: 1, 2692: 2, 3074: 3, 2037: 1, 811: 1, 1394: 1, 1734: 2, 1431: 2, 71: 1, 2450: 1, 843: 8}, 'LFC': {107: -0.026237, 342: -0.040027, 2271: -0.012038, 633: 0.008158, 1506: 0.048663, 3406: -3.7883, 2937: -5.5681, 1326: -0.081906, 77: -0.084469, 1209: -0.25546, 2692: 0.25638, 3074: -3.9149, 2037: -0.26044, 811: 0.27488, 1394: -0.032308, 1734: 0.15721, 1431: 0.035906, 71: -0.97829, 2450: 0.15414, 843: 0.021523}, 'neg_selection_fdr': {107: 0.994505, 342: 0.999867, 2271: 0.999995, 633: 0.999867, 1506: 0.999995, 3406: 0.998009, 2937: 0.998009, 1326: 0.999995, 77: 0.915541, 1209: 0.950881, 2692: 0.999995, 3074: 0.998009, 2037: 0.843819, 811: 0.999867, 1394: 0.99656, 1734: 0.999995, 1431: 0.999995, 71: 0.49404, 2450: 0.999995, 843: 0.999867}, 'pos_selection_fdr': {107: 0.999995, 342: 0.999995, 2271: 0.996152, 633: 0.999995, 1506: 0.999995, 3406: 0.970383, 2937: 0.970383, 1326: 0.999995, 77: 0.999995, 1209: 0.999995, 2692: 0.973244, 3074: 0.987142, 2037: 0.999995, 811: 0.999995, 1394: 0.999995, 1734: 0.999995, 1431: 0.999995, 71: 0.999995, 2450: 0.973244, 843: 0.999995}, 'contrast': {107: 'd1', 342: 'd1', 2271: 'd3', 633: 'd1', 1506: 'd2', 3406: 'd4', 2937: 'd4', 1326: 'd2', 77: 'd1', 1209: 'd2', 2692: 'd3', 3074: 'd4', 2037: 'd3', 811: 'd1', 1394: 'd2', 1734: 'd2', 1431: 'd2', 71: 'd1', 2450: 'd3', 843: 'd1'}}
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])
    



def test_cli_analysis_log(analysis_test_data, tmpdir):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    treat_col,  bline = 'day',  'd0'
    cmd_str = f'mbarq analyze -i {count_file}  -s {sample_file} ' \
              f'-c {controls} --treatment_column {treat_col} ' \
              f' --baseline {bline} ' \
              f' -o {tmpdir} -g Name -n cli_analysis_test1 '
    subprocess.call(shlex.split(cmd_str))
