from mbarq.analysis import CountDataSet, ControlBarcodes, Experiment, SampleData
import pandas as pd
import pytest


# Testing controlBarcodes

WT_BCS = [
    "TACCCAGAGCACACTCA",
    "ATCCGCGTCACCGAAAA",
    "ACAGAGCTCGGGAGTCT",
    "ACTACAAGACTGGTTAA",
    "AGATGCATGACTAGCTA",
    "AGAATGACCCGGAGGCT",
    "AGGAAGGCGACGAAATC",
    "AGTCATCGATGCTATAT",
    "TAAGTCCGGGCTAAGTC",
    "AACAACACGGTAAGCAA",
    "TATAACACCCCCGATTC",
    "CTACGACAGGGACTTAA",
    "GTGTATAGCAGGAACCC",
    "CCGACGACTGATTGTCC",
    "TCTCACGCAGCGTTTCG",
]


def test__read_control_file_3cols(analysis_test_data, tmpdir):
    _, _, r1, _, _, _, _, _ = analysis_test_data
    cntrlBC = ControlBarcodes(r1, tmpdir)
    cntrlBC.read_control_file()
    assert len(cntrlBC.wt_barcodes == len(WT_BCS))
    assert all(cntrlBC.wt_barcodes == WT_BCS)


def test__read_control_file_2cols(analysis_test_data, tmpdir):
    _, r1, _, _, _, _, _, _ = analysis_test_data
    cntrlBC = ControlBarcodes(r1, tmpdir)
    cntrlBC.read_control_file()
    assert len(cntrlBC.wt_barcodes == len(WT_BCS))
    assert all(cntrlBC.wt_barcodes == WT_BCS)


def test__read_control_file_1cols(analysis_test_data, tmpdir):
    r1, _, _, _, _, _, _, _ = analysis_test_data
    cntrlBC = ControlBarcodes(r1, tmpdir)
    cntrlBC.read_control_file()
    assert len(cntrlBC.wt_barcodes == len(WT_BCS))
    assert all(cntrlBC.wt_barcodes == WT_BCS)


def test__read_control_file_3cols_no_wt(analysis_test_data, tmpdir):
    _, _, _, _, _, _, r1, _ = analysis_test_data
    cntrlBC = ControlBarcodes(r1, tmpdir)
    with pytest.raises(SystemExit):
        cntrlBC.read_control_file()


def test__read_control_file_3bc(analysis_test_data, tmpdir):
    _, _, _, _, _, _, _, r1 = analysis_test_data
    cntrlBC = ControlBarcodes(r1, tmpdir)
    with pytest.raises(SystemExit):
        cntrlBC.read_control_file()


# Testing sampleData


def test_sample_data(analysis_test_data, tmpdir):
    _, _, _, _, _, sample_file, _, _ = analysis_test_data
    sd = SampleData(sample_file, "Test", "day", "d0", "experiment", tmpdir)
    sd.read_sample_data_csv()
    expected_sampleIDs = [
        "dnaid1315_10",
        "dnaid1315_17",
        "dnaid1315_18",
        "dnaid1315_19",
        "dnaid1315_20",
        "dnaid1315_28",
        "dnaid1315_40",
        "dnaid1315_42",
        "dnaid1315_50",
        "dnaid1315_52",
        "dnaid1315_53",
        "dnaid1315_63",
        "dnaid1315_66",
        "dnaid1315_81",
        "dnaid1315_90",
        "dnaid1315_92",
        "dnaid1315_94",
        "dnaid1315_96",
        "dnaid1315_107",
        "dnaid1315_117",
        "dnaid1315_124",
        "dnaid1315_128",
        "dnaid1315_129",
        "dnaid1315_131",
        "dnaid1315_136",
        "dnaid1315_137",
        "dnaid1315_139",
    ]
    expected_contrasts = ["d1", "d2", "d3", "d4"]
    assert sd.contrasts == expected_contrasts
    assert sd.sampleIDs == expected_sampleIDs


def test_sample_data_fail(analysis_test_data, tmpdir):
    _, _, _, _, _, sample_file, _, _ = analysis_test_data
    # batch column does not exist
    sd = SampleData(sample_file, "Test", "day", "d0", "batch", tmpdir)
    with pytest.raises(SystemExit):
        sd.read_sample_data_csv()
    # baseline does not exist
    sd = SampleData(sample_file, "Test", "day", "d9", "experiment", tmpdir)
    with pytest.raises(SystemExit):
        sd.read_sample_data_csv()
    sd = SampleData(sample_file, "Test", "cucumber", "d0", "experiment", tmpdir)
    with pytest.raises(SystemExit):
        sd.read_sample_data_csv()


# Testing experiment


def test_get_good_samples(analysis_test_data, tmpdir):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "TestExp"
    exp = Experiment(
        count_file,
        sample_file,
        controls,
        name,
        "Name",
        "day",
        "d0",
        "experiment",
        0.8,
        tmpdir,
    )
    exp._get_good_samples()
    expected_good_samples = [
        "dnaid1315_10",
        "dnaid1315_107",
        "dnaid1315_117",
        "dnaid1315_124",
        "dnaid1315_128",
        "dnaid1315_129",
        "dnaid1315_131",
        "dnaid1315_136",
        "dnaid1315_17",
        "dnaid1315_18",
        "dnaid1315_19",
        "dnaid1315_20",
        "dnaid1315_28",
        "dnaid1315_40",
        "dnaid1315_42",
        "dnaid1315_50",
        "dnaid1315_52",
        "dnaid1315_66",
        "dnaid1315_81",
        "dnaid1315_90",
        "dnaid1315_92",
        "dnaid1315_94",
    ]
    assert sorted(exp.sampleIDs) == sorted(expected_good_samples)


def test_prepare_mageck_dataset(analysis_test_data, tmpdir):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "test_prepare_mageck_dataset"
    exp = Experiment(
        count_file,
        sample_file,
        controls,
        name,
        "Name",
        "day",
        "d0",
        "experiment",
        0.8,
        tmpdir,
    )
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    expected_good_samples = [
        "dnaid1315_10",
        "dnaid1315_107",
        "dnaid1315_117",
        "dnaid1315_124",
        "dnaid1315_128",
        "dnaid1315_129",
        "dnaid1315_131",
        "dnaid1315_136",
        "dnaid1315_17",
        "dnaid1315_18",
        "dnaid1315_19",
        "dnaid1315_20",
        "dnaid1315_28",
        "dnaid1315_40",
        "dnaid1315_42",
        "dnaid1315_50",
        "dnaid1315_52",
        "dnaid1315_66",
        "dnaid1315_81",
        "dnaid1315_90",
        "dnaid1315_92",
        "dnaid1315_94",
    ]
    assert exp.batch_file.is_file()
    out_samples = list(pd.read_table(exp.batch_file)["sampleID"].values)
    assert sorted(out_samples) == sorted(expected_good_samples)


def test_get_contrast_samples(analysis_test_data, tmpdir):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "test_get_contrast_samples"
    exp = Experiment(
        count_file,
        sample_file,
        controls,
        name,
        "Name",
        "day",
        "d0",
        "experiment",
        0.8,
        tmpdir,
    )
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    treatment = "d1"
    controls, treats, contrast_table = exp.get_contrast_samples(treatment)
    assert controls == "dnaid1315_10,dnaid1315_81,dnaid1315_107"
    assert (
        treats
        == "dnaid1315_17,dnaid1315_18,dnaid1315_19,dnaid1315_124,dnaid1315_117,dnaid1315_90"
    )


def test_run_mageck(analysis_test_data, tmpdir):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "test_run_mageck"
    exp = Experiment(
        count_file,
        sample_file,
        controls,
        name,
        "Name",
        "day",
        "d0",
        "experiment",
        0.8,
        tmpdir,
    )
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    exp.write_control_barcodes_to_file()
    treatment = "d1"
    controls, treats, contrast_table = exp.get_contrast_samples(treatment)
    run_name = "test_run_mageck_d1_vs_d0"
    exp.run_mageck(treats, controls, contrast_table, run_name)
    actual_gene_summary = tmpdir.join(f"{run_name}.gene_summary.txt")
    actual_results = pd.read_table(actual_gene_summary).to_dict()
    sample_results = {
        "id": {
            40: "phoB",
            333: "shfB",
            775: "SL1344_2249",
            744: "bigA",
            114: "cib",
            856: "SL1344_2509",
            211: "traA",
            16: "rfaF",
            76: "fimW",
            602: "SL1344_0490",
        },
        "num": {
            40: 1,
            333: 1,
            775: 1,
            744: 3,
            114: 4,
            856: 1,
            211: 4,
            16: 1,
            76: 2,
            602: 2,
        },
        "neg|score": {
            40: 0.034796,
            333: 0.41223,
            775: 0.88245,
            744: 0.85883,
            114: 0.14777,
            856: 0.95078,
            211: 0.26368,
            16: 0.0065831,
            76: 0.094818,
            602: 0.71299,
        },
        "neg|p-value": {
            40: 4.9515e-06,
            333: 0.83111,
            775: 1.0,
            744: 0.99554,
            114: 4.9515e-06,
            856: 1.0,
            211: 4.9515e-06,
            16: 4.9515e-06,
            76: 4.9515e-06,
            602: 0.97201,
        },
        "neg|fdr": {
            40: 4.8e-05,
            333: 1.0,
            775: 1.0,
            744: 1.0,
            114: 4.8e-05,
            856: 1.0,
            211: 4.8e-05,
            16: 4.8e-05,
            76: 4.8e-05,
            602: 1.0,
        },
        "neg|rank": {
            40: 41,
            333: 334,
            775: 776,
            744: 745,
            114: 115,
            856: 857,
            211: 212,
            16: 17,
            76: 77,
            602: 603,
        },
        "neg|goodsgrna": {
            40: 1,
            333: 0,
            775: 0,
            744: 0,
            114: 1,
            856: 0,
            211: 0,
            16: 1,
            76: 1,
            602: 0,
        },
        "neg|lfc": {
            40: -0.81048,
            333: 0.17386,
            775: 0.49644,
            744: 0.35601,
            114: -0.18279,
            856: 0.57625,
            211: -0.25486,
            16: -4.8318,
            76: -0.24614,
            602: 0.24554,
        },
        "pos|score": {
            40: 0.9652,
            333: 0.58777,
            775: 0.11755,
            744: 0.4143,
            114: 0.99775,
            856: 0.049216,
            211: 0.99986,
            16: 0.99342,
            76: 0.96591,
            602: 0.64487,
        },
        "pos|p-value": {
            40: 1.0,
            333: 0.16889,
            775: 4.9515e-06,
            744: 4.9515e-06,
            114: 0.8007,
            856: 4.9515e-06,
            211: 0.98802,
            16: 1.0,
            76: 0.75033,
            602: 0.3058,
        },
        "pos|fdr": {
            40: 1.0,
            333: 0.261893,
            775: 1.1e-05,
            744: 1.1e-05,
            114: 0.885591,
            856: 1.1e-05,
            211: 1.0,
            16: 1.0,
            76: 0.832906,
            602: 0.43025,
        },
        "pos|rank": {
            40: 850,
            333: 514,
            775: 101,
            744: 360,
            114: 896,
            856: 33,
            211: 907,
            16: 884,
            76: 852,
            602: 563,
        },
        "pos|goodsgrna": {
            40: 0,
            333: 0,
            775: 0,
            744: 0,
            114: 0,
            856: 0,
            211: 0,
            16: 0,
            76: 0,
            602: 0,
        },
        "pos|lfc": {
            40: -0.81048,
            333: 0.17386,
            775: 0.49644,
            744: 0.35601,
            114: -0.18279,
            856: 0.57625,
            211: -0.25486,
            16: -4.8318,
            76: -0.24614,
            602: 0.24554,
        },
    }
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])


def test_run_mageck_norm_methods(
    analysis_test_data, tmpdir
):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "test_run_mageck_norm_methods"
    exp = Experiment(
        count_file,
        sample_file,
        controls,
        name,
        "Name",
        "day",
        "d0",
        "experiment",
        0.8,
        tmpdir,
    )
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    exp.write_control_barcodes_to_file()
    treatment = "d1"
    controls, treats, contrast_table = exp.get_contrast_samples(treatment)
    run_name = "test_run_mageck_d1_vs_d0"
    cmd = exp.run_mageck(treats, controls, contrast_table, run_name, run=False)
    assert "--norm-method control" in cmd
    assert "--control-sgrna" in cmd
    cmd = exp.run_mageck(
        treats, controls, contrast_table, run_name, normalize_by="median", run=False
    )
    assert "--norm-method median" in cmd
    assert "--control-sgrna" not in cmd
    cmd = exp.run_mageck(
        treats, controls, contrast_table, run_name, normalize_by="xxx", run=False
    )
    assert "--norm-method median" in cmd


def test_run_all_contrasts(analysis_test_data, tmpdir):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "test_run_all_contrasts"
    exp = Experiment(
        count_file,
        sample_file,
        controls,
        name,
        "Name",
        "day",
        "d0",
        "experiment",
        0.8,
        tmpdir,
    )
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    exp.write_control_barcodes_to_file()
    exp.run_all_contrasts()
    
    actual_gene_summary = tmpdir.join(
        "test_run_all_contrasts_d1_vs_d0.gene_summary.txt"
    )
    actual_results = pd.read_table(actual_gene_summary).to_dict()
    sample_results = {'id': {883: 'SL1344_3514', 371: 'glpG', 611: 'SL1344_3105', 408: 'sseJ', 773: 'SL1344_2152', 430: 'orgA', 573: 'nanA', 8: 'gor', 718: 'ompL', 606: 'SL1344_2547', 685: 'ybaL', 321: 'kdpD', 864: 'SL1344_0036', 461: 'SL1344_2327', 799: 'ydiF', 240: 'pilU', 254: 'fucI', 232: 'gcvT', 115: 'yafD', 431: 'ygcF'}, 'num': {883: 7, 371: 2, 611: 1, 408: 3, 773: 1, 430: 1, 573: 1, 8: 2, 718: 7, 606: 1, 685: 1, 321: 1, 864: 1, 461: 2, 799: 1, 240: 5, 254: 1, 232: 1, 115: 1, 431: 1}, 'neg|score': {883: 0.97665, 371: 0.44663, 611: 0.72382, 408: 0.49816, 773: 0.88056, 430: 0.51881, 573: 0.68871, 8: 0.0002556, 718: 0.83267, 606: 0.71818, 685: 0.79969, 321: 0.39969, 864: 0.95705, 461: 0.55207, 799: 0.90188, 240: 0.30161, 254: 0.31755, 232: 0.28119, 115: 0.14828, 431: 0.51944}, 'neg|p-value': {883: 1.0, 371: 0.88788, 611: 1.0, 408: 0.9641, 773: 1.0, 430: 0.83112, 573: 1.0, 8: 4.9515e-06, 718: 0.99963, 606: 1.0, 685: 1.0, 321: 0.83111, 864: 1.0, 461: 0.972, 799: 1.0, 240: 4.9515e-06, 254: 0.66584, 232: 0.66584, 115: 0.49998, 431: 0.83112}, 'neg|fdr': {883: 1.0, 371: 1.0, 611: 1.0, 408: 1.0, 773: 1.0, 430: 1.0, 573: 1.0, 8: 4.8e-05, 718: 1.0, 606: 1.0, 685: 1.0, 321: 1.0, 864: 1.0, 461: 1.0, 799: 1.0, 240: 4.8e-05, 254: 1.0, 232: 1.0, 115: 1.0, 431: 1.0}, 'neg|rank': {883: 884, 371: 372, 611: 612, 408: 409, 773: 774, 430: 431, 573: 574, 8: 9, 718: 719, 606: 607, 685: 686, 321: 322, 864: 865, 461: 462, 799: 800, 240: 241, 254: 255, 232: 233, 115: 116, 431: 432}, 'neg|goodsgrna': {883: 0, 371: 0, 611: 0, 408: 0, 773: 0, 430: 0, 573: 0, 8: 2, 718: 0, 606: 0, 685: 0, 321: 0, 864: 0, 461: 0, 799: 0, 240: 0, 254: 0, 232: 0, 115: 0, 431: 0}, 'neg|lfc': {883: 0.40284, 371: 0.2605, 611: 0.36349, 408: 0.11409, 773: 0.49383, 430: 0.24054, 573: 0.33439, 8: -2.0582, 718: 0.52408, 606: 0.35242, 685: 0.41556, 321: 0.1727, 864: 0.59549, 461: 0.15437, 799: 0.51573, 240: -0.20292, 254: 0.11926, 232: 0.088957, 115: -0.050174, 431: 0.24144}, 'pos|score': {883: 0.35092, 371: 0.28838, 611: 0.27618, 408: 0.96321, 773: 0.11944, 430: 0.48119, 573: 0.31129, 8: 0.99974, 718: 0.10669, 606: 0.28182, 685: 0.20031, 321: 0.60031, 864: 0.042947, 461: 0.82643, 799: 0.098119, 240: 0.99894, 254: 0.68245, 232: 0.71881, 115: 0.85172, 431: 0.48056}, 'pos|p-value': {883: 4.9515e-06, 371: 4.9515e-06, 611: 4.9515e-06, 408: 0.42228, 773: 4.9515e-06, 430: 0.16888, 573: 4.9515e-06, 8: 1.0, 718: 4.9515e-06, 606: 4.9515e-06, 685: 4.9515e-06, 321: 0.16889, 864: 4.9515e-06, 461: 0.30581, 799: 4.9515e-06, 240: 0.86703, 254: 0.33416, 232: 0.33416, 115: 0.50002, 431: 0.16888}, 'pos|fdr': {883: 1.1e-05, 371: 1.1e-05, 611: 1.1e-05, 408: 0.529577, 773: 1.1e-05, 430: 0.261893, 573: 1.1e-05, 8: 1.0, 718: 1.1e-05, 606: 1.1e-05, 685: 1.1e-05, 321: 0.261893, 864: 1.1e-05, 461: 0.43025, 799: 1.1e-05, 240: 0.953212, 254: 0.43025, 232: 0.43025, 115: 0.603184, 431: 0.261893}, 'pos|rank': {883: 308, 371: 251, 611: 240, 408: 846, 773: 104, 430: 420, 573: 281, 8: 904, 718: 90, 606: 246, 685: 177, 321: 525, 864: 28, 461: 717, 799: 86, 240: 898, 254: 594, 232: 621, 115: 739, 431: 419}, 'pos|goodsgrna': {883: 0, 371: 0, 611: 0, 408: 0, 773: 0, 430: 0, 573: 0, 8: 0, 718: 0, 606: 0, 685: 0, 321: 0, 864: 0, 461: 0, 799: 0, 240: 0, 254: 0, 232: 0, 115: 0, 431: 0}, 'pos|lfc': {883: 0.40284, 371: 0.2605, 611: 0.36349, 408: 0.11409, 773: 0.49383, 430: 0.24054, 573: 0.33439, 8: -2.0582, 718: 0.52408, 606: 0.35242, 685: 0.41556, 321: 0.1727, 864: 0.59549, 461: 0.15437, 799: 0.51573, 240: -0.20292, 254: 0.11926, 232: 0.088957, 115: -0.050174, 431: 0.24144}}
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])

    actual_gene_summary = tmpdir.join("test_run_all_contrasts_d2_vs_d0.gene_summary.txt")
    actual_results = pd.read_table(actual_gene_summary).to_dict()
    sample_results = {'id': {283: 'sopA', 512: 'exoX', 831: 'mgtC', 602: 'SL1344_3010', 757: 'tdcE', 294: 'pilT', 140: 'ybfA', 910: 'glpT', 469: 'SL1344_0502', 834: 'stjC', 158: 'SL1344_P1_0023', 903: 'yhjA', 898: 'ybdZ', 807: 'pagN', 748: 'yidZ', 850: 'SL1344_4132', 809: 'pipB2', 458: 'ydiM', 439: 'SL1344_1477', 219: 'dsdC'}, 'num': {283: 4, 512: 1, 831: 5, 602: 1, 757: 3, 294: 5, 140: 2, 910: 1, 469: 1, 834: 3, 158: 1, 903: 6, 898: 1, 807: 7, 748: 1, 850: 1, 809: 6, 458: 3, 439: 1, 219: 1}, 'neg|score': {283: 0.33083, 512: 0.58929, 831: 0.92508, 602: 0.69643, 757: 0.85132, 294: 0.34276, 140: 0.16717, 910: 0.99843, 469: 0.52788, 834: 0.92552, 158: 0.18578, 903: 0.9944, 898: 0.99091, 807: 0.90367, 748: 0.84305, 850: 0.9433, 809: 0.90425, 458: 0.51552, 439: 0.49154, 219: 0.26222}, 'neg|p-value': {283: 0.93663, 512: 1.0, 831: 1.0, 602: 1.0, 757: 1.0, 294: 0.96898, 140: 0.74754, 910: 1.0, 469: 1.0, 834: 1.0, 158: 0.8311, 903: 1.0, 898: 1.0, 807: 1.0, 748: 1.0, 850: 1.0, 809: 1.0, 458: 0.99555, 439: 1.0, 219: 1.0}, 'neg|fdr': {283: 1.0, 512: 1.0, 831: 1.0, 602: 1.0, 757: 1.0, 294: 1.0, 140: 1.0, 910: 1.0, 469: 1.0, 834: 1.0, 158: 1.0, 903: 1.0, 898: 1.0, 807: 1.0, 748: 1.0, 850: 1.0, 809: 1.0, 458: 1.0, 439: 1.0, 219: 1.0}, 'neg|rank': {283: 284, 512: 513, 831: 832, 602: 603, 757: 758, 294: 295, 140: 141, 910: 911, 469: 470, 834: 835, 158: 159, 903: 904, 898: 899, 807: 808, 748: 749, 850: 851, 809: 810, 458: 459, 439: 440, 219: 220}, 'neg|goodsgrna': {283: 0, 512: 0, 831: 0, 602: 0, 757: 0, 294: 0, 140: 0, 910: 0, 469: 0, 834: 0, 158: 0, 903: 0, 898: 0, 807: 0, 748: 0, 850: 0, 809: 0, 458: 0, 439: 0, 219: 0}, 'neg|lfc': {283: 0.43051, 512: 0.47724, 831: 0.76612, 602: 0.54777, 757: 0.46814, 294: 0.0043855, 140: -0.023426, 910: 2.1313, 469: 0.43301, 834: 0.556, 158: 0.15318, 903: 0.69436, 898: 1.2097, 807: 0.31998, 748: 0.68321, 850: 0.82663, 809: 0.44236, 458: 0.44779, 439: 0.40869, 219: 0.23214}, 'pos|score': {283: 0.65801, 512: 0.41071, 831: 0.0078054, 602: 0.30357, 757: 0.42715, 294: 0.99579, 140: 0.98872, 910: 0.0015664, 469: 0.47212, 834: 0.11814, 158: 0.81422, 903: 0.0039557, 898: 0.0090852, 807: 0.93267, 748: 0.15695, 850: 0.056704, 809: 0.94591, 458: 0.62638, 439: 0.50846, 219: 0.73778}, 'pos|p-value': {283: 4.9461e-06, 512: 4.9461e-06, 831: 4.9461e-06, 602: 4.9461e-06, 757: 4.9461e-06, 294: 4.9461e-06, 140: 0.55609, 910: 4.9461e-06, 469: 4.9461e-06, 834: 4.9461e-06, 158: 0.1689, 903: 4.9461e-06, 898: 4.9461e-06, 807: 4.9461e-06, 748: 4.9461e-06, 850: 4.9461e-06, 809: 4.9461e-06, 458: 4.9461e-06, 439: 4.9461e-06, 219: 4.9461e-06}, 'pos|fdr': {283: 6e-06, 512: 6e-06, 831: 6e-06, 602: 6e-06, 757: 6e-06, 294: 6e-06, 140: 0.606948, 910: 6e-06, 469: 6e-06, 834: 6e-06, 158: 0.197736, 903: 6e-06, 898: 6e-06, 807: 6e-06, 748: 6e-06, 850: 6e-06, 809: 6e-06, 458: 6e-06, 439: 6e-06, 219: 6e-06}, 'pos|rank': {283: 562, 512: 346, 831: 19, 602: 257, 757: 361, 294: 902, 140: 893, 910: 11, 469: 395, 834: 108, 158: 695, 903: 13, 898: 21, 807: 811, 748: 144, 850: 53, 809: 824, 458: 535, 439: 426, 219: 624}, 'pos|goodsgrna': {283: 0, 512: 0, 831: 3, 602: 0, 757: 0, 294: 0, 140: 0, 910: 1, 469: 0, 834: 1, 158: 0, 903: 3, 898: 1, 807: 0, 748: 0, 850: 1, 809: 0, 458: 0, 439: 0, 219: 0}, 'pos|lfc': {283: 0.43051, 512: 0.47724, 831: 0.76612, 602: 0.54777, 757: 0.46814, 294: 0.0043855, 140: -0.023426, 910: 2.1313, 469: 0.43301, 834: 0.556, 158: 0.15318, 903: 0.69436, 898: 1.2097, 807: 0.31998, 748: 0.68321, 850: 0.82663, 809: 0.44236, 458: 0.44779, 439: 0.40869, 219: 0.23214}}
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])
    
    actual_gene_summary = tmpdir.join("test_run_all_contrasts_d3_vs_d0.gene_summary.txt")
    actual_results = pd.read_table(actual_gene_summary).to_dict()
    sample_results =  {'id': {896: 'SL1344_3055', 817: 'basR', 630: 'SL1344_2698', 95: 'SL1344_2881', 726: 'SL1344_1929', 621: 'yicL', 350: 'SL1344_2537', 825: 'ugtL', 659: 'SL1344_1598', 842: 'nagZ', 774: 'yncC', 272: 'dcoC', 668: 'hsdS', 261: 'yafB', 312: 'shfB', 552: 'sciK1', 247: 'STnc440', 459: 'SL1344_2599', 297: 'SL1344_4152', 301: 'pilQ'}, 'num': {896: 2, 817: 1, 630: 5, 95: 1, 726: 5, 621: 5, 350: 1, 825: 2, 659: 1, 842: 1, 774: 2, 272: 1, 668: 1, 261: 3, 312: 1, 552: 1, 247: 1, 459: 1, 297: 1, 301: 2}, 'neg|score': {896: 0.98324, 817: 0.91129, 630: 0.74508, 95: 0.097492, 726: 0.83713, 621: 0.73765, 350: 0.44357, 825: 0.91916, 659: 0.77586, 842: 0.93386, 774: 0.87783, 272: 0.34514, 668: 0.78088, 261: 0.3317, 312: 0.38966, 552: 0.6768, 247: 0.31317, 459: 0.57461, 297: 0.37273, 301: 0.37545}, 'neg|p-value': {896: 0.97247, 817: 1.0, 630: 0.87006, 95: 0.3335, 726: 0.87006, 621: 0.87006, 350: 0.33351, 825: 0.88787, 659: 0.66584, 842: 1.0, 774: 0.74988, 272: 0.33351, 668: 0.66584, 261: 0.70277, 312: 0.33351, 552: 0.49936, 247: 0.33351, 459: 0.49935, 297: 0.33351, 301: 0.55465}, 'neg|fdr': {896: 1.0, 817: 1.0, 630: 0.983194, 95: 0.924097, 726: 0.983194, 621: 0.983194, 350: 0.924097, 825: 0.984393, 659: 0.965085, 842: 1.0, 774: 0.965085, 272: 0.924097, 668: 0.965085, 261: 0.965085, 312: 0.924097, 552: 0.924097, 247: 0.924097, 459: 0.924097, 297: 0.924097, 301: 0.924097}, 'neg|rank': {896: 897, 817: 818, 630: 631, 95: 96, 726: 727, 621: 622, 350: 351, 825: 826, 659: 660, 842: 843, 774: 775, 272: 273, 668: 669, 261: 262, 312: 313, 552: 553, 247: 248, 459: 460, 297: 298, 301: 302}, 'neg|goodsgrna': {896: 0, 817: 0, 630: 0, 95: 1, 726: 0, 621: 0, 350: 0, 825: 0, 659: 0, 842: 0, 774: 0, 272: 0, 668: 0, 261: 1, 312: 0, 552: 0, 247: 0, 459: 0, 297: 0, 301: 0}, 'neg|lfc': {896: 0.43425, 817: 0.50836, 630: -0.17994, 95: -0.79279, 726: -0.070438, 621: -0.31213, 350: -0.060968, 825: 0.4395, 659: 0.27565, 842: 0.59899, 774: 0.29301, 272: -0.18307, 668: 0.27613, 261: -0.41956, 312: -0.13056, 552: 0.17306, 247: -0.22076, 459: 0.067423, 297: -0.14802, 301: -0.28818}, 'pos|score': {896: 0.21019, 817: 0.088715, 630: 0.88743, 95: 0.90251, 726: 0.15135, 621: 0.84077, 350: 0.55643, 825: 0.088843, 659: 0.22414, 842: 0.066144, 774: 0.23012, 272: 0.65486, 668: 0.21912, 261: 0.98061, 312: 0.61034, 552: 0.3232, 247: 0.68683, 459: 0.42539, 297: 0.62727, 301: 0.87783}, 'pos|p-value': {896: 0.30683, 817: 4.9515e-06, 630: 0.96848, 95: 0.6665, 726: 4.9515e-06, 621: 0.96848, 350: 0.66649, 825: 4.9515e-06, 659: 0.33416, 842: 4.9515e-06, 774: 0.30683, 272: 0.66649, 668: 0.33416, 261: 0.96415, 312: 0.66649, 552: 0.50064, 247: 0.66649, 459: 0.50065, 297: 0.66649, 301: 0.88989}, 'pos|fdr': {896: 0.888145, 817: 3.6e-05, 630: 1.0, 95: 0.888145, 726: 3.6e-05, 621: 1.0, 350: 0.888145, 825: 3.6e-05, 659: 0.888145, 842: 3.6e-05, 774: 0.888145, 272: 0.888145, 668: 0.888145, 261: 1.0, 312: 0.888145, 552: 0.888145, 247: 0.888145, 459: 0.888145, 297: 0.888145, 301: 0.999905}, 'pos|rank': {896: 202, 817: 76, 630: 780, 95: 795, 726: 147, 621: 734, 350: 502, 825: 77, 659: 220, 842: 57, 774: 223, 272: 581, 668: 214, 261: 869, 312: 543, 552: 305, 247: 611, 459: 387, 297: 558, 301: 767}, 'pos|goodsgrna': {896: 0, 817: 0, 630: 0, 95: 0, 726: 1, 621: 0, 350: 0, 825: 0, 659: 0, 842: 0, 774: 0, 272: 0, 668: 0, 261: 0, 312: 0, 552: 0, 247: 0, 459: 0, 297: 0, 301: 0}, 'pos|lfc': {896: 0.43425, 817: 0.50836, 630: -0.17994, 95: -0.79279, 726: -0.070438, 621: -0.31213, 350: -0.060968, 825: 0.4395, 659: 0.27565, 842: 0.59899, 774: 0.29301, 272: -0.18307, 668: 0.27613, 261: -0.41956, 312: -0.13056, 552: 0.17306, 247: -0.22076, 459: 0.067423, 297: -0.14802, 301: -0.28818}}
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])
    
    actual_gene_summary = tmpdir.join("test_run_all_contrasts_d4_vs_d0.gene_summary.txt")
    actual_results = pd.read_table(actual_gene_summary).to_dict()
    sample_results = {'id': {258: 'SL1344_2249', 718: 'xapB', 424: 'cybC', 720: 'ssaL', 525: 'SL1344_2702', 893: 'SL1344_3229', 431: 'iolE', 852: 'SL1344_3704', 610: 'ygjT', 312: 'yicL', 536: 'mppA', 698: 'ycdZ', 659: 'yihX', 30: 'hypD', 564: 'ybbM', 499: 'SL1344_3250', 173: 'ychM', 198: 'yeeX', 509: 'uxaC', 894: 'araH'}, 'num': {258: 1, 718: 2, 424: 1, 720: 2, 525: 1, 893: 1, 431: 2, 852: 1, 610: 1, 312: 5, 536: 1, 698: 1, 659: 2, 30: 1, 564: 2, 499: 1, 173: 1, 198: 1, 509: 2, 894: 3}, 'neg|score': {258: 0.24859, 718: 0.81851, 424: 0.46489, 720: 0.82011, 525: 0.60282, 893: 0.98025, 431: 0.47517, 852: 0.94577, 610: 0.70313, 312: 0.32233, 536: 0.61661, 698: 0.79843, 659: 0.74874, 30: 0.011599, 564: 0.64786, 499: 0.56897, 173: 0.14013, 198: 0.17022, 509: 0.58177, 894: 0.98245}, 'neg|p-value': {258: 4.9515e-06, 718: 0.30432, 424: 0.16647, 720: 0.30432, 525: 0.16648, 893: 1.0, 431: 0.0281, 852: 0.83112, 610: 0.33175, 312: 0.034675, 536: 0.16648, 698: 0.49827, 659: 0.30431, 30: 4.9515e-06, 564: 0.30431, 499: 0.16648, 173: 4.9515e-06, 198: 4.9515e-06, 509: 0.0281, 894: 0.70421}, 'neg|fdr': {258: 1.7e-05, 718: 0.420103, 424: 0.257293, 720: 0.420103, 525: 0.257293, 893: 1.0, 431: 0.067883, 852: 0.856305, 610: 0.425362, 312: 0.083329, 536: 0.257293, 698: 0.593285, 659: 0.420103, 30: 1.7e-05, 564: 0.420103, 499: 0.257293, 173: 1.7e-05, 198: 1.7e-05, 509: 0.067883, 894: 0.810111}, 'neg|rank': {258: 259, 718: 719, 424: 425, 720: 721, 525: 526, 893: 894, 431: 432, 852: 853, 610: 611, 312: 313, 536: 537, 698: 699, 659: 660, 30: 31, 564: 565, 499: 500, 173: 174, 198: 199, 509: 510, 894: 895}, 'neg|goodsgrna': {258: 1, 718: 0, 424: 0, 720: 0, 525: 0, 893: 0, 431: 1, 852: 0, 610: 0, 312: 3, 536: 0, 698: 0, 659: 0, 30: 1, 564: 1, 499: 0, 173: 1, 198: 1, 509: 1, 894: 0}, 'neg|lfc': {258: -0.98939, 718: -0.41126, 424: -0.61852, 720: -0.35446, 525: -0.40634, 893: 1.8311, 431: -0.56822, 852: 0.79273, 610: -0.24157, 312: -0.71199, 536: -0.38241, 698: -0.048877, 659: -0.13055, 30: -5.4726, 564: -0.51684, 499: -0.46533, 173: -1.3567, 198: -1.2339, 509: -0.43278, 894: 0.14925}, 'pos|score': {258: 0.75141, 718: 0.61746, 424: 0.53511, 720: 0.517, 525: 0.39718, 893: 0.019749, 431: 0.47335, 852: 0.054232, 610: 0.29687, 312: 0.44169, 536: 0.38339, 698: 0.20157, 659: 0.19788, 30: 0.9884, 564: 0.57526, 499: 0.43103, 173: 0.85987, 198: 0.82978, 509: 0.38827, 894: 0.063497}, 'pos|p-value': {258: 1.0, 718: 0.9719, 424: 0.83353, 720: 0.9719, 525: 0.83352, 893: 4.9515e-06, 431: 0.89032, 852: 0.16888, 610: 0.66825, 312: 0.59773, 536: 0.83352, 698: 0.50173, 659: 0.30581, 30: 1.0, 564: 0.9719, 499: 0.83352, 173: 1.0, 198: 1.0, 509: 0.89031, 894: 0.07396}, 'pos|fdr': {258: 1.0, 718: 1.0, 424: 1.0, 720: 1.0, 525: 1.0, 893: 0.000216, 431: 1.0, 852: 1.0, 610: 1.0, 312: 1.0, 536: 1.0, 698: 1.0, 659: 1.0, 30: 1.0, 564: 1.0, 499: 1.0, 173: 1.0, 198: 1.0, 509: 1.0, 894: 1.0}, 'pos|rank': {258: 656, 718: 550, 424: 473, 720: 458, 525: 359, 893: 19, 431: 421, 852: 48, 610: 271, 312: 398, 536: 347, 698: 187, 659: 183, 30: 886, 564: 506, 499: 387, 173: 741, 198: 716, 509: 352, 894: 58}, 'pos|goodsgrna': {258: 0, 718: 0, 424: 0, 720: 0, 525: 0, 893: 1, 431: 0, 852: 1, 610: 0, 312: 0, 536: 0, 698: 0, 659: 0, 30: 0, 564: 0, 499: 0, 173: 0, 198: 0, 509: 0, 894: 1}, 'pos|lfc': {258: -0.98939, 718: -0.41126, 424: -0.61852, 720: -0.35446, 525: -0.40634, 893: 1.8311, 431: -0.56822, 852: 0.79273, 610: -0.24157, 312: -0.71199, 536: -0.38241, 698: -0.048877, 659: -0.13055, 30: -5.4726, 564: -0.51684, 499: -0.46533, 173: -1.3567, 198: -1.2339, 509: -0.43278, 894: 0.14925}}
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])

def test_process_results(analysis_test_data, tmpdir):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "test_process_results"
    exp = Experiment(
        count_file,
        sample_file,
        controls,
        name,
        "Name",
        "day",
        "d0",
        "experiment",
        0.8,
        tmpdir,
    )
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    exp.write_control_barcodes_to_file()
    exp.run_all_contrasts()
    exp.process_results()

    actual_rra = tmpdir.join("test_process_results_rra_results.csv")
    actual_results = pd.read_csv(actual_rra).to_dict()
    sample_results = {'Name': {1283: 'sciX', 1955: 'SL1344_4468', 2115: 'xylA', 1704: 'yieM', 3475: 'basR', 920: 'tatA', 2017: 'dcoA', 945: 'lrhA', 2361: 'SL1344_2691', 213: 'SL1344_1264', 1356: 'SL1344_1477', 1931: 'SL1344_1567', 1742: 'STM3026', 2251: 'SL1344_2738', 3285: 'ampH', 2865: 'prgH', 21: 'rfbX', 993: 'pspF', 3558: 'SL1344_0019', 99: 'ilvI'}, 'number_of_barcodes': {1283: 1, 1955: 5, 2115: 1, 1704: 1, 3475: 1, 920: 2, 2017: 1, 945: 1, 2361: 1, 213: 1, 1356: 1, 1931: 4, 1742: 1, 2251: 2, 3285: 1, 2865: 1, 21: 1, 993: 1, 3558: 3, 99: 2}, 'LFC': {1283: 0.35757, 1955: -0.57837, 2115: -0.16625, 1704: 0.72942, 3475: 0.05004, 920: -6.9198, 2017: -0.38311, 945: -2.4561, 2361: 0.14491, 213: 0.078651, 1356: 0.40869, 1931: 0.38251, 1742: 0.7928, 2251: 0.15916, 3285: -0.38851, 2865: -4.9191, 21: -2.7876, 993: -0.13799, 3558: -0.41019, 99: -0.118}, 'neg_selection_fdr': {1283: 1.0, 1955: 0.924097, 2115: 0.924097, 1704: 1.0, 3475: 0.593285, 920: 6.1e-05, 2017: 0.924097, 945: 6.1e-05, 2361: 0.924097, 213: 1.0, 1356: 1.0, 1931: 0.924097, 1742: 1.0, 2251: 0.924097, 3285: 0.257293, 2865: 1.7e-05, 21: 4.8e-05, 993: 1.0, 3558: 0.521566, 99: 4.8e-05}, 'pos_selection_fdr': {1283: 6e-06, 1955: 1.0, 2115: 0.888145, 1704: 6e-06, 3475: 1.0, 920: 1.0, 2017: 0.888145, 945: 1.0, 2361: 0.888145, 213: 0.43025, 1356: 6e-06, 1931: 0.888145, 1742: 6e-06, 2251: 3.6e-05, 3285: 1.0, 2865: 1.0, 21: 1.0, 993: 0.550381, 3558: 1.0, 99: 0.646139}, 'contrast': {1283: 'd2', 1955: 'd3', 2115: 'd3', 1704: 'd2', 3475: 'd4', 920: 'd2', 2017: 'd3', 945: 'd2', 2361: 'd3', 213: 'd1', 1356: 'd2', 1931: 'd3', 1742: 'd2', 2251: 'd3', 3285: 'd4', 2865: 'd4', 21: 'd1', 993: 'd2', 3558: 'd4', 99: 'd1'}}
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])

# def test_run_experiment(analysis_test_data, tmpdir, dnaid1315_expected_outcomes):
#     _, _, _, controls, count_file, sample_file, _,_ = analysis_test_data
#     name = "test_run_experiment"
#     exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
#     exp.run_experiment()
#     expected_rra = dnaid1315_expected_outcomes / "test_process_results_rra_results_no_index.csv"
#     actual_rra = tmpdir.join("test_run_experiment_rra_results.csv")
#     assert_files_are_same(actual_rra, expected_rra)

def test_run_experiment_no_batch(
    analysis_test_data, tmpdir
):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "test_run_experiment"
    exp = Experiment(
        count_file, sample_file, controls, name, "Name", "day", "d0", "", 0.8, tmpdir
    )
    exp.run_experiment()
    
    sample_results = {'Name': {1283: 'sciX', 1955: 'SL1344_4468', 2115: 'xylA', 1704: 'yieM', 3475: 'basR', 920: 'tatA', 2017: 'dcoA', 945: 'lrhA', 2361: 'SL1344_2691', 213: 'SL1344_1264', 1356: 'SL1344_1477', 1931: 'SL1344_1567', 1742: 'STM3026', 2251: 'SL1344_2738', 3285: 'ampH', 2865: 'prgH', 21: 'rfbX', 993: 'pspF', 3558: 'SL1344_0019', 99: 'ilvI'}, 'number_of_barcodes': {1283: 1, 1955: 5, 2115: 1, 1704: 1, 3475: 1, 920: 2, 2017: 1, 945: 1, 2361: 1, 213: 1, 1356: 1, 1931: 4, 1742: 1, 2251: 2, 3285: 1, 2865: 1, 21: 1, 993: 1, 3558: 3, 99: 2}, 'LFC': {1283: 0.35757, 1955: -0.57837, 2115: -0.16625, 1704: 0.72942, 3475: 0.05004, 920: -6.9198, 2017: -0.38311, 945: -2.4561, 2361: 0.14491, 213: 0.078651, 1356: 0.40869, 1931: 0.38251, 1742: 0.7928, 2251: 0.15916, 3285: -0.38851, 2865: -4.9191, 21: -2.7876, 993: -0.13799, 3558: -0.41019, 99: -0.118}, 'neg_selection_fdr': {1283: 1.0, 1955: 0.924097, 2115: 0.924097, 1704: 1.0, 3475: 0.593285, 920: 6.1e-05, 2017: 0.924097, 945: 6.1e-05, 2361: 0.924097, 213: 1.0, 1356: 1.0, 1931: 0.924097, 1742: 1.0, 2251: 0.924097, 3285: 0.257293, 2865: 1.7e-05, 21: 4.8e-05, 993: 1.0, 3558: 0.521566, 99: 4.8e-05}, 'pos_selection_fdr': {1283: 6e-06, 1955: 1.0, 2115: 0.888145, 1704: 6e-06, 3475: 1.0, 920: 1.0, 2017: 0.888145, 945: 1.0, 2361: 0.888145, 213: 0.43025, 1356: 6e-06, 1931: 0.888145, 1742: 6e-06, 2251: 3.6e-05, 3285: 1.0, 2865: 1.0, 21: 1.0, 993: 0.550381, 3558: 1.0, 99: 0.646139}, 'contrast': {1283: 'd2', 1955: 'd3', 2115: 'd3', 1704: 'd2', 3475: 'd4', 920: 'd2', 2017: 'd3', 945: 'd2', 2361: 'd3', 213: 'd1', 1356: 'd2', 1931: 'd3', 1742: 'd2', 2251: 'd3', 3285: 'd4', 2865: 'd4', 21: 'd1', 993: 'd2', 3558: 'd4', 99: 'd1'}}

    actual_results = pd.read_csv(tmpdir.join("test_run_experiment_rra_results.csv"))
    for key in sample_results.keys():
        ar = actual_results[key]
        assert all([v == ar[k] for k, v in sample_results[key].items()])


def test_run_experiment_no_control(
    analysis_test_data, tmpdir
):
    _, _, _, controls, count_file, sample_file, _, _ = analysis_test_data
    name = "test_run_experiment"
    exp = Experiment(
        count_file, sample_file, "", name, "Name", "day", "d0", "", 0.8, tmpdir
    )
    exp.run_experiment()
    actual_rra = tmpdir.join("test_run_experiment_rra_results.csv")
    actual_results = pd.read_csv(actual_rra).to_dict()
    sample_results = {'Name': {107: 'gtgA', 342: 'sseI', 2271: 'aroG', 633: 'yjfM', 1506: 'SL1344_3106', 3406: 'SL1344_0330', 2937: 'SL1344_0832', 1326: 'SL1344_0699', 77: 'sipA', 1209: 'tehB', 2692: 'STnc780', 3074: 'SL1344_3750', 2037: 'adi', 811: 'yibK', 1394: 'yrbD', 1734: 'STnc710', 1431: 'SL1344_0702', 71: 'yjgF', 2450: 'ydjM', 843: 'sopD'}, 'number_of_barcodes': {107: 3, 342: 1, 2271: 1, 633: 2, 1506: 2, 3406: 1, 2937: 3, 1326: 3, 77: 4, 1209: 1, 2692: 2, 3074: 3, 2037: 1, 811: 1, 1394: 1, 1734: 2, 1431: 2, 71: 1, 2450: 1, 843: 8}, 'LFC': {107: -0.026237, 342: -0.040027, 2271: -0.012038, 633: 0.008158, 1506: 0.048663, 3406: -3.7883, 2937: -5.5681, 1326: -0.081906, 77: -0.084469, 1209: -0.25546, 2692: 0.25638, 3074: -3.9149, 2037: -0.26044, 811: 0.27488, 1394: -0.032308, 1734: 0.15721, 1431: 0.035906, 71: -0.97829, 2450: 0.15414, 843: 0.021523}, 'neg_selection_fdr': {107: 0.994505, 342: 0.999867, 2271: 0.999995, 633: 0.999867, 1506: 0.999995, 3406: 0.998009, 2937: 0.998009, 1326: 0.999995, 77: 0.915541, 1209: 0.950881, 2692: 0.999995, 3074: 0.998009, 2037: 0.843819, 811: 0.999867, 1394: 0.99656, 1734: 0.999995, 1431: 0.999995, 71: 0.49404, 2450: 0.999995, 843: 0.999867}, 'pos_selection_fdr': {107: 0.999995, 342: 0.999995, 2271: 0.996152, 633: 0.999995, 1506: 0.999995, 3406: 0.970383, 2937: 0.970383, 1326: 0.999995, 77: 0.999995, 1209: 0.999995, 2692: 0.973244, 3074: 0.987142, 2037: 0.999995, 811: 0.999995, 1394: 0.999995, 1734: 0.999995, 1431: 0.999995, 71: 0.999995, 2450: 0.973244, 843: 0.999995}, 'contrast': {107: 'd1', 342: 'd1', 2271: 'd3', 633: 'd1', 1506: 'd2', 3406: 'd4', 2937: 'd4', 1326: 'd2', 77: 'd1', 1209: 'd2', 2692: 'd3', 3074: 'd4', 2037: 'd3', 811: 'd1', 1394: 'd2', 1734: 'd2', 1431: 'd2', 71: 'd1', 2450: 'd3', 843: 'd1'}}
    for key in sample_results.keys():
            ar = actual_results[key]
            assert all([v == ar[k] for k, v in sample_results[key].items()])