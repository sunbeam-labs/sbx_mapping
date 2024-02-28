import os


def test_full_run(run_sunbeam):
    (
        bfragilis_sliding_cov_fp,
        ecoli_sliding_cov_fp,
        bfragilis_filtered_cov_fp,
        ecoli_filtered_cov_fp,
        bfragilis_num_reads_fp,
        ecoli_num_reads_fp,
        benchmarks_fp,
    ) = run_sunbeam

    # Check output
    assert os.path.exists(bfragilis_sliding_cov_fp)
    assert os.path.exists(ecoli_sliding_cov_fp)
    assert os.path.exists(bfragilis_filtered_cov_fp)
    assert os.path.exists(ecoli_filtered_cov_fp)
    assert os.path.exists(bfragilis_num_reads_fp)
    assert os.path.exists(ecoli_num_reads_fp)

    assert os.stat(bfragilis_sliding_cov_fp).st_size != 0
    assert os.stat(ecoli_sliding_cov_fp).st_size != 0
