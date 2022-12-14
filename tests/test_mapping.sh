# Test that mapping reports aligning read pairs for a contrived example.
function test_mapping {
    # Create two read pairs using lines from the human genome fasta.
    r1_1=$TEMPDIR/data_files/PCMP_stub_human_R1.fastq.gz
    r2_1=$TEMPDIR/data_files/PCMP_stub_human_R2.fastq.gz
    r1_2=$TEMPDIR/data_files/PCMP_stub2_human_R1.fastq.gz
    r2_2=$TEMPDIR/data_files/PCMP_stub2_human_R2.fastq.gz
    human=$TEMPDIR/indexes/human.fasta
    (
        echo "@read0"
        sed -n 2p $human
        echo "+"
        sed -n 's:.:G:g;2p' $human
    ) | gzip > $r1_1
    (
        echo "@read0"
        sed -n 2p $human | rev | tr '[ACTG]' '[TGAC]'
        echo "+"
        sed -n 's:.:G:g;2p' $human
    ) | gzip > $r2_1
    (
        echo "@read0"
        sed -n 15p $human | cut -c 1-40
        echo "+"
        sed -n 's:.:G:g;2p' $human | cut -c 1-40
    ) | gzip > $r1_2
    (
        echo "@read0"
        sed -n 15p $human | cut -c 1-40 | rev | tr '[ACTG]' '[TGAC]'
        echo "+"
        sed -n 's:.:G:g;2p' $human | cut -c 1-40
    ) | gzip > $r2_2
    # Run sunbeam mapping rules with these two samples defined.
    (
	    echo "stub_human,$r1_1,$r2_1"
	    echo "stub2_human,$r1_2,$r2_2"
    ) > $TEMPDIR/samples_test_mapping.csv
    sunbeam config modify --str 'all: {samplelist_fp: "samples_test_mapping.csv"}' \
        $TEMPDIR/tmp_config.yml > $TEMPDIR/test_mapping_config.yml

    # Move human host files to top-level, since we're using that genome for
    # mapping in this test and shouldn't decontaminate using it as well.
    for file in $TEMPDIR/hosts/human*; do
        mv $file $TEMPDIR/hosts_${file##*/}
    done
    sunbeam run --configfile $TEMPDIR/test_mapping_config.yml all_mapping
    # Move human host files back to original location
    for file in $TEMPDIR/hosts_*; do
        mv $file ${file/hosts_/hosts\//}
    done
    # After the header line, there should be two lines in the human and phix
    # coverage summaries, with two reads mapping for human and none for phix.
    # The lines should be sorted in standard alphanumeric order; stub2_human
    # will come before stub_human.
    (
	    csv_human=$TEMPDIR/sunbeam_output/mapping/human/coverage.csv
	    csv_phix=$TEMPDIR/sunbeam_output/mapping/phix174/coverage.csv
	    function col3 { cut -f3 -d, | tr '\n' : ; }
	    function col5 { cut -f5 -d, | tr '\n' : ; }
	    test "Sample:stub2_human:stub_human:" == $(col3 < "$csv_human")
	    test "Max:2:2:" == $(col5 < "$csv_human")
	    test "Sample:stub2_human:stub_human:" == $(col3 < "$csv_phix")
	    test "Max:0:0:" == $(col5 < "$csv_phix")
    ) || (
	    echo "Unexpected coverage.csv content from mapping rules" > /dev/stderr
        cat $csv_human
        cat $csv_phix
	    false
	)

    echo "test_mapping passed" >> test_results
}