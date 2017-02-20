#!/bin/bash
source 00.sh

set -e

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Profiling raw tRNA sequences"
trna-profile $files/raw_tRNA_sequences.fa -o $output_dir/test_tRNA_profile.db -s test_sample

INFO "Printing stats"
trna-get-db-info -p $output_dir/test_tRNA_profile.db

INFO "Gen anti-codon profile"
trna-gen-anticodon-profile $output_dir/test_tRNA_profile.db -o $output_dir/anti-codon-profile.txt --percent-norm

INFO "Get tRNA sequences"
trna-get-sequences -p $output_dir/test_tRNA_profile.db -o $output_dir/tRNA_sequences.fa
