#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Profiling raw tRNA sequences"
trna-profile $files/raw_tRNA_sequences.fa -o $output_dir/test_tRNA_profile.db -n Test
