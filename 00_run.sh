#!/bin/bash

# for testing purposes get the first 10,000 lines from both R1 and R2
gzip -cd ../TP-10/TP-10_TTAGGC_L003_R1_001.fastq.gz | head -n 10000 > r1
gzip -cd ../TP-10/TP-10_TTAGGC_L003_R2_001.fastq.gz | head -n 10000 > r2

# generate a test config ini file from samples.txt
iu-gen-configs samples.txt

# merge r1 and r2, retain only the overlapping parts, and look for both directions (`--marker-gene-stringent`)
iu-merge-pairs test.ini --retain-only-overlap --marker-gene-stringent

# filter all the reads with more than 0 mismatch to minimize the impact of sequencing errors:
iu-filter-merged-reads test_MERGED -m 0 -o test_MERGED_FINAL
