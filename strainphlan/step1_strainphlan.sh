#!/bin/bash
metaphlan2.py --input_type fastq <(zcat $1 $2) --samout $3
sample2markers.py --ifn_samples $3 --input_type sam --output_dir $4
