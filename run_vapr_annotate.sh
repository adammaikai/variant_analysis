#!/bin/bash

# Run VAPr to annotate and output MAF

python3 vapr_annotate.py \
	--vcf_dir /mnt/data1/adam/jamieson/jane/apobec_dna/vcf/vapr_test \
	--output_dir /mnt/data1/adam/jamieson/jane/apobec_dna/maf/vapr_test_out \
	--output_maf Isquith_vapr_test.maf \
	--mongodb_name Isquith \
	--collection_name test \
	--annovar_path /mnt/data1/adam/software/annovar \
	--remove_collection False \
	--gzipped False \
	--threads 12
