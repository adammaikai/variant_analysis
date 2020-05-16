#!/bin/bash

# Run VAPr to annotate and output MAF

python3 vapr_annotate.py \
	--vcf_dir vcf_dir \
	--output_dir output_dir \
	--output_maf maf \
	--mongodb_name my_mongodb \
	--collection_name collection1 \
	--annovar_path /path/to/annovar \
	--remove_collection False \
	--gzipped False \
	--threads n_threads
