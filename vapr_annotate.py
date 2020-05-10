# VAPr Demo

from VAPr import vapr_core
from VAPr.formatting import maf_formatter, create_whole_dataset_list
from pymongo import MongoClient
import pandas as pd
import os
import sys
import time
import argparse
import re

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='Annotate VCF file/s with VAPr and output as MAF file.')
parser.add_argument('--vcf_dir', required=True, help='Directory where VCF files live.')
parser.add_argument('--output_dir', required=True, help='Output directory.')
parser.add_argument('--output_maf', required=True, help='Output MAF name (file name only.')
parser.add_argument('--mongodb_name', required=True, help='MongoDB name. Suggested: PI name.')
parser.add_argument('--collection_name', required=True, help='MongoDB collection name. May be date, version, or experiment name.')
parser.add_argument('--annovar_path', required=True, help='Path to Annovar installation.')
parser.add_argument('--remove_collection', type=str2bool, nargs='?', const=True, default=False, help='Must be one of True/False. If True, remove MongoDB collection after printing MAF file.')
parser.add_argument('--gzipped', type=str2bool, nargs='?', const=True, default=False, help='Must be one of True/False. Are VCF files gzipped?')
parser.add_argument('--threads', default=1, help='Number of threads to parallelize annotation.')


# command line arguments
args = parser.parse_args()
print("\n\nInput Arguments:")
for arg in vars(args):
    print(arg, ":", getattr(args, arg))
print("\n\n")

vcf_dir = args.vcf_dir
output_dir = args.output_dir
output_maf = args.output_maf
mongodb_name = args.mongodb_name
collection_name = args.collection_name
annovar_path = args.annovar_path
remove_collection = args.remove_collection
gzipped=args.gzipped
threads = int(args.threads)

# Instantiate VapAnnotator Object
annotator = vapr_core.VaprAnnotator(input_dir=vcf_dir,
										output_dir=output_dir,
										mongo_db_name=mongodb_name,
									    mongo_collection_name=collection_name,
									    build_ver='hg19',
									    vcfs_gzipped=gzipped,
									    annovar_install_path=annovar_path)

# VAPr annotate
dataset = annotator.annotate(num_processes=threads)
ds = create_whole_dataset_list(mongodb_name, collection_name)

# # Convert and write to MAF format
maf = maf_formatter(ds2)
maf.to_csv(os.path.join(output_dir, output_maf), index=False, sep="\t", chunks=10000)

if remove_collection:
	client = MongoClient()
	db = getattr(client, mongodb_name)
	collection = getattr(db, collection_name)
	collection.remove()


