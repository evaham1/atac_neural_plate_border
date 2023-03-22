#!/usr/bin/env python

# load libraries
import os
import argparse
import sys
import pandas as pd
import numpy as np

import scanpy as sc
#import pyranges as pr
import warnings

def parse_args(args=None):
    Description = " "
    Epilog = " "

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help="Input path.", metavar='')
    parser.add_argument('-c', '--ncores', type=str, help="Number of CPUs.", metavar='')
    return parser.parse_args(args)

def main(args=None):

    # read in command line arguments
    args = parse_args(args)

    # set output paths
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)
    if not os.path.exists(rds_path):
        os.mkdir(rds_path)

    #####   1) Cell-level metadata (includes seacell assignments)
    # Load data
    ad_full = sc.read(args.input + '/AnnData_metacells_assigned.h5ad')
    print(ad_full)
    # Export data
    ad_full.obs.to_csv(os.path.join(rds_path,'Cell_metadata.csv'))

    #####   2) Feature-level metadata
    # Export data
    ad_full.var.to_csv(os.path.join(rds_path,'Feature_metadata.csv'))

    #####   3) Summarised count matrix
    # Load data
    ad_sum = sc.read(args.input + '/AnnData_summarised_by_metacells.h5ad')
    print(ad_sum)
    # Export data
    ad_sum.to_df(layer="raw").to_csv(os.path.join(rds_path,"Summarised_by_metacells_counts.csv"))

if __name__ == '__main__':
    sys.exit(main())