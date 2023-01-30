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

    #####   1) export for the full anndata object, need cell_ids:seacell_ids dictionary
    # Load data
    ad_full = sc.read('AnnData_metacells_assigned.h5ad')
    print(ad_full)
    # Export data
    ad.obs.to_csv('AnnData_metacells_assigned_cell_metadata.csv')

    #####   2) export for the summarised anndata object, need summarised peak count matrix
    ad_sum = sc.read('AnnData_summarised_by_metacells.h5ad')
    print(ad_sum)
    # Export data

if __name__ == '__main__':
    sys.exit(main())