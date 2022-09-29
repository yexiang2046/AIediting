#!/usr/bin/env python3
#Small Script to annotate sailor output files with feature and geneID of called editing sites

import argparse
import sys
import os
import pandas as pd
import pybedtools as bedtools
from datetime import datetime

if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('--gtf', required=True, help='GTF file containing human annotations)
    ap.add_argument('--fwd', required=True, help='fwd bed file from sailor output')
    ap.add_argument('--rev', required=True, help='rev bed file from sailor output')
    ap.add_argument('--o', required=True, help='path for output file')
    args = ap.parse_args()

    startTime = datetime.now()

    #reads in files
    #merged = bedtools.BedTool(args.bed)
    gtf = bedtools.BedTool(args.gtf)
    fwd = bedtools.BedTool(args.fwd)
    rev = bedtools.BedTool(args.rev)
    gtf = bedtools.BedTool.sort(gtf)


    #merges the fwd and rev bed files
    merged = bedtools.BedTool.cat(fwd, rev, postmerge = False)
    merged = bedtools.BedTool.sort(merged)


    #peforms intersection
    merged_gtf = merged.closest(gtf, s = True, D = 'ref', loj = True)

    #converts to pandas dataframe
    df = bedtools.BedTool.to_dataframe(merged_gtf, header = None, names = ['chr', 'pos-1', 'pos', 'coverage', 'conf',
                                                                           'strand', 'chrGTF', 'source', 'feature',
                                                                           'pos_1_GTF', 'pos_2_GTF', 'dot', 'strandGTF',
                                                                           'dot2', 'info', 'distance'])

    #removes unneeded columns
    df_cut = df.drop(['chrGTF', 'source', 'pos_1_GTF', 'pos_2_GTF', 'dot', 'strandGTF', 'dot2'], axis = 1)
    df_cut['absValDist'] = df_cut['distance'].abs()

    #adds numerical value to feature types
    df_cut.loc[df_cut['feature'] == 'transcript', 'numerical'] = 0
    df_cut.loc[df_cut['feature'] == 'gene', 'numerical'] = 1
    df_cut.loc[df_cut['feature'] == 'exon', 'numerical'] = 2
    df_cut.loc[df_cut['feature'] == 'UTR', 'numerical'] = 3
    df_cut.loc[df_cut['feature'] == 'CDS', 'numerical'] = 4

    #sorts the dataframe by numerical for removing duplicates
    df_cut_sorted = df_cut.sort_values(by = ['absValDist', 'numerical'], ascending = False)

    #removes duplicates and keeps the highest value
    df_nodup = df_cut_sorted.drop_duplicates(subset = ['chr', 'pos-1', 'pos', 'coverage', 'conf', 'strand'], keep = 'first')

    #splits the info column
    df_nodup_expand = df_nodup['info'].str.split(';', expand = True)

    df_nodup_expand = df_nodup_expand.drop([5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24], axis=1)
    df_nodup_expand['geneID'] = ''
    df_nodup_expand['bioType'] = ''
    df_nodup_expand['geneName'] = ''

    for index, row in df_nodup_expand.iterrows():
        if row[2] != None and row[2].startswith(' gene_type'):
            df_nodup_expand.loc[[index], ['geneID', 'bioType', 'geneName']] = [row[0][9:-1], row[2][12:-1], row[3][12:-1]]
        elif row[1] != None and row[1].startswith(' gene_type'):
            df_nodup_expand.loc[[index], ['geneID', 'bioType', 'geneName']] = [row[0][9:-1], row[1][12:-1], row[2][12:-1]]
        else:
            pass

    df_test = df_nodup_expand

    #removes the unneeded info
    df_nodup_expand = df_nodup_expand.drop([0, 1, 2, 3, 4], axis = 1)

    #combines the two dataframes
    concat = pd.concat([df_nodup, df_nodup_expand], axis = 1)

    #removes all of the unnecessary stuff from the info column

    concat_clean = concat.drop(['info', 'numerical', 'absValDist'], axis = 1)
    concat_clean = concat_clean.replace({'.' : 'intergenic', '' : 'intergenic'})
    concat_clean = concat_clean.replace({'feature' : 'gene'}, 'intron')
    concat_clean['feature_distance'] = ''

    for index, row in concat_clean.iterrows():
        if row['distance'] == 0:
            concat_clean.loc[[index], ['feature_distance']] = 'Within a' + ' ' + row['bioType'] + ' ' + row['feature']

        elif 0 < row['distance'] < 2000:
            concat_clean.loc[[index], ['feature_distance']] = 'Upstream of a' + ' ' + row['bioType'] + ' ' + row['feature']

        elif 0 > row['distance'] > -2000:
            concat_clean.loc[[index], ['feature_distance']] = 'Downstream of a' + ' ' + row['bioType'] + ' ' + row['feature']

        elif row['distance'] < -2000 or row['distance'] > 2000:
            concat_clean.loc[[index], ['feature_distance']] = 'intergenic'

    concat_clean = concat_clean.drop(['feature'], axis = 1)

    print('\nTime elasped: ', datetime.now() - startTime)

    #writes the file to output csv file
    concat_clean.to_csv(args.o)