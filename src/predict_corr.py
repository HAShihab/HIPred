#!/usr/bin/python -u

import pandas

#
if __name__ == '__main__':
    
    df   = pandas.read_csv("./tmp/prediction_matrix.tsv", sep="\t", header=0, index_col=0)
    
    data = df.corr(method="spearman")
    data.to_csv('./tmp/correlation_matrix.tsv', sep="\t", float_format='%.4f')
    