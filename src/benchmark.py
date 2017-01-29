#!/usr/bin/python -u

import os
import math
import numpy
import random
import pandas

import matplotlib.pyplot as plt

from sklearn import *

#
if __name__ == '__main__':
    
    numpy.random.seed(1)
    random.seed(1)   
    
    predmatrix = pandas.read_csv("./tmp/prediction_matrix.tsv", sep="\t", header=0, index_col=0)
    
    
    # benchmark the algorithms using the training data ...
    print 
    print "performance on training data"
    print
    
    Train = {}
    rm    = []
    
    with open("./data/Dang.tsv", "r") as fin:
        fin.next()
        for record in fin:
            record = record.strip().split("\t")
            
            Train[record[0]] = 1
    
    with open("./data/1kg_LoFT.tsv", "r") as fin:
        for record in fin:
            record = record.strip()
            
            if Train.has_key(record):
                rm.append(record); continue
            Train[record] = -1
    
    print "  # Pos (raw): {0}".format(Train.values().count(1))
    print "  # Neg (raw): {0}".format(Train.values().count(-1))
    
    # ignore inconsistent data
    print "found {0} inconsistent record(s) ...".format(len(rm))
    
    for x in rm:
        print "  {0}".format(x)
        
        del Train[x]
    
    print
    print "  # Pos (processed): {0}".format(Train.values().count(1))
    print "  # Neg (processed): {0}".format(Train.values().count(-1))
    print
    
    
    for fn in predmatrix.columns:
        # ignore our performance (report cv statistic instead)
        if fn == "HIPred":
            continue 
            
        df = predmatrix[fn]
        df = df.ix[ Train.keys() ]
        
        p = []
        q = []
        for x in Train:
            if not pandas.isnull(df[x]):
                p.append(Train[x])
                q.append(df[x])
        
        if fn in [ "Petrovski", "Rackham" ]:
            # reverse score
            q = numpy.array(q) * -1

        if fn in [ "Petrovski", "Rackham" ]:
             pred   = [ 0 if x < 0.0 else 1 for x in q ]
        else:
             pred   = [ 0 if x < 0.5 else 1 for x in q ]

        matrix = pandas.crosstab(numpy.array(p), numpy.array(pred), rownames=['Truth'], colnames=['Predicted'], margins=True)
        
        print matrix

        AUC = metrics.roc_auc_score( numpy.array(p), numpy.array(q) )
        
        print fn
        print "   n             : {0}".format(len(q))
        print "  AUC            : {0}".format(round(AUC, 4))
        print

        tp = float(matrix[1][1])
        fp = float(matrix[1][0])
        tn = float(matrix[0][0])
        fn = float(matrix[0][1])

        Accuracy    = (tp + tn) / (tp + fn + tn + fn)
        Sensitivity = tp / (tp + fn)
        Specificity = tn / (tn + fp)
        PPV         = tp / (tp + fp)
        NPV         = tn / (tn + fn)

        print
        print "  tp             : {0}".format(tp)
        print "  fp             : {0}".format(fp)
        print "  tn             : {0}".format(tn)
        print "  fn             : {0}".format(fn)
        print
        print "  Accuracy       : {0}".format(round(Accuracy, 4))
        print "  Sensitivity    : {0}".format(round(Sensitivity, 4))
        print "  Specificity    : {0}".format(round(Specificity, 4))
        print "  Precision (PPV): {0}".format(round(PPV, 4))
        print "  NPV            : {0}".format(round(NPV, 4))
        print "  AUC            : {0}".format(round(AUC, 4))
        print
    
    
    
    # now, perform the benchmark using an independent test set ... 
    print 
    print "performance on validation data"
    print
    
    benchmark = pandas.read_csv("./data/Petrovski_DatasetS1.csv", sep="\t", header=0, index_col=0)
    
    # match negative instances based on gene length ...
    genesize  = {}
    
    with open("./data/Annotations/meta.txt", "r") as fin:
        fin.next()
        for record in fin:
            record = record.strip().split(",")
            
            genesize[record[0]] = math.fabs( int(record[4]) - int(record[3]) )
    
    for dataset in [ "OMIM Haploinsufficiency", "OMIM de novo & Haploinsuficciency", "MGI Lethality orthologs", "MGI Seizure orthologs", 'ASD1', 'ASD2' ]:

        if dataset == "ASD1":
            positive = [ x.strip() for x in open("./data/ASD1.txt", "r") if len(x) > 0 ]
        elif dataset == "ASD2":
            positive = [ x.strip() for x in open("./data/ASD2.txt", "r") if len(x) > 0 ]
        else:
            positive = list(benchmark[ benchmark[dataset] > 0 ].index)
        negative = []
        
        
        # remove training data from positive dataset
        for x in Train:
            if x in positive:
                positive.remove(x)
        
        # match each positive with a negative based on gene size
        for x in positive:
            if genesize.has_key(x):
                p = genesize[x]
                q = dict(genesize)
                
                del q[x] # remove exact match so it isn't chosen
                
                while True:
                    x = min(q, key=lambda x:abs(q[x] - p))
                    
                    # append to negative if not in positive (just in case) ...
                    if not x in positive:
                        negative.append( x ); break
                    del q[x]
                    
        # remove ambiguous data and assign labels
        bmark = {}
        for q in set(negative) - set(positive):
            bmark[q] = -1
        for q in set(positive):
            bmark[q] = 1
        
        
        print dataset
        for x in bmark:
            print "{0},{1}".format(x, bmark[x])
        print
        print
        
        if dataset == "OMIM Haploinsufficiency":
            dataset = "OMIM HI"
        if dataset == "OMIM de novo & Haploinsuficciency":
            dataset = "OMIM HI de novo"
        if dataset == "MGI Lethality orthologs":
            dataset = "MGI Lethality"
        if dataset == "MGI Seizure orthologs":
            dataset = "MGI Seizure"
            
        print
        print "# dataset :", dataset
        print "# Pos: {0}".format(bmark.values().count(1))
        print "# Neg: {0}".format(bmark.values().count(-1))
        print
        
        perf  = {}
        for method in predmatrix.columns:
            df = predmatrix[method]
            df = df.ix[ bmark.keys() ]
            
            p = []
            q = []
            for x in bmark:
                if not pandas.isnull(df[x]):
                    p.append(bmark[x])
                    q.append(df[x])
            
            if method in [ "RVIS", "EvoTol" ]:
                # reverse order ...
                q = numpy.array(q) * -1
            
            fpr, tpr, _ = metrics.roc_curve( numpy.array(p), numpy.array(q) )
            auc         = metrics.auc(fpr, tpr)
            
            perf[method] = { 'auc': auc, 'tpr': tpr, 'fpr': fpr }
            
            
            if method in [ "RVIS", "EvoTol" ]:
                pred   = [ 0 if x < 0.0 else 1 for x in q ]
            else:
                pred   = [ 0 if x < 0.5 else 1 for x in q ]

            matrix = pandas.crosstab(numpy.array(p), numpy.array(pred), rownames=['Truth'], colnames=['Predicted'], margins=True)

            tp = float(matrix[1][1])
            fp = float(matrix[1][0])
            tn = float(matrix[0][0])
            fn = float(matrix[0][1])

            Accuracy    = "%.4f" % ((tp + tn) / (tp + fn + tn + fn))
            Sensitivity = "%.4f" % (tp / (tp + fn))
            Specificity = "%.4f" % (tn / (tn + fp))
            PPV         = "%.4f" % (tp / (tp + fp))
            NPV         = "%.4f" % (tn / (tn + fn))
            
            print "     " + " & ".join(map(str, [ method, Accuracy, Sensitivity, Specificity, PPV, NPV, "%.4f" % auc ]))

        
        fig   = plt.figure()
        plt.plot([0, 1], [0, 1], 'k--')
        
        for x in sorted(perf.items(), key=lambda x: x[1]['auc']):
            plt.plot(x[1]['fpr'], x[1]['tpr'], label="{0} (AUC:{1})".format(x[0], "%.2f" % x[1]['auc']))
        
        plt.xlim([0.0, 1.05])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(dataset)
        plt.legend(loc="lower right")
        
        fig.savefig("./tmp/{0}.png".format(dataset), dpi=fig.dpi)
    
    print
    