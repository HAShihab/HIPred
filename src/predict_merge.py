#!/usr/bin/python -u

import pandas

#
if __name__ == '__main__':
    
    Petrovski   = {} # RVIS
    Khurana     = {} # Indispensability Scores
    Rackham     = {} # EvoTol
    Huang_NoImp = {} # HIS
    Huang_Imp   = {} # HIS (Imputed)
    Steinberg   = {} # GHIS
    HIPred      = {}
    
    # Petrovski
    rm = []
    with open("./data/Petrovski.txt", "r") as fin:
        fin.next()
        for record in fin:
            record = record.strip().split("\t")
            
            if Petrovski.has_key(record[0]) and not Petrovski[record[0]] == record[1]:
                rm.append(record[0]); continue
            Petrovski[record[0]] = record[1]
    
    for x in set(rm):
        del Petrovski[x]
    
    # Khurana 
    rm = []
    with open("./data/Khurana.txt", "r") as fin:
        fin.next()
        for record in fin:
            record = record.strip().split("\t")
            
            if Khurana.has_key(record[0]) and not Khurana[record[0]] == record[-1]:
                rm.append(record[0]); continue
            Khurana[record[0]] = record[-1]
    
    for x in set(rm):
        del Khurana[x]
        
    # Rackham 
    rm = []
    with open("./data/Rackham.txt", "r") as fin:
        fin.next()
        for record in fin:
            record = record.strip().split("\t")
            
            if not record[1] == "1":
                break
            
            if Rackham.has_key(record[2]) and not Rackham[record[2]] == record[3]:
                rm.append(record[2]); continue
            Rackham[record[2]] = record[3]
            
    for x in set(rm):
        del Rackham[x]
        
    # Huang_NoImp
    rm = []
    with open("./data/Huang_NoImp.txt", "r") as fin:
        fin.next()
        for record in fin:
            record  = record.strip().split("\t")
            g, s, p = record[3].split("|")
            
            if Huang_NoImp.has_key(g) and not Huang_NoImp[g] == s:
                rm.append(g); continue
            Huang_NoImp[g] = s
    
    for x in set(rm):
        del Huang_NoImp[x]
    
    # Huang_Imp
    rm = []
    with open("./data/Huang_Imp.txt", "r") as fin:
        fin.next()
        for record in fin:
            record  = record.strip().split("\t")
            g, s, p = record[3].split("|")
            
            if Huang_Imp.has_key(g) and not Huang_Imp[g] == s:
                rm.append(g); continue
            Huang_Imp[g] = s
    
    for x in set(rm):
        del Huang_Imp[x]
    
    # Steinberg
    rm = []
    with open("./data/Steinberg.txt", "r") as fin:
        fin.next()
        for record in fin:
            record = record.strip().split("\t")
            
            if Steinberg.has_key(record[1]) and not Steinberg[record[1]] == record[2]:
                rm.append(record[1]); continue
            Steinberg[record[1]] = record[2]
    
    for x in set(rm):
        del Steinberg[x]
    
    # HIPred 
    rm = []
    with open("../HIPred.tsv", "r") as fin:
        fin.next()
        for record in fin:
            record = record.strip().split("\t")
            
            if HIPred.has_key(record[0]) and not HIPred[record[0]] == record[1]:
                rm.append(record[0]); continue
            HIPred[record[0]] = record[1]
    
    for x in set(rm):
        del HIPred[x]
    
    #
    
    keys = set(Petrovski.keys() + Khurana.keys() + Rackham.keys() + Huang_NoImp.keys() + Huang_Imp.keys() + Steinberg.keys() + HIPred.keys())
    df   = pandas.DataFrame(index=keys, columns=['RVIS', 'IS', 'EvoTol', 'HIS', 'HIS - Imputed', 'GHIS', 'HIPred'])
    
    for x in keys:
        if Petrovski.has_key(x):
            df.ix[x]['RVIS']   = float(Petrovski[x])
        if Khurana.has_key(x):
            df.ix[x]['IS']     = float(Khurana[x])
        if Rackham.has_key(x):
            df.ix[x]['EvoTol'] = float(Rackham[x])
        if Huang_NoImp.has_key(x):
            df.ix[x]['HIS']    = float(Huang_NoImp[x])
        if Huang_Imp.has_key(x):
            df.ix[x]['HIS - Imputed'] = float(Huang_Imp[x])
        if Steinberg.has_key(x):
            df.ix[x]['GHIS']   = float(Steinberg[x])
        if HIPred.has_key(x):
            df.ix[x]['HIPred'] = float(HIPred[x])
    
    df.to_csv('./tmp/prediction_matrix.tsv', sep="\t")
    