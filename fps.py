# coding: utf8
# -*- coding: utf-8 -*-
#
#   Copyright 2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import os, csv, sys, itertools, numpy, glob, pickle, multiprocessing
from pyroc import *
from rdkit.ML.Scoring import Scoring
fptypes = (
    "OpenEye-MACCS166",
    "OpenEye-Path",
    "OpenEye-Circular",
    "OpenEye-Tree",
    "RDKit-MACCS166",
    "RDKit-Fingerprint",
    "RDKit-Morgan",
    "RDKit-Torsion",
    "RDKit-AtomPair",
    "OpenBabel-FP2/1",
    "OpenBabel-FP3",
    "OpenBabel-FP4",
    "OpenBabel-MACCS",
    "ChemFP-Substruct-OpenEye/1",
    "ChemFP-Substruct-RDKit/1",
    "ChemFP-Substruct-OpenBabel/1",
    )
combo_list = (
              ("Combo-ag", [
                      "RDKit-Morgan"
                      , "OpenEye-Path"
                      , "RDKit-Torsion"
                      , "RDKit-AtomPair"
                      , "OpenBabel-FP2/1"
                      , "OpenEye-Tree"]
                  , 6)
              , ("Combo-ag2", [
                      "RDKit-Morgan"
                      , "OpenEye-Path"
                      , "RDKit-Torsion"
                      , "RDKit-AtomPair"
                      , "OpenBabel-FP2/1"
                      , "OpenEye-Tree"]
                  , 2)
              , ("Combo-free", [
                       "RDKit-Morgan"
                      , "RDKit-Torsion"
                      , "RDKit-AtomPair"
                      , "OpenBabel-FP2/1"]
                  , 4)
              , ("Combo-free2", [
                       "RDKit-Morgan"
                      , "RDKit-Torsion"
                      , "RDKit-AtomPair"
                      , "OpenBabel-FP2/1"]
                  , 2)
              , ("Combo-RDKit", [
                       "RDKit-Morgan"
                      , "RDKit-Torsion"
                      , "RDKit-AtomPair"]
                  , 3)
              , ("Combo-OE", [
                       "OpenEye-Circular"
                      , "OpenEye-Path"
                      , "OpenEye-Tree"]
                  , 3)
              , ("Combo-Morgan-Path", [
                       "RDKit-Morgan"
                      , "OpenEye-Path"]
                  , 2)
              , ("Combo-Torsion-Path", [
                       "RDKit-Torsion"
                      , "OpenEye-Path"]
                  , 2)
              , ("Combo-Torsion-FP2", [
                       "RDKit-Torsion"
                      , "OpenBabel-FP2/1"]
                  , 2)
              , ("Combo-Morgan-FP2", [
                       "RDKit-Morgan"
                      , "OpenBabel-FP2/1"]
                  , 2)
              , ("Combo-Torsion-Circular", [
                       "OpenEye-Circular"
                      , "RDKit-Torsion"]
                  , 2)
              , ("Combo-AtomPair-Tree", [
                       "OpenEye-Tree"
                      , "RDKit-AtomPair"]
                  , 2)
              , ("Combo-AtomPair-Path", [
                       "OpenEye-Path"
                      , "RDKit-AtomPair"]
                  , 2)
              )
#combo_list = []#TODO
INFIX = "_allgood"


def parsereport(report, n=5):#Per validació i prou!!!!!!!!!
    r = open(report)
    result = ['', ''] + [0]*n
    header = []
    for line in csv.reader(r, delimiter="\t"):
        if 'Name' not in line:
            for n in xrange(len(result)):
                if not result[n] and n > 1:
                    v = float(line[n].strip())
                    if v > result[n]:
                        result[n] = v
                elif n <= 1:
                    result[n] = line[n]
        else:
            header = line
    r.close()
    result.pop(2)
    header.pop(2)
    return result, header

def build_results(target):
    results_dict = {}
    target_dir = os.path.join(root, target)
    print target_dir
    for fptype in  fptypes:
        csvfn = "".join(c for c in fptype if c.isalnum() or c in (' ','.','_')).rstrip() + '_new.csv'
        print csvfn
        try:
            nlines = 0
            for line in csv.reader(open(os.path.join(target_dir, csvfn))):
                nlines +=1
                try:
                    if len(results_dict[line[0]]) == len(fptypes):
                        continue
                    results_dict[line[0]].append(float(line[1]))
                except KeyError:
                    print line
                    results_dict[line[0]] = [float(line[1])]
            print nlines, len(results_dict)
        except IOError:
            return None
    print "Results done"
    return results_dict
def separe_methods(results_dict, colnames=[]):
    methods = []
    total = float(len(results_dict))
    n = 0
    for key,  values in results_dict.items():
        #values = [values[0]]
        values = values[:]
        if not methods:
            methods = [[] for v in values]
        for x in xrange(len(values)):
            score =  (total-n)/total
            methods[x].append([key, values[x]])
        n += 1
    for method in methods:
        method.sort(key=lambda row: row[1], reverse=True)
    return methods

def zscore(results_dict, nbest=6, colnames=[]):
    methods = separe_methods(results_dict, colnames=colnames)
    col_idx = []
    if colnames:
        col_idx = [fptypes.index(coln) for coln in colnames]
    if col_idx:
        methods = [m for m in methods if methods.index(m) in col_idx]
    #methods.sort(key=lambda method: method[0][1], reverse=True)
    zscore_dict = {}
    for method in methods:
        values = numpy.array([float(v[1]) for v in method])
        average = values.mean()
        std = values.std()
        for k, v in method:
            v = (v - average)/std
            try:
                zscore_dict[k].append(v)
                zscore_dict[k].sort(reverse=True)
            except:
                zscore_dict[k] = [v]
    for k,  zscores in zscore_dict.items():
        slice = zscores[:nbest]
        zscore_dict[k] = [float(sum(slice)/len(slice))]
    return zscore_dict

def parallel_selection(results_dict):
    methods = separe_methods(results_dict)
    methods.sort(key=lambda method: method[0][1], reverse=True)
    total = float(len(results_dict))
    methods = [iter(method) for method in methods]
    methods = itertools.cycle(methods)
    n = 0
    selected = {}
    for method in methods:
        if n>=total:
            break
        for key, value in method:
            if key not in selected:
                selected[key] = [(total-n)]
                n += 1
                break
    return selected


def for_pyroc(results_dict):
    pyroc_samples = []
    print "results: %s" % len(results_dict)
    for key, values in results_dict.items():
        #print values
        if "CHEMBL" in key or "ligands_" in key or "actives_" in key:
             active = 1
        else:
            active = 0
        if not pyroc_samples:
            for x in xrange(len(values)):
                pyroc_samples.append([])
        #print len(pyroc_samples)
        for x in xrange(len(values)):
            #print x
            score = values[x]
            pyroc_samples[x].append((active, score))
    return [ROCData(pyroc_sample) for pyroc_sample in pyroc_samples]

def get_EF(roc, percent, relative=False):
    """
    DUD-E version
    """
    if not relative: #Previous version
        percent = float(percent)
        total = len(roc.data)
        total_actives = len([v for v in roc.data if v[0]])
        found = 0
        #print "Total actives: %s out of %s" % (total_actives,  len(roc.data))
        stop = percent*total/100.0
        n = 1
        while n <= stop and n < total:
            if roc.data[n-1][0]:
                found +=1
            n +=1
        ef = (100.0*found/total_actives)
        if relative:
            maxfound = min(int(stop), total_actives)
            maxef = (100.0*maxfound/total_actives)
            print "EF:%s,found %s" % (ef, found)
            print "max EF:%s, max found %s\n" % (maxef, maxfound)
            return ef/maxef
        return ef
    else: #DUD-E version
        percent = float(percent)
        total = len(roc.data)
        print total
        total_actives = len([v for v in roc.data if v[0]])
        total_decoys = total - total_actives
        found = 0
        #print "Total actives: %s out of %s" % (total_actives,  len(roc.data))
        stop = percent*total_decoys/100.0
        print "%s%% of %s decoys is %s" % (percent,total_decoys, stop)
        n = 1
        decoys = 0
        while decoys <= stop and found < total_actives:
            if roc.data[n-1][0]:
                found +=1
            else:
                decoys += 1
            n +=1
        ef = (100.0*found/total_actives)
        print "found %s actives of a total of %s, EF%s = %s" % (found, total_actives, percent, ef)
        return ef


def plot_rocs(rocs, labels, target):
    print len(rocs)
    print len(labels)
    print rocs
    print labels
    labels = [labels[x]+ '; AUC: %s; EF1: %s' % (str(rocs[x].auc())[:4], str(get_EF(rocs[x], 1))[:4]) for x in xrange(len(rocs))]
    plot_multiple_roc(rocs,title=target, labels=labels,include_baseline=True, equal_aspect=True)

def build_rocs(target, plot=True, results_dict=None):
    labels = []
    if not results_dict and target:
        results_dict = build_results(target)
    if not results_dict: return None, None, None

    #Methods
    rocs = for_pyroc(results_dict)
    #ShapeTanimoto ColorTanimoto	EON_ET_pb	EON_ET_coul EON_ShapeTanimoto
    labels += list(fptypes)
    n = len(labels)
    #Parallel selection
    rocs += for_pyroc(parallel_selection(results_dict))
    labels.append("Parallel")
    #Zscore
    rocs += for_pyroc(zscore(results_dict, nbest=6))
    labels.append("ZScore")
    #Zscore x best
#    for x in xrange(1, n):
#        rocs += for_pyroc(zscore(results_dict, nbest=x))
#        labels.append("Z%s" % x)
    #Zscore combos
    for combo in combo_list:
        rocs += for_pyroc(zscore(results_dict, nbest=combo[2], colnames = combo[1]))
        labels.append(combo[0])
    if plot:
        return plot_rocs(rocs, labels, target)
    return (labels, rocs, results_dict)



if __name__=="__main__":
#    DB = "DUD"
#    DB = "DUD-E"
#    DB = "dekois"
#    DB = "common"
    DBs = ["DUD", "DUD-E", "dekois", "common"]
    DBs = ["ZINC-DUD-E", "ZINC-dekois"]
    DBs = ["ZINC-DUD-E"]#TODO
    cwd = os.getcwd()
    for DB in DBs:
        if DB == "DUD-E":
            actives = "actives"
            root = cwd + "/Targets_test"
            alltargets = ['XIAP','AA2AR', 'ABL1', 'ACE', 'ACES', 'ADA', 'ADA17', 'ADRB1', 'ADRB2', 'AKT1', 'AKT2', 'ALDR', 'AMPC', 'ANDR', 'AOFB', 'BACE1', 'BRAF', 'CAH2', 'CASP3', 'CDK2', 'COMT', 'CP2C9', 'CP3A4', 'CSF1R', 'CXCR4', 'DEF', 'DHI1', 'DPP4', 'DRD3', 'DYR', 'EGFR', 'ESR1', 'ESR2', 'FA10', 'FA7', 'FABP4', 'FAK1', 'FGFR1', 'FKB1A', 'FNTA', 'FPPS', 'GCR', 'GLCM', 'GRIA2', 'GRIK1', 'HDAC2', 'HDAC8', 'HIVINT', 'HIVPR', 'HIVRT', 'HMDH', 'HS90A', 'HXK4', 'IGF1R', 'INHA', 'ITAL', 'JAK2', 'KIF11', 'KIT', 'KITH', 'KPCB', 'LCK', 'LKHA4', 'MAPK2', 'MCR', 'MET', 'MK01', 'MK10', 'MK14', 'MMP13', 'MP2K1', 'NOS1', 'NRAM', 'PA2GA', 'PARP1', 'PDE5A', 'PGH1', 'PGH2', 'PLK1', 'PNPH', 'PPARA', 'PPARD', 'PPARG', 'PRGR', 'PTN1', 'PUR2', 'PYGM', 'PYRD', 'RENI', 'ROCK1', 'RXRA', 'SAHH', 'SRC', 'TGFR1', 'THB', 'THRB', 'TRY1', 'TRYB1', 'TYSY', 'UROK', 'VGFR2', 'WEE1']
            badtargets = ["FGFR1"]
        elif DB == "ZINC-DUD-E":
            actives = "actives"
            root = cwd + "/Targets_test/ZINC"
            alltargets = ['XIAP','AA2AR', 'ABL1', 'ACE', 'ACES', 'ADA', 'ADA17', 'ADRB1', 'ADRB2', 'AKT1', 'AKT2', 'ALDR', 'AMPC', 'ANDR', 'AOFB', 'BACE1', 'BRAF', 'CAH2', 'CASP3', 'CDK2', 'COMT', 'CP2C9', 'CP3A4', 'CSF1R', 'CXCR4', 'DEF', 'DHI1', 'DPP4', 'DRD3', 'DYR', 'EGFR', 'ESR1', 'ESR2', 'FA10', 'FA7', 'FABP4', 'FAK1', 'FGFR1', 'FKB1A', 'FNTA', 'FPPS', 'GCR', 'GLCM', 'GRIA2', 'GRIK1', 'HDAC2', 'HDAC8', 'HIVINT', 'HIVPR', 'HIVRT', 'HMDH', 'HS90A', 'HXK4', 'IGF1R', 'INHA', 'ITAL', 'JAK2', 'KIF11', 'KIT', 'KITH', 'KPCB', 'LCK', 'LKHA4', 'MAPK2', 'MCR', 'MET', 'MK01', 'MK10', 'MK14', 'MMP13', 'MP2K1', 'NOS1', 'NRAM', 'PA2GA', 'PARP1', 'PDE5A', 'PGH1', 'PGH2', 'PLK1', 'PNPH', 'PPARA', 'PPARD', 'PPARG', 'PRGR', 'PTN1', 'PUR2', 'PYGM', 'PYRD', 'RENI', 'ROCK1', 'RXRA', 'SAHH', 'SRC', 'TGFR1', 'THB', 'THRB', 'TRY1', 'TRYB1', 'TYSY', 'UROK', 'VGFR2', 'WEE1']
            badtargets = ["FGFR1"]
            alltargets = ["PTN1"]#TODO
        elif DB == "ZINC-dekois":
            actives = "actives"
            root = cwd + "/Targets_test/ZINC/dekois"
            alltargets = ['HDAC2', 'JNK1', 'THROMBIN', 'ACE', 'GBA', 'DHFR', 'HDAC8', 'TPA', 'ALR2', 'EGFR', 'SARS-HCoV', 'RXRa', 'A2A', 'ER-beta', 'KIF11', 'PPARa', 'ADAM17', 'VEGFR2', 'PYGL-in', 'PDE4b', 'JNK2', 'PR', 'MK2', 'JNK3', 'TP', 'COX2', 'NA', 'COX1', 'CDK2', 'QPCT', 'AURKA', 'AURKB', 'TS', 'P38-alpha', 'HIV1RT', '17betaHSD1', 'JAK3', 'ERBB2', 'SIRT2', 'IGF1R', 'CTSK', 'GR', 'EPHB4', 'uPA', 'HIV1PR', 'PARP-1', 'TK', 'SRC', 'TIE2', 'HSP90', 'AR', 'CYP2A6', 'PRKCQ', 'ACHE', 'MDM2', 'FXA', 'FGFR1', 'INHA', 'PNP', 'ITK', 'AKT1', 'PIM-1', 'PIM-2', 'PI3Kg', 'CATL', 'ROCK1', 'GSK3B', 'PDK1', 'PPARg', 'ADRB2', 'PYGL-out', 'BRAF', 'FKBP1A', '11betaHSD1', 'ACE2', 'PDE5', 'VEGFR1', 'MMP2', 'LCK', 'HMGR', 'BCL2']
            badtargets = ['TPA', 'PRKCQ', 'MMP2', 'HIV1PR', 'MDM2']
        elif DB == "DUD":
            actives = "ligands"
            root = cwd + "/Targets_test/DUD"
            alltargets = ['ace', 'ache', 'ada', 'alr2', 'ampc', 'ar', 'cdk2', 'comt', 'cox1', 'cox2', 'dhfr', 'egfr', 'er_agonist', 'er_antagonist', 'fgfr1', 'fxa', 'gart', 'gpb', 'gr', 'hivpr', 'hivrt', 'hmga', 'hsp90', 'inha', 'mr', 'na', 'p38', 'parp', 'pde5', 'pdgfrb', 'pnp', 'ppar_gamma', 'pr', 'rxr_alpha', 'sahh', 'src', 'thrombin', 'tk', 'trypsin', 'vegfr2']
            badtargets = ["er_antagonist"]
        elif DB == "dekois":
            actives = "actives"
            root = cwd + "/Targets_test/dekois"
            alltargets = ['HDAC2', 'JNK1', 'THROMBIN', 'ACE', 'GBA', 'DHFR', 'HDAC8', 'TPA', 'ALR2', 'EGFR', 'SARS-HCoV', 'RXRa', 'A2A', 'ER-beta', 'KIF11', 'PPARa', 'ADAM17', 'VEGFR2', 'PYGL-in', 'PDE4b', 'JNK2', 'PR', 'MK2', 'JNK3', 'TP', 'COX2', 'NA', 'COX1', 'CDK2', 'QPCT', 'AURKA', 'AURKB', 'TS', 'P38-alpha', 'HIV1RT', '17betaHSD1', 'JAK3', 'ERBB2', 'SIRT2', 'IGF1R', 'CTSK', 'GR', 'EPHB4', 'uPA', 'HIV1PR', 'PARP-1', 'TK', 'SRC', 'TIE2', 'HSP90', 'AR', 'CYP2A6', 'PRKCQ', 'ACHE', 'MDM2', 'FXA', 'FGFR1', 'INHA', 'PNP', 'ITK', 'AKT1', 'PIM-1', 'PIM-2', 'PI3Kg', 'CATL', 'ROCK1', 'GSK3B', 'PDK1', 'PPARg', 'ADRB2', 'PYGL-out', 'BRAF', 'FKBP1A', '11betaHSD1', 'ACE2', 'PDE5', 'VEGFR1', 'MMP2', 'LCK', 'HMGR', 'BCL2']
            badtargets = ['TPA', 'PRKCQ', 'MMP2']
        elif DB == "common":
            actives = "actives"
            alltargets = ['ACE', 'ACES', 'ANDR', 'CDK2', 'DYR', 'EGFR', 'FA10', 'GCR', 'HIVPR', 'HIVRT', 'HMDH', 'HS90A', 'INHA', 'KITH', 'MK14', 'NRAM', 'PARP1', 'PDE5A', 'PGH1', 'PGH2', 'PNPH', 'PPARG', 'PRGR', 'RXRA', 'SRC', 'THRB', 'VGFR2']
            badtargets = []
            root = cwd + "/Targets_test/common"
        try:
            targets = [sys.argv[1].upper()]
            writefile = False
        except IndexError:
            if DB not in ("common", "ZINC-DUD-E", "ZINC-dekois"):
                targets = [target for target in alltargets if  os.path.isfile(os.path.join(root,target,"PDB/ligands%s.sdf" % INFIX))]
            else:
                targets =  alltargets
            targets = [t for t in targets if t not in badtargets] #Target with no decoys? DUDE bug?
            print targets
            print root
            #targets = ['trypsin']
            writefile = True
        header = []
        auc_rows = []
        ef1_rows = []
        ef1r_rows = []
        bedroc_rows = []
        all_results_dict = {}
        ###
        def process_results(target):
            header = []
            print "###Checking target %s..." % target
            labels, rocs, results_dict = build_rocs(target, plot=not(writefile))
            if labels == rocs == None: return None
            print "Data loaded"
            if not header:
                header = labels
            auc_row = [target]
            ef1_row = [target]
            ef1r_row = [target]
            bedroc_row = [target]
            for roc in rocs:
                auc_row.append(roc.auc())
                ef1_row.append(get_EF(roc, 10, relative=True))
                ef1r_row.append(get_EF(roc, 1, relative=True))
                bedroc_row.append(Scoring.CalcBEDROC(roc.data, 0, 20))
            return header, auc_row, ef1_row, ef1r_row, bedroc_row, target

        pool = multiprocessing.Pool()
        for r in pool.imap(process_results, targets):
        #for r in (process_results(t) for t in targets):
            if not r: continue
            header, auc_row, ef1_row, ef1r_row, bedroc_row, target = r
            if not header:
                header = labels
            auc_rows.append(auc_row)
            ef1_rows.append(ef1_row)
            ef1r_rows.append(ef1r_row)
            bedroc_rows.append(bedroc_row)
        def build_average(n, n_rows):
            letters_orig = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            letters = list(letters_orig) + [letters_orig[0] + c for c in letters_orig]
            l = 0
            averagerow = ['Average:']
            for column in auc_rows[0][1:]:
                l +=1
                letter = letters[l]
                formula = "=AVERAGE(" + letter + str(n+1) + ":" +  letter + str(n + n_rows) + ")"
                averagerow.append(formula)
            return averagerow
        n = 1
        n_rows = len(auc_rows)

        res = open("resum_fps_%s.csv" % DB, 'w')#TODO
        w = csv.writer(res)
        w.writerow(['AUC'] + header)
        #
        w.writerows(auc_rows)
        averagerow = build_average(n, n_rows)
        n += n_rows +3
        w.writerow(averagerow)
        w.writerow([''])
        #
        w.writerow(['EF10'] + header)
        w.writerows(ef1_rows)
        averagerow = build_average(n, n_rows)
        n += n_rows +3
        w.writerow(averagerow)
        w.writerow([''])
        #
        w.writerow(['EF1 (DUD-E)'] + header)
        w.writerows(ef1r_rows)
        averagerow = build_average(n, n_rows)
        n += n_rows +3
        w.writerow(averagerow)
        w.writerow([''])
        #
        w.writerow(['BEDROC'] + header)
        w.writerows(bedroc_rows)
        averagerow = build_average(n, n_rows)
        n += n_rows +3
        w.writerow(averagerow)
        w.writerow([''])

        res.close()
