# coding: utf8
# -*- coding: utf-8 -*-
#
#   Copyright 2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import os, glob, urllib2, multiprocessing
from cinfony import pybel, rdk
pybel.ob.obErrorLog.SetOutputLevel(-1)

URL = "http://dude.docking.org/targets/%s/%s_final.sdf.gz"

targets = ['AA2AR', 'ABL1', 'ACE', 'ACES', 'ADA', 'ADA17', 'ADRB1', 'ADRB2', 'AKT1', 'AKT2', 'ALDR', 'AMPC', 'ANDR', 'AOFB', 'BACE1', 'BRAF', 'CAH2', 'CASP3', 'CDK2', 'COMT', 'CP2C9', 'CP3A4', 'CSF1R', 'CXCR4', 'DEF', 'DHI1', 'DPP4', 'DRD3', 'DYR', 'EGFR', 'ESR1', 'ESR2', 'FA10', 'FA7', 'FABP4', 'FAK1', 'FGFR1', 'FKB1A', 'FNTA', 'FPPS', 'GCR', 'GLCM', 'GRIA2', 'GRIK1', 'HDAC2', 'HDAC8', 'HIVINT', 'HIVPR', 'HIVRT', 'HMDH', 'HS90A', 'HXK4', 'IGF1R', 'INHA', 'ITAL', 'JAK2', 'KIF11', 'KIT', 'KITH', 'KPCB', 'LCK', 'LKHA4', 'MAPK2', 'MCR', 'MET', 'MK01', 'MK10', 'MK14', 'MMP13', 'MP2K1', 'NOS1', 'NRAM', 'PA2GA', 'PARP1', 'PDE5A', 'PGH1', 'PGH2', 'PLK1', 'PNPH', 'PPARA', 'PPARD', 'PPARG', 'PRGR', 'PTN1', 'PUR2', 'PYGM', 'PYRD', 'RENI', 'ROCK1', 'RXRA', 'SAHH', 'SRC', 'TGFR1', 'THB', 'THRB', 'TRY1', 'TRYB1', 'TYSY', 'UROK', 'VGFR2', 'WEE1', 'XIAP']

def prepare(target):
        for typ in ('actives', 'decoys'):
            url = URL % (target.lower(), typ)
            d = os.path.join("Targets_test", target, typ)
            if not os.path.isdir(d):
                os.makedirs(d)
            molfn = os.path.join(d, os.path.basename(url))
            if not os.path.isfile(molfn):
                print "Downloading %s" % url
                i = urllib2.urlopen(url)
                o = open(molfn, 'wb')
                o.write(i.read())
                i.close()
                o.close()
            #########Only 7 RB #########
            if not os.path.isfile(molfn.replace('.sdf.gz', '_7_RB.sdf')):
                o = pybel.Outputfile('sdf', molfn.replace('.sdf.gz', '_7_RB.sdf'))
                n = 0
                for mol in pybel.readfile('sdf', molfn):
                    if mol.OBMol.NumRotors() <= 7:
                        o.write(mol)
                        n += 1
                o.close()
                if not n:
                    try:
                        os.remove(molfn.replace('.sdf.gz', '_7_RB.sdf'))
                    except:
                        pass
            #########################
            clean(validate(molfn, typ))
        return target
#            print "Splitting %s" % molfn
#            for mol in pybel.readfile('sdf', molfn):
#                if mol.OBMol.NumRotors() > 7: continue
#                mol.write('sdf', os.path.join(d, mol.title + '_final.sd'), overwrite=True)

def validate(fn, typ):
    ofn = fn.replace('.sdf.gz', '_sanitized.sdf')
    if not os.path.isfile(ofn):
        good = 0
        total = 0
        o = rdk.Outputfile('sdf', ofn, overwrite=True)
        n = 0
        for mol in pybel.readfile('sdf', fn):
            mol.title = typ + "_" + mol.title + "_" + str(n)
            total +=1
            n +=1
            try:
                o.write(rdk.readstring('mol', mol.write('mol')))
                good +=1
            except Exception, e:
                print e
                pass
        print "%s %% good mols (%s/%s)" % ((good*100./total), good, total)
        o.close()
    return ofn

def clean(fn):
    ligands = os.path.join(os.path.dirname(fn), "..", "PDB","ligands_allgood.sdf")
    if not os.path.isfile(ligands):
        return
    inchikeys = set([mol.write('inchikey') for mol in pybel.readfile('sdf', ligands)])
    ofn = fn.replace('.sdf', '_filtered.sdf')
    if not os.path.isfile(ofn):
        good = 0
        total = 0
        o = rdk.Outputfile('sdf', ofn, overwrite=True)
        for mol in pybel.readfile('sdf', fn):
            total +=1
            if mol.write('inchikey') not in inchikeys:
                try:
                    o.write(rdk.readstring('mol', mol.write('mol')))
                    good +=1
                except Exception, e:
                    #print e
                    pass
        o.close()
        print "%s %% remaining mols (%s/%s)" % ((good*100./total), good, total)


pool = multiprocessing.Pool()
for t in pool.imap(prepare, targets):
    print t
#for target in targets:
#    t = prepare(target)
#    print t
