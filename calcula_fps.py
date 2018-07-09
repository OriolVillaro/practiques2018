# -*- coding: utf-8 -*-
#
#   Copyright 2013 Adrià Cereto Massagué <adrian.cereto@.urv.cat>
#
import chemfp, os, csv, multiprocessing

INFIX = "_allgood"
DB = "DUD"
#DB = "DUD-E"
#DB = "dekois"
#DB = "common"
pre_root = os.path.abspath(os.path.dirname(__file__))
DECOYS = None
#DECOYS = "ZINC"
decoystdir = pre_root + "/Targets_test/ZINC/"
if DB == "DUD-E":
    actives = "actives"
    root = pre_root + "/Targets_test"
    decoystdir = pre_root + "/Targets_test/ZINC/"
    alltargets = ['XIAP','AA2AR', 'ABL1', 'ACE', 'ACES', 'ADA', 'ADA17', 'ADRB1', 'ADRB2', 'AKT1', 'AKT2', 'ALDR', 'AMPC', 'ANDR', 'AOFB', 'BACE1', 'BRAF', 'CAH2', 'CASP3', 'CDK2', 'COMT', 'CP2C9', 'CP3A4', 'CSF1R', 'CXCR4', 'DEF', 'DHI1', 'DPP4', 'DRD3', 'DYR', 'EGFR', 'ESR1', 'ESR2', 'FA10', 'FA7', 'FABP4', 'FAK1', 'FGFR1', 'FKB1A', 'FNTA', 'FPPS', 'GCR', 'GLCM', 'GRIA2', 'GRIK1', 'HDAC2', 'HDAC8', 'HIVINT', 'HIVPR', 'HIVRT', 'HMDH', 'HS90A', 'HXK4', 'IGF1R', 'INHA', 'ITAL', 'JAK2', 'KIF11', 'KIT', 'KITH', 'KPCB', 'LCK', 'LKHA4', 'MAPK2', 'MCR', 'MET', 'MK01', 'MK10', 'MK14', 'MMP13', 'MP2K1', 'NOS1', 'NRAM', 'PA2GA', 'PARP1', 'PDE5A', 'PGH1', 'PGH2', 'PLK1', 'PNPH', 'PPARA', 'PPARD', 'PPARG', 'PRGR', 'PTN1', 'PUR2', 'PYGM', 'PYRD', 'RENI', 'ROCK1', 'RXRA', 'SAHH', 'SRC', 'TGFR1', 'THB', 'THRB', 'TRY1', 'TRYB1', 'TYSY', 'UROK', 'VGFR2', 'WEE1']
    badtargets = ["FGFR1"]
    #alltargets = ['DPP4']
elif DB == "DUD":
    actives = "ligands"
    root = pre_root + "/Targets_test/DUD"
    alltargets = ['ace', 'ache', 'ada', 'alr2', 'ampc', 'ar', 'cdk2', 'comt', 'cox1', 'cox2', 'dhfr', 'egfr', 'er_agonist', 'er_antagonist', 'fgfr1', 'fxa', 'gart', 'gpb', 'gr', 'hivpr', 'hivrt', 'hmga', 'hsp90', 'inha', 'mr', 'na', 'p38', 'parp', 'pde5', 'pdgfrb', 'pnp', 'ppar_gamma', 'pr', 'rxr_alpha', 'sahh', 'src', 'thrombin', 'tk', 'trypsin', 'vegfr2']
    badtargets = ["er_antagonist"]
elif DB == "dekois":
    actives = "actives"
    decoystdir = pre_root + "/Targets_test/ZINC/dekois"
    root = pre_root + "/Targets_test/dekois"
    alltargets = ['HDAC2', 'JNK1', 'THROMBIN', 'ACE', 'GBA', 'DHFR', 'HDAC8', 'TPA', 'ALR2', 'EGFR', 'SARS-HCoV', 'RXRa', 'A2A', 'ER-beta', 'KIF11', 'PPARa', 'ADAM17', 'VEGFR2', 'PYGL-in', 'PDE4b', 'JNK2', 'PR', 'MK2', 'JNK3', 'TP', 'COX2', 'NA', 'COX1', 'CDK2', 'QPCT', 'AURKA', 'AURKB', 'TS', 'P38-alpha', 'HIV1RT', '17betaHSD1', 'JAK3', 'ERBB2', 'SIRT2', 'IGF1R', 'CTSK', 'GR', 'EPHB4', 'uPA', 'HIV1PR', 'PARP-1', 'TK', 'SRC', 'TIE2', 'HSP90', 'AR', 'CYP2A6', 'PRKCQ', 'ACHE', 'MDM2', 'FXA', 'FGFR1', 'INHA', 'PNP', 'ITK', 'AKT1', 'PIM-1', 'PIM-2', 'PI3Kg', 'CATL', 'ROCK1', 'GSK3B', 'PDK1', 'PPARg', 'ADRB2', 'PYGL-out', 'BRAF', 'FKBP1A', '11betaHSD1', 'ACE2', 'PDE5', 'VEGFR1', 'MMP2', 'LCK', 'HMGR', 'BCL2']
    badtargets = ['TPA', 'PRKCQ', 'MMP2']
elif DB == "common":
    actives = "actives"
    alltargets = ['ACE', 'ACES', 'ANDR', 'CDK2', 'DYR', 'EGFR', 'FA10', 'GCR', 'HIVPR', 'HIVRT', 'HMDH', 'HS90A', 'INHA', 'KITH', 'MK14', 'NRAM', 'PARP1', 'PDE5A', 'PGH1', 'PGH2', 'PNPH', 'PPARG', 'PRGR', 'RXRA', 'SRC', 'THRB', 'VGFR2']
    badtargets = []
    root = pre_root + "/Targets_test/common"
fptypes = (
    "OpenEye-MACCS166",
    "OpenEye-Path",
    "OpenEye-Circular",
    "OpenEye-Tree",
    "RDKit-MACCS166",
    "RDKit-Fingerprint",
    "RDKit-Morgan",#Falla als decoys de XIAP
    "RDKit-Torsion",
    "RDKit-AtomPair",
    "OpenBabel-FP2/1",
    "OpenBabel-FP3",
    "OpenBabel-FP4",#Falla als decoys de XIAP
    "OpenBabel-MACCS",
    "ChemFP-Substruct-OpenEye/1",
    "ChemFP-Substruct-RDKit/1",
    "ChemFP-Substruct-OpenBabel/1",
    )

def fingerprint_screen(target):
    target_dir = os.path.join(root, target)
#    activesfn = os.path.join(target_dir, "actives", "actives_final_7_RB.sdf")
#    decoysfn = os.path.join(target_dir, "decoys", "decoys_final_7_RB.sdf")
    if DB == "DUD":
        t = target + "_"
        f = ""
    elif DB == "DUD-E":
        t = ""
        f = "_final"
    elif DB in ("dekois", "common"):
        t = ""
        f = ""
    activesfn = os.path.join(target_dir, actives, t + actives + f +"_sanitized_filtered.sdf")
    if DECOYS == "ZINC":
        decoysfn = os.path.join(decoystdir, target, "decoys", "decoys.sdf")
    else:
        decoysfn = os.path.join(target_dir, "decoys", t + "decoys" + f + "_sanitized_filtered.sdf")

    if not (os.path.isfile(activesfn) and os.path.isfile(decoysfn)):
        print activesfn
        return "Missing actives or decoys for target %s" % target

    for fptype in fptypes:
        fpname =  "".join(c for c in fptype if c.isalnum() or c in (' ','.','_')).rstrip()
        def load_fps(strucfile):
            fp_fname = os.path.splitext(strucfile)[0] + "_%s_new.fps" % fpname
            if os.path.isfile(fp_fname):
                #print "%s exists!" % fp_fname
                return chemfp.load_fingerprints(fp_fname)
            else:
                fparena = chemfp.load_fingerprints(chemfp.read_structure_fingerprints(
                                                                         fptype,
                                                                         source=strucfile
                                                                                             ))
                fparena.save(fp_fname)
                print "%s saved!" % fp_fname
                return fparena
        activefps = load_fps(activesfn)
        decoyfps = load_fps(decoysfn)
        crystalpath = os.path.join(target_dir, "PDB", "ligands%s.sdf" % INFIX)
        if DB == "common":
            crystalpath = crystalpath.replace("common" + os.sep, "")
        crystalfps = load_fps(crystalpath)

        a_search = chemfp.knearest_tanimoto_search(activefps, crystalfps, k=1, threshold=0)
        d_search = chemfp.knearest_tanimoto_search(decoyfps, crystalfps, k=1, threshold=0)
        print target
        #print  [(query_id, hits.get_ids_and_scores(), hits.get_ids_and_scores()) for query_id, hits in a_search]# if hits.get_ids_and_scores()]
        results = [(query_id, hits.get_ids_and_scores()[0][1], hits.get_ids_and_scores()[0][0]) for query_id, hits in a_search if hits.get_ids_and_scores()]
        results += [(query_id, hits.get_ids_and_scores()[0][1], hits.get_ids_and_scores()[0][0]) for query_id, hits in d_search if hits.get_ids_and_scores()]
        results.sort(key=lambda x: x[1], reverse=True)
        if not results:
            print "No results for %s and %s" % (target, fptype)
        if DECOYS:
            outfn = os.path.join(decoystdir, target, fpname + '_new.csv')
        else:
            outfn = os.path.join(target_dir, fpname + '_new.csv')
        output = open(outfn, 'w')
        csv.writer(output).writerows(results)
        #print "Saving results to %s" % outfn
        output.close()
    return "%s done" % target

if DB != "common":
    targets = [target for target in alltargets if  os.path.isfile(os.path.join(root,target,"PDB/ligands%s.sdf" % INFIX))]
else:
    targets =  alltargets
targets = [t for t in targets if t not in badtargets] #Target with no decoys? DUDE bug?
#targets = ['KITH']
print targets

n = len(targets)
m = 0

pool = multiprocessing.Pool()
for r in pool.imap_unordered(fingerprint_screen, targets):
    m += 1
    print "(%s/%s) %s" % (m, n, r)
#for target in targets:
#    r = fingerprint_screen(target)
#    m += 1
#    print "(%s/%s) %s" % (m, n, r)
