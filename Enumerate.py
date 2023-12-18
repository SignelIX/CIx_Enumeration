from rdkit import Chem
from rdkit.Chem import rdChemReactions
import json
import toml
import pandas as pd
from numpy import random
from multiprocessing.pool import ThreadPool as Pool
import threading
import MolDisplay
import ChemUtilities
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from rdkit.Chem import SaltRemover
import pathlib
from tqdm import tqdm
import re
import gc
import os, psutil
import csv
import time
import copy
import argparse
import yaml
import numexpr
import math

CPU_COUNT = os.cpu_count()
NUM_WORKERS = CPU_COUNT * 2
chksz = 50000
numexpr.set_num_threads(numexpr.detect_number_of_cores())
if "NUMEXPR_MAX_THREADS" not in os.environ:
    os.environ["NUMEXPR_MAX_THREADS"] = '16'
else:
    NUM_WORKERS = int(os.environ["NUMEXPR_MAX_THREADS"] )


class Enumerate:
    smilescol = '__SMILES'
    idcol = '__BB_ID'
    rxncol = '__RXNNAME'
    rxndict = {}
    named_reactions = None
    bb_info_dfs = None
    bb_rxn_dict = None
    LookupReactants = False
    smiles_colnames = ['SMILES', 'smiles']
    bbid_colnames = ['BB_ID', 'id']
    rxn_colnames = ['reaction_name', 'reaction']
    reagent_prefix = 'r'
    rxnscheme = None
    scheme_reactants = None
    scheme_name = None
    scheme_source = None

    #region Helpers
    @staticmethod
    def Deduplicate( outpath, outsuff):
        df = pd.read_csv(outpath + '.' + outsuff + '.all.csv')
        df = df.drop_duplicates(keep='first', subset=['full_smiles'])
        df.to_csv(outpath + '.' + outsuff + '.all.dedup.csv')
        return outpath + '.' + outsuff + '.all.dedup.csv'

    #endregion Helpers

    #region Reactions
    @staticmethod
    def React_Molecules (r1, r2, SM_rxn, showmols):

        if showmols == True:
            print ('Reaction Params:')
            print(r1, r2, SM_rxn)
        if isinstance (SM_rxn,str):
            SM_rxn = [SM_rxn]

        if r2 == 'None' or r2 is None or r2 == 'NOSMILES' or r2 == "SKIPCYC":
            if not isinstance(r1, Chem.rdchem.Mol):
                m1 = Chem.MolFromSmiles(r1)
            else:
                m1 = r1
            reacts = [m1]
        else:
            if not isinstance(r1, Chem.rdchem.Mol):
                m1 = Chem.MolFromSmiles(r1)
            else:
                m1 = r1
            if not isinstance(r2, Chem.rdchem.Mol):
                m2 = Chem.MolFromSmiles(r2)
            else:
                m2 = r2
            reacts = (m1, m2)

        prod_ct = 0
        for r in SM_rxn:
            if r == "product":
                smi = r1
            else:
                rxn = rdChemReactions.ReactionFromSmarts(r)
                products = rxn.RunReactants(reacts)
                if len(products) == 0:
                    product = None
                else:
                    if len(products) > 1:

                        proddict= {}

                        for cx in range (0, len (products)):
                            if Chem.MolToSmiles(products[cx][0]) not in proddict:
                                proddict [Chem.MolToSmiles(products[cx][0])] = 1
                                prod_ct += 1
                            else:
                                proddict[Chem.MolToSmiles(products[cx][0])] += 1

                    product = products [0][0]

                try:
                    if product is  None:
                        smi = None
                    else:
                        smi = Chem.MolToSmiles(product)
                        break
                except:

                    smi = None
        if showmols == True:
            if smi is not None:
                MolDisplay.ShowMols([smi, r1, r2])
            else:
                MolDisplay.ShowMols([ r1, r2])
        return smi, prod_ct, products

    def ReadRxnScheme (self, fname_json, schemename, verbose = True, FullInfo = False, Store = False):

        try:
            data = {}
            data[schemename] = json.loads(fname_json)
        except Exception as e:
            try:
                f = open (fname_json)
                data = json.load (f)
                f.close ()
            except Exception as E:
                if verbose:
                    print ('json read error:', str(E))
                if Store:
                    self.scheme_reactants = None
                    self.rxnscheme = None
                return None, None
        if schemename not in data:
            if verbose:
                print('schemename error:', schemename)
            if Store:
                self.scheme_reactants = None
                self.rxnscheme = None
            return None, None
        if 'reactants' not in data[schemename]:
            reactants = None
        else:
            reactants = data[schemename]["reactants"]

        if Store:
            self.scheme_reactants = reactants
            self.rxnscheme = data[schemename]
            self.schemename = schemename
            self.scheme_source = fname_json
        if FullInfo:
            return data[schemename], reactants
        else:
            return data[schemename]["steps"], reactants
    def Load_NamedReactions (self, rxnfiles):
        self.named_reactions = {}
        namedrxnsjson, etc = self.ReadRxnScheme(self.scheme_source, 'Named_Reactions', FullInfo=True, Store=False)
        self.named_reactions.update(namedrxnsjson)
        if rxnfiles is not None:
            for rxl in rxnfiles:
                if '/' not in rxl:
                    rxl = str(pathlib.Path(__file__).parent.parent) + '/' + rxl
                with open(rxl, "r") as jsonFile:
                    data = json.load(jsonFile)
                self.named_reactions.update(data)
    def Add_Library_BBLists (self, librarypath, library ):
        bbfiles = self.Get_BBFiles(librarypath + '/' + library)
        dfs = self.load_BBlists(bbfiles.values())
        dflist = sorted(dfs.items())
        dflist = [x[1] for x in dflist]
        self.Set_BB_Info(dflist)
    def getReactants (self, rx, reactants, prior_product, in_reactants, step, rxtants):
        if step["Reactants"][rx] == 'None':
            reactants.append ('None')
        elif step["Reactants"][rx] == 'NOSMILES':
            reactants.append ('None')
        elif step["Reactants"][rx] == 'SKIPCYC':
            reactants.append('None')
        elif step["Reactants"][rx] == 'p':
            reactants.append(prior_product)
        elif rxtants != None and step["Reactants"][rx] in rxtants:
            idx = rxtants.index (step["Reactants"][rx])
            reactants.append(in_reactants[idx])
        elif step["Reactants"][rx].startswith (self.reagent_prefix):
            idx = int(step["Reactants"][rx][len (self.reagent_prefix):])
            reactants.append(in_reactants[idx])
        else:
            otherrxtidx = abs(rx - 1)
            if step['Reactants'][otherrxtidx].startswith (self.reagent_prefix):
                ridx = int(step['Reactants'][otherrxtidx][len(self.reagent_prefix):])
                mappedcol =step["Reactants"][rx]
                if mappedcol in self.bb_info_dfs[ridx].columns:
                    bbdf = self.bb_info_dfs[ridx]
                    otherrxtntsmiles = in_reactants[ridx]
                    if otherrxtntsmiles != 'C':
                        colvalsmiles = bbdf[bbdf[self.smilescol] == otherrxtntsmiles ].iloc [0][ mappedcol]
                    else:
                        colvalsmiles = bbdf.iloc[0][mappedcol]
                    reactants.append(colvalsmiles)
                    return reactants
            reactants.append(step["Reactants"][rx])
        return reactants

    def Set_BB_Info (self, bbdfs):
        self.bb_info_dfs = bbdfs
        for bx in range(len(self.bb_info_dfs)):
            self.bb_info_dfs [bx] = self.Reset_BBDF_ColNames (self.bb_info_dfs[bx])
        self.bb_rxn_dict = {}
        for bx in range(len(self.bb_info_dfs)):
            if self.rxncol in self.bb_info_dfs [bx].columns:
                self.bb_rxn_dict [bx] = {}
                for rx, row in self.bb_info_dfs [bx].iterrows ():
                    smiles = row [self.smilescol]
                    if not type(smiles) == str:
                        smiles = 'None'
                    if type (row [self.rxncol]) == str:
                        self.bb_rxn_dict[bx][smiles] = row [self.rxncol]

    def Lookup_Rxns_From_Reactants(self, rxtnt_names, smiles):
        rxdict = {}
        preflen = len(self.reagent_prefix)
        if self.bb_info_dfs is not None and self.LookupReactants == True:
            for rxt in range(len(rxtnt_names)):
                name =rxtnt_names [rxt]
                smilesval = smiles [rxt]
                if name.startswith(self.reagent_prefix):
                    rxidx = int(name[preflen:])
                    if rxidx in self.bb_rxn_dict:
                        if smilesval in self.bb_rxn_dict[rxidx]:
                            rxdict[rxidx] = self.bb_rxn_dict [rxidx][smilesval]
        return rxdict

    def Pull_ReactionDetails (self, prior_product, in_reactants, step, rxtants):
        reactants = []
        reactants = self.getReactants(0, reactants, prior_product, in_reactants, step, rxtants)
        reactants = self.getReactants(1, reactants, prior_product, in_reactants, step, rxtants)
        reactions_dict = self.Lookup_Rxns_From_Reactants (step ['Reactants'], reactants)


        if reactions_dict != {} :
            reaction = reactions_dict [list(reactions_dict.keys ())[0]]
        elif "default" in step["Rxns"] and len(step["Rxns"]) == 1:
            reaction = step["Rxns"]["default"]
        else:
            reaction = 'unassigned'
            for rxn in step["Rxns"]:
                if rxn != 'default':
                    if rxn == 'product':
                        reaction = 'product'
                        break
                    if rxn == 'NOSMILES':
                        if reactants [0] == 'NOSMILES' or reactants [1] == 'NOSMILES':
                            reaction = step["Rxns"][rxn]
                            break
                        else:
                            reaction = 'FAIL'
                            continue
                    if rxn == 'C' or rxn == '[C;H4]':
                        if reactants[0] == 'C' or reactants [1] == 'C':
                            reaction = step["Rxns"][rxn]
                            break
                    if rxn == 'SKIPCYC':
                        if reactants[0] == 'SKIPCYC' or reactants[1] == 'SKIPCYC':
                            reaction = step["Rxns"][rxn]
                            break
                        else:
                            reaction = 'FAIL'
                            continue
                    if reactants [1] == 'None':
                        if not isinstance(reactants [0], Chem.rdchem.Mol):
                            m = Chem.MolFromSmiles(reactants [0])
                        else:
                            m = reactants [0]
                    else:
                        if not isinstance(reactants[1], Chem.rdchem.Mol):
                            m = Chem.MolFromSmiles(reactants[1])
                        else:
                            m = reactants[1]
                    pattern = Chem.MolFromSmarts(rxn)
                    if m is None:
                        reaction = 'FAIL'
                    else:
                        if m.HasSubstructMatch(pattern):
                            reaction =  step["Rxns"][rxn]
                            break
                        else:
                            reaction ='FAIL'
            if reaction == 'FAIL' and "default" in step["Rxns"] :
                reaction = step["Rxns"]["default"]
        return reactants, reaction

    def RunRxnScheme (self, in_reactants, schemefile_jsontext, schemename, showmols, schemeinfo = None):
        if schemeinfo is not None:
            scheme, rxtants = [schemeinfo[0], schemeinfo [1]]
            if 'steps' in  scheme:
                scheme = scheme ['steps']
        else:
            scheme, rxtants = self.ReadRxnScheme(schemefile_jsontext, schemename)
        if scheme is None:
            return 'NOSCHEMEFAIL',0, None
        p=in_reactants [0]
        prod_ct = 1
        intermeds = []

        for step in scheme:
            stepname = step
            if type (scheme) == dict:
                step = scheme [step]
            reactants, reaction = self.Pull_ReactionDetails (p, in_reactants, step, rxtants)
            if reaction == 'FAIL':
                p = 'FAIL'
                break
            if reaction == 'terminate':
                break
            if reaction == 'product':
                continue
            else:

                if self.named_reactions is not None:
                    rlist = []
                    if type (reaction) == str:
                        reaction = [reaction]
                    for r in reaction:
                        if r in self.named_reactions:
                            if type (self.named_reactions [r] ) == str:
                                rlist.append (self.named_reactions [r])
                            else:
                                rlist.extend (self.named_reactions [r])
                        else:
                            rlist.append (r)
                    reaction = rlist

                try:
                    p, outct, products = self.React_Molecules(reactants [0], reactants [1],  reaction,  showmols)
                except Exception as e:
                    p = None
                    outct = 0
                if outct > 0:
                    prod_ct *= outct
            if showmols:
                print(stepname, p, reactants)
            if p is None:
                p = 'FAIL'
                break
            intermeds.append (p)

        return p, prod_ct, [scheme, rxtants, intermeds]
    #endregion Reactions

    def Reset_BBDF_ColNames (self, bbdf):
        changecoldict = {}
        if self.smiles_colnames is not None and len(self.smiles_colnames) > 0:
            for smc in self.smiles_colnames:
                if smc in bbdf.columns:
                    changecoldict[smc] = self.smilescol
                    break
        if self.bbid_colnames  is not None and len(self.bbid_colnames ) > 0:
            for bbidc in self.bbid_colnames :
                if bbidc in bbdf.columns:
                    changecoldict[bbidc] = self.idcol
                    break
        if self.rxn_colnames is not None and len(self.rxn_colnames) > 0:
            for rxnidc in self.rxn_colnames:
                if rxnidc in bbdf.columns:
                    changecoldict[rxnidc] = self.rxncol
                    break
        if len(changecoldict) > 0:
            bbdf = bbdf.rename(columns=changecoldict)
        return bbdf

    #region BBs
    def load_BBlists(self, bblists):
        cycs = []
        bbdict = {}

        if type (bblists) is dict:
            bbdict = bblists
        else:
            for l in bblists:
                cyclesplit = l.split('.')
                cyc = cyclesplit[len(cyclesplit) - 2]
                cycs.append(cyc)
                if type(l) is str:
                    if '.' not in l:
                        bbdict [cyc] = None
                    else:
                        bbdict[cyc] = pd.read_csv(l)
        for cyc in bbdict.keys():
            if bbdict [cyc] is not None:
                bbdict[cyc] = self.Reset_BBDF_ColNames (bbdict [cyc])
                bbdict[cyc] = bbdict[cyc].drop_duplicates(subset=self.idcol, keep="first").reset_index(drop=True)
        return bbdict

    def pull_BBs(self, inpath, idcol, smilescol):
        globlist = pathlib.Path(inpath).glob('*.csv')
        full_df = None
        for f in tqdm(globlist):
            if full_df is None:
                full_df = pd.read_csv(f)
                full_df = full_df[[idcol, smilescol]]
                full_df = full_df.rename(columns={idcol: self.idcol, smilescol: self.smilescol})
            else:
                df = pd.read_csv(f)
                df = df.rename(columns={idcol: self.idcol, smilescol: self.smilescol})
                try:
                    df = df[[self.idcol, self.smilescol]].dropna()
                except:
                    print('exiting:', f)
                    exit()
                full_df = full_df.append(df, ignore_index=True)
        return full_df
    @staticmethod
    def Get_BBFiles ( bb_specialcode, lib_subfolder, enumsfolder, libname):
        libspecprefix = ''
        if lib_subfolder != '' and lib_subfolder is not None:
            libspecprefix = '/' + lib_subfolder
        bbpath = enumsfolder + libname + libspecprefix + '/BBLists'
        if bb_specialcode is not None and bb_specialcode != '':
            srchstr = libname + '.' + bb_specialcode +  '.BB?.csv'
            flist = pathlib.Path(bbpath).glob(srchstr)
        else:
            srchstr = libname + '.BB?.csv'
            flist = pathlib.Path(bbpath).glob(srchstr)
        infilelist = []

        for f in flist:
            c = str(f)
            result = re.search(r'\.BB[0-9]+\.', c)
            if result is not None:
                infilelist.append(c)
        infilelist.sort()
        if len(infilelist) == 0:
            print('FAIL: No BB files found with format ' + srchstr + ' found in ' + bbpath, infilelist)
            print ('Inputs:BBcode:' , bb_specialcode, ' lib_subfldr:', lib_subfolder, ' enum folder:', enumsfolder, ' lib:', libname)
            return ''

        return infilelist

    @staticmethod
    def Get_BBFiles( inpath):
        flist = pathlib.Path(inpath).glob('*.csv')
        bbfiles = {}
        for f in flist:
            c = str(f)
            result = re.search(r'\.BB[0-9]+\.', c)
            if result is not None:
                bbfiles[result.group(0)[3:4]] = c
        return bbfiles
    # endregion BBs

    #region Enumerations
    def EnumFromBBFiles(self, libname, bb_specialcode, out_specialcode, enumsfolder, lib_subfolder,
                        num_strux, rxschemefile, picklistfile=None, SMILEScolnames = [], BBcolnames = [],
                        rem_dups = False, returndf = False, write_fails_enums = True, retIntermeds = False, overrideBB ={}):

        #find and load BB inputs
        infilelist = self.Get_BBFiles (bb_specialcode, lib_subfolder, enumsfolder, libname)
        if type (infilelist) == str:
            return

        #load default reaction scheme file if not specified
        if rxschemefile is None:
            rxschemefile = enumsfolder + 'RxnSchemes.json'

        #set up output file
        samplespath = enumsfolder  + libname + '/' + lib_subfolder + '/Samples/'
        if not os.path.exists(samplespath):
            os.makedirs(samplespath)
        outpath = samplespath + libname
        if  out_specialcode != '' and out_specialcode is not None:
            out_specialcode += '.' + out_specialcode

        if returndf is True:
            outpath = None

        #Enumerate
        outfile = self.enumerate_library_strux(libname, rxschemefile, infilelist, outpath, num_strux, picklistfile,
                                               SMILEScolnames=SMILEScolnames, BBIDcolnames=BBcolnames,
                                               removeduplicateproducts=rem_dups, write_fails_enums=write_fails_enums,
                                               retIntermeds = retIntermeds, overrideBB = overrideBB)
        return outfile

    def TestReactionScheme(self,schemename, rxtnts, rxnschemefile, retIntermeds = False):
        for r in range (len(rxtnts)):
            if  (type (rxtnts[r]) is not str) :
                rxtnts [r] = 'None'
        res = self.RunRxnScheme(rxtnts, rxnschemefile, schemename, False)
        rxnslist = []
        for k in res [2][0]:
            rxnslist.append([str(k['Reactants']) + '<p>' + str(k['Rxns'])])
        if res[1] > 1:
            if not retIntermeds:
                return 'FAIL--MULTIPLE PRODUCTS'
            else:
                return 'FAIL--MULTIPLE PRODUCTS', None, None
        if not retIntermeds:
            return res[0]
        else:
            if res [2] is not None:
                return res[0], res[2][2], rxnslist
            else:
                return res[0], None, None

    def rec_bbpull(self, bdfs, level, cycct, bbslist, ct, reslist, fullct, hdrs,
                   libname, rxschemefile, outpath, rndct, removeduplicateproducts,
                   currct=0, appendmode=False, retIntermeds=False):
        if reslist is None:
            reslist = [[]] * min(chksz, fullct)

        for i, b in bdfs[level].iterrows():
            bb = b[self.idcol]
            bbs = b[self.smilescol]
            rxn = b[self.rxncol]
            if level == 0:
                bbslist = []

            bblevellist = copy.deepcopy(bbslist)
            bblevellist.append(bb)
            bblevellist.append(bbs)
            bblevellist.append(rxn)

            if level < cycct - 1:
                ct, reslist, currct, appendmode = self.rec_bbpull(bdfs, level + 1, cycct, bblevellist, ct, reslist,
                                                            fullct, hdrs,
                                                            libname, rxschemefile, outpath, rndct,
                                                            removeduplicateproducts,
                                                            currct=currct, appendmode=appendmode,
                                                            retIntermeds=retIntermeds)
            else:
                reslist[currct] = bblevellist

                ct += 1
                currct += 1

                if currct == chksz or ct == fullct:
                    enum_in = pd.DataFrame(reslist, columns=hdrs)
                    # Enumerate
                    self.DoParallel_Enumeration(enum_in, hdrs, libname, rxschemefile, outpath, cycct, rndct,
                                                removeduplicateproducts, appendmode=appendmode,
                                                retIntermeds=retIntermeds)
                    reslist = [[]] * min(chksz, fullct - ct)
                    currct = 0
                    appendmode = True
                    gc.collect()

        return ct, reslist, currct, appendmode
    def enumerate_library_strux(self, libname, rxschemefile, infilelist, outpath, rndct=-1, bblistfile=None,
                                SMILEScolnames = [], BBIDcolnames = [], removeduplicateproducts = False,
                                outtype = 'filepath', write_fails_enums = True, retIntermeds=False, overrideBB = {}):


        #set up default values
        if SMILEScolnames is None:
            SMILEScolnames = []
        if BBIDcolnames is None:
            BBIDcolnames = []

        infilelist.sort ()
        cycct = len (infilelist)

        if type(infilelist) is list:
            cycdict = self.load_BBlists(infilelist)
            bdfs = [None] * cycct
            if bblistfile is not None:
                picklistdf = pd.read_csv(bblistfile)
                for ix in range (0, cycct):
                    bdfs [ix] = picklistdf[picklistdf['Cycle'] == 'BB' + str (ix + 1)]
                    bdfs = bdfs.merge(cycdict['BB' + str (ix +1)], on=[self.idcol], how='left')
            else:
                for ix in range(0, cycct):
                    bdfs[ix] = cycdict['BB' + str(ix+1)]

            fullct = 1
            for ix in range(0, cycct):
                fullct *= len(bdfs[ix])

            if rndct > fullct:
                rndct = -1

            hdrs = []
            for ix in range(0, cycct):
                hdrs.append('bb' + str(ix + 1))
                hdrs.append('bb' + str(ix + 1) + '_smiles')

            appendmode = False

            if rndct == -1 :
                #generate the full combinatorial library
                hdrs = []
                for ix in range(0, cycct):
                    hdrs.append('bb' + str(ix + 1))
                    hdrs.append('bb' + str(ix + 1) + '_smiles')
                with open(outpath + ".EnumList.csv", "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(hdrs)
                self.rec_bbpull (bdfs, 0, cycct, [], 0, None, fullct, hdrs,
                            libname, rxschemefile, outpath, rndct, removeduplicateproducts,
                            appendmode = appendmode,
                            retIntermeds=retIntermeds)
            else:
                reslist = [[]] * min(rndct, chksz)
                ct = 0
                currct = 0

                #generate a random set of bbs to enumerate
                while ct < rndct:
                    bblist = []
                    for ix in range (0, cycct):
                        if str (ix+1) in overrideBB :
                            bb = overrideBB [str(ix+1)]
                            bbs = overrideBB [str(ix+1)]
                            bblist.append(bb)
                            bblist.append(bbs)
                        else:
                            ri = random.randint (0, len (cycdict['BB' + str (ix + 1)]))
                            b = cycdict['BB' + str (ix + 1)].iloc[ri]
                            bb = b[self.idcol]
                            bbs = b[self.smilescol]
                            bblist.append (bb)
                            bblist.append (bbs)

                    #Enumerate if not already found
                    if bblist not in reslist:
                        reslist[currct] = bblist
                        ct += 1
                        currct += 1
                        if currct ==  chksz  or ct == rndct:
                            enum_in = pd.DataFrame(reslist, columns=hdrs)
                            #Enumerate
                            outpath = self.DoParallel_Enumeration(enum_in, hdrs, libname, rxschemefile, outpath, cycct,
                                    rndct, removeduplicateproducts, appendmode=appendmode, retIntermeds=retIntermeds)
                            reslist = [[]] * min (chksz, rndct - ct)
                            currct = 0
                            appendmode = True
        else:
            #enumerate based on input dataframe
            resdf = infilelist
            enum_in = resdf
            hdrs = None

            # Enumerate
            outpath = self.DoParallel_Enumeration(enum_in, hdrs, libname, rxschemefile, outpath, cycct, rndct,
                       removeduplicateproducts, appendmode = False, write_fails_enums=write_fails_enums,
                       retIntermeds=retIntermeds)

        return outpath

    def DoParallel_Enumeration (self, enum_in, hdrs, libname, rxschemefile, outpath, cycct, rndct=-1,
                                removeduplicateproducts = False, appendmode = False, write_fails_enums=True,
                                retIntermeds = False):
        def taskfcn(row, libname, rxschemefile, showstrux, schemeinfo, cycct, retIntermeds = False):
            rxtnts = []
            for ix in range(0, cycct):
                rxtnts.append(row['bb' + str(ix + 1) + '_smiles'])
            try:
                res, prodct, schemeinfo = self.RunRxnScheme(rxtnts, rxschemefile, libname, showmols=showstrux,
                                                            schemeinfo=schemeinfo)
                if prodct > 1:
                    print ('FAIL--MULTIPLE PRODUCTS', res)
                    res =  ['FAIL--MULTIPLE PRODUCTS']
                else:
                    res = [res]
            except Exception as e:
                print(str(e))
                res =  ['FAIL']

            if retIntermeds:
                for k in range (0, len(schemeinfo[0])):
                    res.append(schemeinfo[0][k])
                    if  k < len(schemeinfo[2]):
                        res.append (schemeinfo[2][k])
                    else:
                        res.append (None)
                return res
            return res

        def processchunk (resdf, df, outpath, retIntermeds, schemeinfo):
            pbar = ProgressBar()
            pbar.register()
            ddf = dd.from_pandas(resdf, npartitions=CPU_COUNT * 10)
            if len (ddf) <  5000:
                wrkrs = 1
            else:
                wrkrs  = NUM_WORKERS
            metaser = [(0, object)]
            colnames = {0: 'full_smiles'}
            if retIntermeds:
                rxnct = len(schemeinfo [0])
                for r in range(0,rxnct):
                    metaser.append ((2*r + 1, str))
                    metaser.append((2 * r + 2, str))
                    colnames [2*r + 1] = 'step' + str(r) + '_rxn'
                    colnames [2*r + 2] = 'step' + str(r) + '_intermed'
            res = ddf.apply(taskfcn, axis=1, result_type='expand',
                            args=(libname, rxschemefile, rndct == 1, schemeinfo, cycct, retIntermeds),
                            meta=metaser).compute(scheduler='processes',  num_workers=wrkrs)
            pbar.unregister()
            moddf = resdf.merge(res, left_index=True, right_index=True)
            moddf = moddf.rename(columns=colnames)

            if outpath is not None:
                enumdf = moddf[~moddf['full_smiles'].isin( ['FAIL','FAIL--MULTIPLE PRODUCTS'])]
                enumdf.to_csv(outpath + '.' + outsuff + '.enum.csv', mode='a', index=False, header=False)
                faildf = moddf[moddf['full_smiles'].isin(['FAIL','FAIL--MULTIPLE PRODUCTS'])]
                if write_fails_enums:
                    faildf.to_csv(outpath + '.' + outsuff + '.fail.csv', mode='a', index=False, header=False)
                    moddf.to_csv(outpath + '.' + outsuff + '.all.csv', mode='a', index=False, header=False)
                gc.collect()
                return None
            else:
                if df is None:
                    df = moddf
                else:
                    df.append(moddf, ignore_index=True)
                gc.collect()
                return df

        #read in reaction scheme
        outsuff = 'full'
        if rndct != -1:
            outsuff = str(rndct)
        hdrstr = ','.join (hdrs)
        schemeinfo = self.ReadRxnScheme(rxschemefile, libname, False)
        if (schemeinfo is None):
            print ('Failed, scheme not detected')
            raise(13)

        #set up for writing to file/files
        if outpath is not None:
            if not write_fails_enums:
                flist = [outpath + '.' + outsuff + '.all.csv',None, None]
            else:
                flist = [outpath + '.' + outsuff + '.all.csv', outpath + '.' + outsuff + '.fail.csv', outpath + '.' + outsuff + '.enum.csv']
            if appendmode == False:
                for fname in flist:
                    f = open(fname, 'w')
                    f.write(hdrstr +',full_smiles')
                    rxnct = len(schemeinfo[0])
                    if retIntermeds:
                        for r in range(0, rxnct):
                            f.write( ',step' + str(r) + '_rxn')
                            f.write(',step' + str(r) + '_intermed')
                    f.write ('\n')
                    f.close()

        #process file in chunks or process enum_in as a dataframe
        if type (enum_in) is str:
            reader = pd.read_csv(enum_in, chunksize=chksz)
            df = None
            cct = 0
            for chunk in reader:
                resdf = pd.DataFrame (chunk, columns = hdrs)
                df = processchunk(resdf, df, outpath , schemeinfo=schemeinfo, retIntermeds=retIntermeds)
                cct += 1
        else:
            df = None
            df = processchunk(enum_in, df, outpath, schemeinfo=schemeinfo,retIntermeds=retIntermeds)

        #save and return filename or return dataframe of results
        if outpath is not None:
            path = outpath + '.' + outsuff + '.all.csv'
            if removeduplicateproducts:
                path = self.Deduplicate (outpath, outsuff)
            return path
        else:
            if removeduplicateproducts:
                df = df.drop_duplicates(keep='first', subset = ['full_smiles'])
            return df
    #endregion Enumerations

class EnumerationCLI :
    @staticmethod
    def Run_CLI (SMILEScolnames = None, BBcolnames = None):
        paramdefaults = [ ('rxnschemefile', './RxnSchemes.json'), ('schemepath','.'), ('scheme',''), ('schemespec',''), ('numstrux', 5000), ('removedups', False)]
        parser = argparse.ArgumentParser(description='Enumeration Options')
        parser.add_argument('-p', '--paramfile', nargs='?', default=None, type=str,
                            help='optional .yaml file for commandline paramaters')
        parser.add_argument('-r', '--rxnschemefile', nargs='?', default=None, type=str,
                            help='Rxnschemes.json file path')
        parser.add_argument('-sp', '--schemepath', nargs='?', default=None, type=str, help='Enumerations folder path')
        parser.add_argument('-s', '--scheme', nargs='?', default=None, type=str, help='Scheme Name')
        parser.add_argument('-sx', '--schemespec', nargs='?', default=None, type=str, help='sub-scheme Specifier')
        parser.add_argument('-bx', '--bbspecialcode', nargs='?', default=None, type=str, help='BB files special code')
        parser.add_argument('-n', '--numstrux', nargs='?', default=None, type=int,
                            help='number of structures to enumerate (-1 for all)')
        parser.add_argument('-rd', '--removedups', nargs='?', default=None, type=str,
                            help='Remove duplicate structures (True/False)')
        parser.add_argument('-wfe', '--write_fails_enums', nargs='?', default=None, type=str,
                            help='Write Fails and Enumerated Molecules in separate files (True/False)')
        parser.add_argument('-rxntest', '--rxntest', nargs=2, default=None, type=str,
                            help='2 reactants SMILES')
        parser.add_argument('-rxn', '--rxn', nargs=1, default=None, type=str,
                            help='Reaction SMARTS')
        args = vars(parser.parse_args())

        if args['paramfile'] is not None:
            with open(args['paramfile'], 'r') as stream:
                try:
                    prms = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print(exc)
            for k in prms.keys ():
                if k not in args or args [k] is None or args[k] == '':
                    args [k] = prms[k]
            for x in paramdefaults:
                if x[0] not in args or args [x[0]] is None or args[x[0]] == '':
                    args [x[0]] = x[1]
        specstr = args['schemespec']
        addspec = ''
        if specstr != '' and specstr is not None:
            addspec = '/' + specstr

        if 'rxntest' in args:
            enum = Enumerate()
            res = enum.React_Molecules(args ['rxntest'] [0],args ['rxntest'] [1] , args['rxn'], True)
            exit ()

        if 'removedups' in args:
            if args['removedups'] not in [True, False]:
                if args['removedups'] == 'True':
                    rd = True
                elif args['removedups'] == 'False':
                    rd = False
                else:
                    'Fail: rd must be True or False'
                    exit()
            else:
                rd = args["removedups"]

        if 'write_fails_enums' in args:
            if args['write_fails_enums'] not in [True, False]:
                if args['write_fails_enums'] == 'True':
                    write_fails_enums = True
                elif args['write_fails_enums'] == 'False':
                    write_fails_enums = False
                else:
                    print ('Fail: write_fails_enums must be True or False')
                    exit()
            else:
                write_fails_enums = args["write_fails_enums"]

        enum = Enumerate()

        outf = enum.EnumFromBBFiles(args['scheme'], args['bbspecialcode'], args['schemespec'], args['schemepath'],
                             addspec, args['numstrux'],
                             args['rxnschemefile'], SMILEScolnames=SMILEScolnames, BBcolnames=BBcolnames, rem_dups=rd, write_fails_enums=write_fails_enums)


if __name__=="__main__":
    EnumerationCLI.Run_CLI()




