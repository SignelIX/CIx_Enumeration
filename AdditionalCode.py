#class HelperFunctions:
    # def Sample_Library (self, BBlistpath, outfilepath, schemepath, scheme_name,  rndmsample_ct, ShowMols, saltstrip = True, CountOnly=False ):
    #     def recurse_library (prev_rvals, prev_idvals, cyc_num, cycct, enum_list, outfile, ct):
    #         for idx, row in cyc[cycles[cyc_num]].iterrows ():
    #             rvals = prev_rvals.copy ()
    #             rvals.append (row [self.smilescol])
    #             idvals = prev_idvals.copy ()
    #             idvals.append(str(row ['ID']))
    #
    #             if cyc_num == cycct -1:
    #                 try:
    #                     enum_res, prod_ct = self.RunRxnScheme(rvals, schemepath, scheme_name, ShowMols, in_reactions=None)
    #                 except:
    #                     enum_res = 'FAIL'
    #                 outfile.write(enum_res)
    #                 for ir in range(0, len(rvals)):
    #                     outfile.write(',' + str(rvals[ir]) + ',' + str(idvals[ir]))
    #                 outfile.write('\n')
    #                 ct += 1
    #
    #                 if ct %1000 ==0:
    #                     print (ct)
    #             else:
    #                 ct = recurse_library (rvals, idvals, cyc_num + 1, cycct, enum_list, outfile, ct)
    #         return ct
    #
    #     def EnumFullMolecule (cycles, cycvals, r_list):
    #         rvals = []
    #         idvals = []
    #
    #         for c_i in range(0, len(cycles)):
    #             rnum = random.randint (0,len(cycvals [c_i]))
    #             rcx = cycvals [c_i].iloc [rnum]
    #             rvals.append (rcx[self.smilescol])
    #             idvals.append (rcx ['ID'])
    #         string_ids = [str(int) for int in idvals]
    #         idstr = ','.join(string_ids)
    #         if not  idstr  in r_list:
    #              r_list[idstr] = True
    #              enum_res, prod_ct = self.RunRxnScheme(rvals, schemepath, scheme_name, ShowMols, in_reactions=None)
    #              resline = enum_res
    #              for ir in range(0, len(rvals)):
    #                   resline +=',' + rvals[ir] + ',' + str(idvals[ir])
    #         else:
    #             return None
    #         return resline
    #
    #     class CompleteEnum ():
    #         tic = None
    #         block = ''
    #         lock = None
    #         outfile = None
    #         def __init__ (self, outf):
    #             self.tic = time.perf_counter()
    #             self.lock = threading.Lock()
    #             self.outfile = outf
    #         enumct = 0
    #
    #         def CompleteEnumAsync (self, res):
    #             if res is not None:
    #                 self.lock.acquire()
    #                 self.enumct = self.enumct + 1
    #                 self.block += res
    #                 self.block += '\n'
    #                 if self.enumct % 1000 == 0:
    #                     print(self.enumct)
    #                     toc = time.perf_counter()
    #                     print(f"Time Count {toc - self.tic:0.4f} seconds")
    #                     self.WriteBlock ()
    #                 self.lock.release()
    #
    #         def WriteBlock (self):
    #             self.outfile.write(self.block)
    #             self.block = ''
    #
    #     df = pd.read_csv (BBlistpath)
    #     df = df.fillna(0)
    #
    #     scheme, reactants= self.ReadRxnScheme(schemepath, scheme_name)
    #     if reactants is None:
    #         cycles = df.Cycle.unique ()
    #         cycles.sort ()
    #     else:
    #         cycles = reactants
    #
    #     cyc = {}
    #     if saltstrip == True:
    #         df = ChemUtilities.SaltStripMolecules(df)
    #
    #     for c in cycles:
    #         cyc[c] = df[df.Cycle.isin([c])]
    #
    #     if CountOnly:
    #         sz = 1
    #         for cx in cyc:
    #             sz *= len (cyc[cx])
    #             print (cx, str(len (cyc[cx])))
    #         print ('libsize = ', sz)
    #         return ()
    #
    #     ct = 0
    #     outfile = open(outfilepath, 'w')
    #     enum_list = []
    #
    #     outfile = open(outfilepath, 'w')
    #
    #     outfile.write('product_SMILES')
    #
    #     for c_i in range(0, len(cycles)):
    #         outfile.write(',' + str(cycles [c_i]) + ',ID_' +  str(cycles [c_i]) )
    #
    #     outfile.write('\n')
    #
    #     if rndmsample_ct == -1:
    #         rvals = []
    #         idvals = []
    #         recurse_library (rvals, idvals, 0, len(cycles), enum_list, outfile, 0)
    #     else:
    #         cycvals = []
    #         for c_i in range (0,len(cycles)):
    #             cycvals.append (cyc[cycles[c_i]])
    #
    #         r_list = {}
    #         block = ''
    #         self.enumct = 0
    #         CE = CompleteEnum (outfile)
    #         pool_size = NUM_WORKERS
    #         pool = Pool(pool_size)
    #         for r in range (0, rndmsample_ct):
    #             pool.apply_async(EnumFullMolecule, args=(cycles, cycvals, r_list), callback=CE.CompleteEnumAsync)
    #         pool.close ()
    #         pool.join ()
    #         CE.WriteBlock()
    #
    #     outfile.close ()

    # def Deprotect(self, infile, deprotect_specfile, dp_outfile, smilescol, replace):
    #     if deprotect_specfile is None:
    #         return
    #     df = pd.read_csv(infile)
    #     dp_list = pd.DataFrame(columns=df.columns)
    #     for idx, row in df.iterrows():
    #         m = Chem.MolFromSmiles(row[smilescol])
    #         try:
    #             res, deprotected = ChemUtilities.Deprotect(m, deprotect_specfile, True)
    #             if deprotected == True:
    #                 if replace == True:
    #                     l = row[smilescol]
    #                     df.at[idx, smilescol] = Chem.MolToSmiles(res[0])
    #                 else:
    #                     r2 = row.copy()
    #                     r2[smilescol] = Chem.MolToSmiles(res[0])
    #                     dp_list = dp_list.append(r2)
    #         except:
    #             deprotected = False
    #             res = [m]
    #     if replace:
    #         df.to_csv(dp_outfile)
    #     else:
    #         df = df.append(dp_list)
    #         df.to_csv(dp_outfile)

    # def FilterBBs(self, bbdict, filterfile):
    #     filtersdf = pd.read_csv(filterfile)
    #     patterndict = {}
    #     removeidxs = {}
    #     for ix, r in filtersdf.iterrows():
    #         pattern = r['smarts']
    #         patterndict[pattern] = Chem.MolFromSmarts(pattern)
    #
    #     for dx in bbdict.keys():
    #         removeidxs[dx] = []
    #         for idx, row in bbdict[dx].iterrows():
    #             try:
    #                 m = Chem.MolFromSmiles(row[self.smilescol])
    #             except Exception as e:
    #                 print(filterfile)
    #                 print(row[self.idcol], row[self.smilescol])
    #                 raise (e)
    #             for v in patterndict.values():
    #                 if m.HasSubstructMatch(v) == True:
    #                     removeidxs[dx].append(idx)
    #                     continue
    #         bbdict[dx].drop(removeidxs[dx], axis=0, inplace=True)
    #         bbdict[dx] = bbdict[dx].reset_index()
    #     return bbdict

    # def Get_MoleculesFromSMARTSFilters(self, infile, outfile, inclSMARTS, exclSMARTS):
    #     df = pd.read_csv(infile)
    #     keep_idxs = []
    #     filters_dict = {}
    #     filtername = 'f1'
    #     filters_dict[filtername] = {}
    #     if inclSMARTS is None:
    #         inclSMARTS = []
    #     filters_dict[filtername]['include'] = inclSMARTS
    #     if exclSMARTS is None:
    #         exclSMARTS = []
    #     filters_dict[filtername]['exclude'] = exclSMARTS
    #     for idx, row in df.iterrows():
    #         input_mol = Chem.MolFromSMILES(row[self.smilescol])
    #         if (self.PassFilters(input_mol, filters_dict, filtername)):
    #             keep_idxs.append(idx)
    #     df = df.iloc[keep_idxs].reset_index()
    #     print(df)
    #     df.to_csv(outfile)

    # def generate_BBlists(self, inpath, outpath, libname, rxnschemefile):
    #     f = open(rxnschemefile)
    #     data = json.load(f)
    #     f.close()
    #     names_dict = data[libname]['filters']['names']
    #     listnames = names_dict.keys()
    #     filters_dict = data[libname]['filters']['BB_filters']
    #     if 'GeneralFilters' in data:
    #         GenFilters_dict = data['GeneralFilters']
    #     else:
    #         GenFilters_dict = {}
    #
    #     for n in names_dict:
    #         if n not in filters_dict:
    #             if n not in GenFilters_dict:
    #                 print('Error: Filter not found')
    #             else:
    #                 filters_dict[n] = GenFilters_dict[n]
    #
    #     print('pull bbs')
    #     df = self.pull_BBs(inpath)
    #     print('pull bbs complete')
    #     for ln in listnames:
    #         df[ln] = None
    #
    #     saltstrip = SaltRemover.SaltRemover()
    #     print('loop start')
    #     for idx, row in df.iterrows():
    #         try:
    #             input_mol = Chem.MolFromSmiles(row[self.smilescol])
    #             input_mol = saltstrip.StripMol(input_mol)
    #             df.iloc[idx][self.smilescol] = Chem.MolToSmiles(input_mol)
    #         except:
    #             input_mol = None
    #         if input_mol is not None:
    #             for ln in listnames:
    #                 if self.PassFilters(input_mol, filters_dict, ln):
    #                     df.iloc[idx][ln] = 'Y'
    #                 else:
    #                     df.iloc[idx][ln] = 'N'
    #     print('loop done')
    #     if not os.path.exists(outpath):
    #         os.makedirs(outpath)
    #
    #     for ln in listnames:
    #         cdf = df[df[ln] == 'Y'].drop_duplicates(subset=[self.idcol]).reset_index(drop=True)
    #         if type(names_dict[ln]) == list:
    #             for ix in range(0, len(names_dict[ln])):
    #                 cdf.to_csv(outpath + '/' + libname + '.' + names_dict[ln][ix] + '.csv', index=False)
    #         else:
    #             cdf.to_csv(outpath + '/' + libname + '.' + names_dict[ln] + '.csv', index=False)
    # def Enumerate_Dummy_Scaffold(self, rxnschemefile, schemename, bbsmiles, rxtntnum, in_reactions=None):
    #     scheme, rxtnts = self.ReadRxnScheme(rxnschemefile, schemename, FullInfo=True)
    #     if 'scaffold_dummy_structures' not in scheme:
    #         return 'FAIL'
    #     else:
    #         inrxtnts = scheme['scaffold_dummy_structures']
    #         if bbsmiles is not None and rxtntnum is not None:
    #             inrxtnts[rxtntnum] = bbsmiles
    #         p, prod_ct, [scheme, rxtants, intermeds] = self.RunRxnScheme(inrxtnts, rxnschemefile, schemename, False,
    #                                                                      in_reactions=None)
    #         return p
    # def Get_CycleList(self, rxnschemefile, schemename):
    #     scheme, rxtnts = self.ReadRxnScheme(rxnschemefile, schemename, verbose=True, FullInfo=True)
    #     return scheme['BB_Cycles']
    #
    #
    # def Get_LibCycle_BBCriteria(self, rxnschemefile, schemename, cycle, fieldname='filters'):
    #     scheme, rxtnts = self.ReadRxnScheme(rxnschemefile, schemename, FullInfo=True)
    #     if fieldname in scheme:
    #         if cycle in scheme[fieldname]:
    #             return scheme[fieldname][cycle]
    #         else:
    #             if 'names' in scheme[fieldname]:
    #                 for n in scheme[fieldname]['names']:
    #                     if scheme[fieldname]['names'][n] == cycle:
    #                         return scheme[fieldname]['BB_filters'][n]
    #                 return None
    #     else:
    #         return None
    # def PassFilters(self, input_mol, filters_dict, filtname):
    #     filters = filters_dict[filtname]
    #     incl_list = filters['include']
    #     if 'exclude' in filters:
    #         excl_list = filters['exclude']
    #     else:
    #         excl_list = []
    #     for ifl in incl_list:
    #         if type(ifl) == list:
    #             hassub = False
    #             for orx in ifl:
    #                 pattern = Chem.MolFromSmarts(orx)
    #                 if (input_mol.HasSubstructMatch(pattern) == True):
    #                     hassub = True
    #             if hassub == False:
    #                 return False
    #         else:
    #             pattern = Chem.MolFromSmarts(ifl)
    #             if (input_mol.HasSubstructMatch(pattern) == False):
    #                 return False
    #     for efl in excl_list:
    #         pattern = Chem.MolFromSmarts(efl)
    #         if (input_mol.HasSubstructMatch(pattern) == True):
    #             return False
    #     return True
    # def Read_ReactionDatabase(self, rxnfile):
    #     rxndb = toml.load(rxnfile)
    #     self.rxndict = {}
    #     for rtype in rxndb:
    #         subtypes = rxndb[rtype]
    #         for subtype in subtypes:
    #             self.rxndict[subtype] = subtypes[subtype]['smarts']
