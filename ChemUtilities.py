import pandas as pd
import tqdm
from rdkit import Chem
from rdkit.Chem import SaltRemover
import Enumerate

def SaltStripMolecules (molecules: pd.DataFrame, smilescol='SMILES', neutralize = True):
    tqdm.pandas ()
    molecules[smilescol] = molecules[smilescol].progress_apply(lambda smi:SaltStrip(smi, neutralize))
    return molecules

def SaltStrip (molec, neutralize = True):
    try:
        saltstrip = SaltRemover.SaltRemover()
        if type(molec) != Chem.Mol:
            m = Chem.MolFromSmiles(molec)
        else:
            m = molec
        m=  saltstrip.StripMol(m)
        smi = Chem.MolToSmiles(m, kekuleSmiles = False)
        smi = Neutralize(smi)
    except:
        return molec
    return smi

def Neutralize (smi):
    try:
        m = Chem.MolFromSmiles(smi)
    except:
        return smi
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = m.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    try:
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = m.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()

        smi = Chem.MolToSmiles(m, kekuleSmiles=False)
        return smi
    except:
        return smi

def Deprotect (compound, rxnscheme, deprotect_name=None, iterate = False, retsmiles = False):
    deprotected = False
    Enumerator = Enumerate.Enumerate ()
    first_round = True

    has_scheme = Enumerator.ReadRxnScheme(rxnscheme, 'Deprotection',  verbose=True)
    schemeinfo = [[{"Reactants":["r0","None"],"Rxns":{"default":[]}}], ["r0","None"] ]
    if deprotect_name is None:
        for k,v in has_scheme[0][0]['Rxns'].items ():
            if type(v) is list:
                for l in v:
                    schemeinfo[0][0]['Rxns']["default"].append(l)
            else:
                schemeinfo[0][0]['Rxns']["default"].append(v)
    else:
        if type (deprotect_name) is str:
            deprotect_name = [deprotect_name]
        for d in deprotect_name:
            v = has_scheme[0][0]['Rxns'] [d]
            if type(v) is list:
                for l in v:
                    schemeinfo[0][0]['Rxns']["default"].append(l)

    if not has_scheme[0] is  None:
        last_dp = False
        while first_round or  (iterate and last_dp == True):
            res, prod_ct, resinfo = Enumerator.RunRxnScheme([compound, 'None'], None,None, False, schemeinfo )
            if (res != 'FAIL' and res != 'None'):
                if type(res) == str:
                    compound = Chem.MolFromSmiles(res)
                deprotected = True
                last_dp = True
            else:
                last_dp = False

            first_round = False

        first_round = True
        last_dp = False
        res_list = [compound]

    has_scheme = Enumerator.ReadRxnScheme(rxnscheme, 'Deprotect_KeepBoth', False)
    if has_scheme[0] is not None:
        while first_round or (iterate and last_dp == True):
            res, prod_ct = Enumerator.RunRxnScheme([compound, 'None'], rxnscheme, 'Deprotect_KeepBoth', False)
            if res == 'NOSCHEMEFAIL':
                break
            if (res != 'NOSCHEMEFAIL' and res != 'FAIL' and res != 'None'):
                res_list.append (Chem.MolFromSmiles(res))
                compound = Chem.MolFromSmiles(res)
                deprotected = True
                last_dp = True
            else:
                last_dp = False
            first_round = False


    res = res_list
    if retsmiles:
        if deprotected:
            return Chem.MolToSmiles(res[0])
        else:
            return compound
    return res, deprotected