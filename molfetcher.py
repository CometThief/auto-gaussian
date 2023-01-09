import pubchempy as pcp
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
import pandas as pd

def getdf(smiles, directory):
    #m = pcp.get_cids('2-nonenal', 'name', 'substance', list_return='flat')
    #pcp.get_compounds('CC', searchtype='superstructure', listkey_count=3)
    #C1=C(C2=CC(=C(C=C2OC1=O)O)Br)CO
    #C1=CC=C2C=CC=CC2=C1

    #m = pcp.get_compounds('C1=CC=C2C=CC=CC2=C1', namespace = 'smiles',
    #searchtype = 'substructure')

    #m = pcp.get_cids('CC', namespace = 'smiles', searchtype='superstructure', 
    #listkey_count=3)

    if smiles[-2:] == '\n':
        smiles = smiles[:-2]

    mol = Chem.MolFromSmiles(smiles)

    # draw the modecule -----------------------------------------------------------
    #Draw.MolToFile(mol, 'molecule.png')

    # draw the molecule with property ---------------------------------------------
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    
    picdirectory = './gauss_files/output/' + directory + '/'
    os.makedirs(picdirectory, exist_ok=True)
    mol = Chem.AddHs(mol)
    Draw.MolToFile(mol, picdirectory + smiles + '.png')

    #print(mol.GetAtomWithIdx(3).GetIsAromatic())
    #print(mol.GetBondBetweenAtoms(3,6).GetIsAromatic())

    # print the atoms of the molecule ---------------------------------------------
    isaromatic = []
    for atom in mol.GetAtoms():
        isaromatic.append(atom.GetIsAromatic())
        '''
        print(atom.GetIdx(),',',
            atom.GetAtomicNum(),',',
            atom.GetIsAromatic(),',',
            atom.GetSymbol())
        '''

    # print the bonds of the molecule ---------------------------------------------
    '''
    for bond in mol.GetBonds():
        print(bond.GetBeginAtomIdx(),',',
            bond.GetEndAtomIdx(),',',
            bond.GetBondType())
    '''

    

    AllChem.EmbedMolecule(mol,randomSeed=0xf00d)
    w = Chem.MolToMolBlock(mol)

    #with open('Output.sdf', 'w') as file:
        #file.write(w)

    w = w.split('\n')

    templist = [] 

    #if i want to change the try/except stuff later down, use numatoms
    #numatoms = int(w[3].split()[0])

    for i in range(4,len(w)-1):
        if not w[i].startswith('M'):
            templist.append(w[i].split()[:4])

    coordsect = []
    bondsect = []
    nonHs = 0
    for i in templist:
        try:
            temp = [int(x) for x in i]
            bondsect.append(temp)
        except:
            temp = [float(x) for x in i[:-1]]
            temp.append(i[-1])
            coordsect.append(temp)
            if i[-1] != 'H':
                nonHs += 1

    final = pd.DataFrame(coordsect, columns=['x','y','z','atom'])
    for i in range(len(isaromatic), len(coordsect)):
        isaromatic.append(False)

    final['Aromaticity'] = isaromatic

    templist = pd.DataFrame(bondsect)
    bonding = []
    bondtypes = []

    for i in range(1, nonHs+1):
        col1 = templist.loc[templist[0] == i].index.tolist()

        oneline = []
        twoline =[]
        for i in col1:
            oneline.append(templist.iloc[i, 1])
            twoline.append(templist.iloc[i, 2])

        bonding.append(oneline)
        bondtypes.append(twoline)

    for i in range(nonHs+1, len(coordsect)+1):
        col2 = templist.loc[templist[1] == i].index.tolist()

        oneline = []
        twoline =[]
        for i in col2:
            oneline.append(templist.iloc[i, 0])
            twoline.append(1)

        bonding.append(oneline)
        bondtypes.append(twoline)

    final['Bonded To'] = bonding
    final['Bond Orders'] = bondtypes
    
    #print(templist)
    #print(final)
    return final

'''
smiles = 'C1=C(C)C=C2C=CN=C(S)C2=C1' 
a = getdf(smiles)
print(a)
'''