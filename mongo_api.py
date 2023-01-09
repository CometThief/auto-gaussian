import pymongo as pm
from pymongo import MongoClient
from datetime import date
import os
import json
import typing
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDistGeom as molDG
from tqdm.auto import tqdm
import numpy as np
import time
from collections import defaultdict
import sys

DOMAIN = '172.21.0.1'
PORT = 27017

client = MongoClient(
    host = [ str(DOMAIN) + ":" + str(PORT) ],
    serverSelectionTimeoutMS = 3000, # 3 second timeout
    username = "dna",
    password = "12345",
)

def insert(DB, COLLECTION, documents):
    DB = client[DB]
    COLLECTION = DB[COLLECTION]
    inserted = COLLECTION.insert_many(documents)
    return inserted

def rm_db(DB):
    client.drop_database(DB)

def rm_collection(DB, COLLECTION):
    DB = client[DB]
    COLLECTION = DB[COLLECTION]
    deleted = COLLECTION.drop()
    return deleted

def visualize_contents(DB = None, COLLECTION = None):
    
    col_df = None

    if DB != None:
        dbname = DB
        DB = client[DB]

    if (DB != None) and (COLLECTION != None):
        collection = DB[COLLECTION]
        all = collection.find({}, {'_id': False})
        col_df = pd.DataFrame(all)

        print('\nTotal number of documents in ({collection}): {amount}'.format(
            collection=COLLECTION, amount=collection.count_documents({})) )

    elif (DB != None) and (COLLECTION == None):
        print('\nCollections in {}:\n'.format(dbname))
        allcollections = DB.list_collection_names()
        for x in allcollections:
            print(x)
        print('\n')

    else:
        print('\n')
        for db in client.list_databases():
            print(str(db)[1:-1])
        print('\n')

    return col_df

def grab_smiles(DB_NAME, COLLECTION = 'smiles ' + time.strftime("%Y%m%d-%H%M%S"), maxemptyiterations=1000):

    #first we fetch all files in the download directory
    allfiles = set()
    to_remove = set()
    iterations = 0
    dir = './tranches/'
    
    #each file is a tranche, so we gotta iterate through them
    counter = 0
    batch = list()

    while(iterations <= maxemptyiterations):
        if not allfiles:
            for file in os.listdir(dir):
                if file.endswith(".smi"):
                    allfiles.add(os.path.join(dir, file))
            if not allfiles:
                print('The smiles directory is empty')
                iterations += 1
                time.sleep(5)
        else:
            for i in tqdm(allfiles):
                tranche = i[-10:]
                with open(i, 'r') as file:
                    for e in file:
                        #this should be generalized with regex to be able to handle more file types (?)
                        #only works with zinc15 smiles files as is
                        if not e.startswith('smiles'):
                            smiles_chain = e.split()[0]
                            zinc_id = e.split()[1]
                            mol = Chem.MolFromSmiles(smiles_chain)

                            #distance matrix stuff
                            distmatrix = molDG.GetMoleculeBoundsMatrix(mol)
                            maxdist = np.max(distmatrix)

                            #atoms present in molecule by atomic number
                            atomic_count = defaultdict(lambda : 0)
                            for atom in mol.GetAtoms():
                                atomic_count[atom.GetAtomicNum()] += 1
                            atomspresent = [i for i in atomic_count if atomic_count[i]!=atomic_count.default_factory()]

                            document = {
                                'smiles_chain' : smiles_chain,
                                'zinc_id' : zinc_id,
                                'tranche' : tranche,
                                'maxdist' : maxdist,
                                'atomspresent' : atomspresent
                            }
                            
                            batch.append(document)
                            counter += 1
                iterations = 0
                #print('Finished tranche: {ctranche}'.format(ctranche=i))
                os.remove(i)
                to_remove.add(i)
        
        if to_remove:
            for i in to_remove:
                allfiles.remove(i)
        to_remove = set()
                
        if batch:
            insert(DB_NAME, COLLECTION, batch)
            batch = list()

# out of order
def query(DB, COLLECTION, params = {'_id':False}):

    COLLECTION = pm.collection.Collection(makedb(DB), COLLECTION)

    for idx, x in enumerate(COLLECTION.find(projection = params)):
        print(idx + 1, x,'\n')

#def query():

def subsearch(DB_NAME, COLLECTION, substructs = pd.DataFrame(
        {
            'name': ['carbonyl'],
            'smiles': ['C=O']
        }
    ) , showall=False):

    newname = 'subsearch ' + time.strftime("%Y%m%d-%H%M%S")
    substructs['rdkit molecule'] = substructs['smiles'].apply(Chem.MolFromSmiles)

    #mols = list(substructs['rdkit molecule'])
    #names = list(substructs['name'])
    #output image with substructures?
    #Chem.Draw.MolsToGridImage(
    # mols,
    # legends=name
    # molsPerRow=5
    # )

    col = visualize_contents(DB_NAME, COLLECTION)
    matches = []
    #no_matches = []

    for index, row in tqdm(col.iterrows(), total=col.shape[0]):
        mol = Chem.MolFromSmiles(row['smiles_chain'])
        match = False
        for _, substruct in substructs.iterrows():
            if mol.HasSubstructMatch(substruct['rdkit molecule']):
                matches.append(
                    {
                        'zinc_id': row['zinc_id'],
                        'smiles': row['smiles_chain'],
                        'substructure match': substruct['name']
                    }
                )
                match = True
        #if not match:
            #no_matches.append(index)

    if matches and showall==True:
        with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
            print(matches)
    elif matches:
        insert(DB_NAME, newname, matches)
        matches=pd.DataFrame(matches)
        #no_matches=pd.DataFrame(no_matches)
        print('{0} molecules match the filters'.format(len(matches)))
        #print('{0} molecules dont'.format(len(no_matches)))
    else:
        print('\nNo values matched the filters\n')


def sizefilter(DB_NAME, COLLECTION, min = 3, max = 7):

    db = client[DB_NAME]
    collection = db[COLLECTION]
    filtered = []
    newname = 'Size filter from ({collection}): > {min} ; < {max}  {date}'.format(
        collection=COLLECTION, min=min, max=max, date=time.strftime("%Y%m%d-%H%M%S"))

    filtercursor = collection.find( {'$and': [
        { 'maxdist': { '$gt': min } },
        { 'maxdist': { '$lt': max } } 
        ] } )
    
    for x in filtercursor:
        filtered.append(x)
    
    insert(DB_NAME, newname, filtered)
    print('Total number of molecules that fit the size criteria: ', len(filtered))


def atomfilter(DB_NAME, COLLECTION, contains=(6,7,8)):
    #, exclude=range(min,max)

    db = client[DB_NAME]
    collection = db[COLLECTION]
    newcollection = '{contains} Atom filter from ({COLLECTION})  {DATE}'.format(
        contains=contains, COLLECTION=COLLECTION, DATE=time.strftime("%Y%m%d-%H%M%S")
    )

    print('\nFiltering from {DB_NAME} - {COLLECTION}:'.format(
        DB_NAME=DB_NAME, COLLECTION=COLLECTION
    ))

    allvalues = collection.find()
    results = []

    for value in allvalues:
        atoms = value['atomspresent']
        check = list(set(atoms) - set(contains))
        
        if not check:
            results.append(value)

    if results:
        print('{length} values met the criteria'.format(length = len(results)))
        insert(DB_NAME, newcollection, results)
    else: 
        print('\nNo values met the criteria')



def to_gaussinput(DB_NAME, COLLECTION):

    db = client[DB_NAME]
    collection = db[COLLECTION]
    allvalues = pd.DataFrame(collection.find())
    directory = './gauss_files/input/input.smi'
    content = '\n'.join((allvalues['smiles'].tolist()))

    with open(directory, 'w') as f:
        f.write(content)

DB = 'dna'
collname = 'smiles ' + str(date.today())

#db = client[DB]
#print(db)
#print(db.collname.count_documents({}))
#print(db.collname.count_documents(typing.Mapping[str, str]))
#input()
#delete_collection(DB, collname)
#grab_smiles(DB)
#visualize_contents(DB, collname, showall = True)


'''
cProfile.run('main(DB, collname)', 'restats')

import pstats
from pstats import SortKey
p = pstats.Stats('restats')
p.strip_dirs().sort_stats(-1).print_stats()
p.sort_stats(SortKey.NAME)
p.print_stats()

'''



#grab_smiles(DB)
#visualize_contents(DB, collname)
#delete(DB, collname, {})
#visualize_contents(DB, collname)

#client = pm.MongoClient("localhost", 27017, maxPoolSize=50)

