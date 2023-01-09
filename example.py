import database_dl as dd
import mongo_api as mapi
import smiles_gaussian

# download from zinc15
#dd.fetch_zinc()

# grab downloaded files from zinc15 and send them to db
#mapi.rm_db('Zinc15')
#mapi.grab_smiles('Zinc15', 'Zinc15 Universe')

# write smiles file from specific collection
#mapi.to_gaussinput('Zinc15', 'subsearch 20221213-183536')

# run transformation from 
smiles_gaussian.gauss_gen()

# delete a DB or collection
#mapi.rm_db('Zinc15')
#mapi.rm_collection('Zinc15', 'test')

# query - not working
#mapi.query('dna','smiles 2022-11-04')

# size filter
#mapi.sizefilter('Zinc15', 'Zinc15 Universe', min=3, max=7)

# atom filtering
#mapi.atomfilter('Zinc15', 'Size filter from (Zinc15 Universe): > 3 ; < 7  20221213-181709')

# substructure search
#mapi.subsearch('Zinc15', 'Zinc15 Universe')

# visualization
collection = mapi.visualize_contents('Zinc15')
#print(collection)




