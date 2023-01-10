from re import L
import requests
import time
import os
import zlib
import json
import sys
from tqdm import tqdm
import mongo_api as mapi
import multiprocessing

def pause():
    input('@@@ Press enter to continue: ')

#fetches all tranches
def fetch_zinc(rep="3D", since="", db_r="", format="smi",using="uri", DB_NAME = 'Zinc15', COLLECTION='Zinc15 Universe'):
    all_3d_url = "https://zinc.docking.org/tranches/all3D.json"
    download_url = "https://zinc.docking.org/tranches/download"
    r = requests.get(all_3d_url)
    tranches = r.json()

    data = {
        'representation': rep,
        'tranches': ' '.join(x['name'] for x in tranches),
        'format': format,
        'using': using
    }

    r = requests.post(download_url, data=data, stream=True)
    
    #multiprocessing the download and then set up of local database
    #print(r.text)
    process_1 = multiprocessing.Process(target=open_mol2, args=(r, DB_NAME, COLLECTION))
    process_1.start()

    #open_mol2(r)
    #mapi.grab_smiles('Zinc15', 'Zinc15 Universe')
    

#takes list of tranches and treats them individually
def open_mol2(r, DB_NAME, COLLECTION, WORKING_DIR='./tranches/'):

    os.makedirs(WORKING_DIR, exist_ok=True)

    #the '0_' is so the 2 files line up at the top of the 

    with open(WORKING_DIR + '0_alltranches', 'w') as a:
        total = r.text
        a.write(total)
        total = len(str(total).split('\n'))
    
    try:
        with open(WORKING_DIR + '0_finished_tranches', 'r') as all_file:
            done_list = set(all_file.read().split())
    except FileNotFoundError:
        done_list = set()

    print('Downloading tranches')

    process_2 = multiprocessing.Process(target=mapi.grab_smiles, args=(DB_NAME, COLLECTION))
    process_2.start()

    for x in tqdm(r.iter_lines(), total=total):

        if (str(x)) not in done_list:

            #this bit is only to visually confirm current tranche in the terminal
            #sys.stdout.write("\033[K")
            #print('Current tranche: ', x[36:42],
            #'Status: downloading', end='\r')

            #getting each tranche and writing the gz file
            resp = requests.get(x)

            #apparently zinc15 has tons of empty tranches, so this conditional is necessary
            #if server is down then this will assume all tranches are empty
            if(resp.ok):
                current_file_name = os.path.join(WORKING_DIR, x[36:].decode('utf-8'))
                with open (current_file_name, 'wb+') as f:
                    f.write(resp.content)

                #gz_files()
                #zlib test, doesn't work, use gz
                '''with open(current_file_name, 'rb') as f:
                    current_tranche = zlib.decompress(f.read(), zlib.MAX_WBITS|16)
                current_tranche = current_tranche.decode('utf-8').split('\n')
                
                with open ('tmpfile', 'w+') as e:
                    json.dump(current_tranche, e, indent = 4)'''

            #writes to a file with a list of correctly downloaded tranches
            with open(WORKING_DIR + '0_finished_tranches', 'a') as f:
                f.write(str(x) + '\n')
        
        #else:
            #print('Current tranche: ', x[36:42],
            #'Status: already downloaded', end="\r")
    
    done_list.close()

#future function for gz file handling
#def gz_files(current_file):
