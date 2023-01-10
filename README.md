# auto-gaussian
A variety of chemical tools to facilitate the use of Gaussian


<div align="center">

# auto-gaussian

auto-gaussian is a small set of tools designed to facilitate certain repetitive tasks related to [Gaussian][Gaussian].

It handles downloads from certain online databases, makes local database configuration easy, and parses inputs & outputs to/from Gaussian.<br />

[Getting Started](#getting-started) •
[Usage](#usage)

</div>

## Getting started

First, make sure you have [Docker][Docker] installed and configured on your machine. These scripts are configured to work exclusively with Docker because the constant setting up of rdkit, mongodb, etc. workflows was a pain, especially for chemists with minimal experience with these tools. 

The .yml file is included, so simply navigate to the local directory where you cloned this repo and run:

```bash 
docker-compose up -d
```

This will create 3 containers; one for running python scripts and two for isolated databases. Please note that the parent directory of the .yml file will be mounted as a volume on the python container.

On some older distributions you might need to change the ‘version’ in the .yml file to an older, compatible version. I’m using ver 3.9.

Now find the mongo_db container’s ip address using:

```bash 
docker inspect mongo_db
```

Near the end of the information printed you should see something like:

```text
                    "IPAddress": "172.21.0.4",
                    "IPPrefixLen": 16,
                    "IPv6Gateway": "",
                    "GlobalIPv6Address": "",
                    "GlobalIPv6PrefixLen": 0,
                    "MacAddress": "02:42:ac:15:00:04",
                    "DriverOpts": null
                }
            }
        }
    }
]
```

Copy the “IPAddress” field and paste it into the mongo_api.py global variable DOMAIN so it looks something like this:

```python
DOMAIN = '172.21.0.4'
```

Now you should be able to execute an interactive shell inside the python container using:

```bash
docker exec -it coding_area bash
```


## Usage

Some basic functionalities are exemplified here.

[Zinc15 database download](#zinc15-database-download) •
[Filtering](#filtering) •
[Deleting](#deleting-databases-and/or-collections) •

## Zinc15 database download

First, let’s visualize the contents of all MongoDB databases out-of-the-box by creating a local file **example.py** in the same directory as the cloned files:

```python
import mongo_api as mapi

mapi.visualize_contents()
```

#### Output

```bash
python3.py example.py
```
```text

'name': 'admin', 'sizeOnDisk': 102400, 'empty': False
'name': 'config', 'sizeOnDisk': 73728, 'empty': False
'name': 'local', 'sizeOnDisk': 73728, 'empty': False
```

Now I want to clone the entire Zinc15 database to my local hard drive, so I’ll use the fetch_zinc() function:

```python
import database_dl as dd

dd.fetch_zinc()
```

This will send a request to the Zinc15 database and automatically download & organize all 3D molecules as smiles strings in a local MongoDB dabatase. Now let’s visualize all databases again: 

```python
import mongo_api as mapi

mapi.visualize_contents()
```

#### Output

```text
'name': 'Zinc15', 'sizeOnDisk': XXXXX, 'empty': False
'name': 'admin', 'sizeOnDisk': 102400, 'empty': False
'name': 'config', 'sizeOnDisk': 110592, 'empty': False
'name': 'local', 'sizeOnDisk': 73728, 'empty': False
```

I can also view the collections (tables in MySQL) inside of a specific database by adding an optional argument – the name of the collection:

```python
import mongo_api as mapi

mapi.visualize_contents('Zinc15')
```

#### Output

```text
Collections in Zinc15:

Zinc15 Universe
```

Right now there’s only one collection. If we want to view how many elements are in said collection, we can add another optional argument, the collection name:

```python
import mongo_api as mapi

mapi.visualize_contents('Zinc15', 'Zinc15 Universe')
```

#### Output

```text
Total number of documents in (Zinc15 Universe): XXXXX
```

The optional arguments are nested, so it doesn’t make sense to call a collection without also calling the parent database. Now, when calling a collection the function actually returns said collection as a pandas dataframe for ease of use, so you could do something like:

```python
import mongo_api as mapi

collection = mapi.visualize_contents('Zinc15', 'Zinc15 Universe')
print(collection)
```

#### Output

```text
Total number of documents in (Zinc15 Universe): 84116

                                        smiles_chain    zinc_id     tranche   maxdist atomspresent
0      N[C@@H](C(=O)[O-])[C@@H](O)[C@@H](O)C(=O)[O-]   14419527  AAAEHL.smi  7.516452    [7, 6, 8]
1                   Nc1nc(N)c(CCC(=O)[O-])c(=O)[n-]1  214078112  AAAEHL.smi  9.219045    [7, 6, 8]
2                 Nc1[nH]c(=O)[n-]c(=O)c1CCC(=O)[O-]  214141302  AAAEHL.smi  9.054797    [7, 6, 8]
3       N[C@@H](C(=O)[O-])[C@@H](O)[C@H](O)C(=O)[O-]   13351315  AAAEHL.smi  7.516452    [7, 6, 8]
4       C[C@@H](C(=O)[O-])[C@@H](O)[C@H](N)C(=O)[O-]   13353484  AAAEHL.smi  7.516452    [6, 8, 7]
...                                              ...        ...         ...       ...          ...
84111          NC(=O)[C@H]1CCC[N@@H+]1C1CC([NH3+])C1   87653181  ABAAMP.smi  8.491547    [7, 6, 8]
84112           NC(=O)[C@H]1CCC[N@H+]1C1CC([NH3+])C1   87653181  ABAAMP.smi  8.491547    [7, 6, 8]
84113                      CC(C)([NH3+])C1=[NH+]CCN1  104653403  ABAAMP.smi  4.794812       [6, 7]
84114         NC(=O)[C@@H]1CCC[N@@H+]1C1CC([NH3+])C1   87653182  ABAAMP.smi  8.491547    [7, 6, 8]
84115          NC(=O)[C@@H]1CCC[N@H+]1C1CC([NH3+])C1   87653182  ABAAMP.smi  8.491547    [7, 6, 8]

[84116 rows x 5 columns]
```

In reality the number will be much higher than 84116, of course. I also added a few commonly utilized values, such as the atoms contained in each molecule (by atomic number) and a **very** rough approximation of molecule size (simply the longest vector distance between any 2 atoms in the 3D molecule in Å). The zinc_id and tranche columns are simply indicators native to the Zinc15 database. 

## Filtering

These filters start from any collection (chemical space or subspace) and create a new collection (subspace). They can be mixed and matched per the specific needs of the project.
[By size](#filter-by-approximated-size) •
[By atoms present](#filter-by-atoms-present) •
[By substructure present](#filter-by-substructure) •

#### Filter by approximated size

There are three available filtering functions, the first of which is a size filter (keep in mind the size calculations done by default are very approximated). Please note you must specify which collection you are filtering from. Say I want only the Zinc15 molecules between (3 – 7) Å, I would filter with:

```python
import mongo_api as mapi1

mapi.sizefilter('Zinc15', 'Zinc15 Universe', min=3, max=7)
mapi.visualize_contents('Zinc15')
```

#### Output

```text
Collections in Zinc15:

Size filter from (Zinc15 Universe): > 3 ; < 7  20230110-014247
Zinc15 Universe
```
Now a new collection within the db was created, where the numbers are the date in YYYYMMDD-hhmmss format. You could also specify the collection name, if you wanted, with the optional argument:

```python
import mongo_api as mapi

mapi.sizefilter('Zinc15', 'Zinc15 Universe', min=3, max=7, newname = 'Test01')
collection = mapi.visualize_contents('Zinc15', 'Test01')
print(collection)
```

#### Output

```text
Total number of molecules that fit the size criteria:  16616

Total number of documents in (Test01): 16616
                                   smiles_chain    zinc_id     tranche   maxdist   atomspresent
0         O=c1cc(P(=O)([O-])[O-])[nH]c(=O)[nH]1    1604752  AAAEHL.smi  6.449067  [8, 6, 15, 7]
1            O=C([O-])C[C@@H](O)P(=O)([O-])[O-]  140177892  AAAEHL.smi  6.660458     [8, 6, 15]
2                 N=C1C(=O)C(=O)[N-]C(=O)N1[O-]   39059935  AAAEHL.smi  5.588478      [7, 6, 8]
3      O=c1[n-]c(=O)c2c(=O)[n-]c(=O)[nH]c2[nH]1    5670124  AAAEHL.smi  6.958927      [8, 6, 7]
4                      Cn1nc([S-])c(=O)[n-]c1=O    5893136  AAAEHL.smi  5.844241  [6, 7, 16, 8]
...                                         ...        ...         ...       ...            ...
16611         CC(C)[C@@H]([NH3+])C1(O)C[NH2+]C1  958528614  ABAAMP.smi  5.947813      [6, 7, 8]
16612      CC(C)(C)[C@@H]([NH3+])C1(O)C[NH2+]C1  958528467  ABAAMP.smi  5.947813      [6, 7, 8]
16613         [NH3+][C@@H](C1CC1)C1(O)C[NH2+]C1  958528617  ABAAMP.smi  6.006393      [7, 6, 8]
16614            C[N@H+]1CC[C@H](C(N)=[NH2+])C1   84298158  ABAAMP.smi  6.201146         [6, 7]
16615                 CC(C)([NH3+])C1=[NH+]CCN1  104653403  ABAAMP.smi  4.794812         [6, 7]

[16616 rows x 5 columns]
```

Thus reducing the chemical space from ~84k molecules to ~17k.

Also note that if you use the same collection name twice, the data will not be overwritten, but rather appended onto the end of the collection, so accidentally executing the same command twice with the same custom name will result in a collection where the entire second half is a duplicate of the first half.

#### Filter by atoms present

This time I want to filter not from the entire Zinc15 Universe but rather from the ‘Test01’ subspace I created before, and I’d like only molecules that contain Carbon, Nitrogen, and Oxygen. For this I’d pass my new collection name (optional) and a tuple containing the atomic numbers I want to the following function:

```python
import mongo_api as mapi

mapi.atomfilter('Zinc15', 'Test01', contains = (6,7,8), newname = 'Test02')
collection = mapi.visualize_contents('Zinc15', 'Test02')
print(collection)
```

#### Output

```text
13187 values met the criteria

Total number of documents in (Test02): 13187
                                   smiles_chain    zinc_id     tranche   maxdist atomspresent
0                 N=C1C(=O)C(=O)[N-]C(=O)N1[O-]   39059935  AAAEHL.smi  5.588478    [7, 6, 8]
1      O=c1[n-]c(=O)c2c(=O)[n-]c(=O)[nH]c2[nH]1    5670124  AAAEHL.smi  6.958927    [8, 6, 7]
2                 O=c1[n-]c(=O)n([O-])c2nccnc12    6119146  AAAEHL.smi  6.303920    [8, 6, 7]
3               [O-]N=c1nc[n-]c2no[n+]([O-])c12  687727564  AAAEHL.smi  5.821749    [8, 7, 6]
4                Nc1c(NC=O)c(=O)[n-]c(=O)n1[O-]  687728024  AAAEHL.smi  6.861276    [7, 6, 8]
...                                         ...        ...         ...       ...          ...
13182         CC(C)[C@@H]([NH3+])C1(O)C[NH2+]C1  958528614  ABAAMP.smi  5.947813    [6, 7, 8]
13183      CC(C)(C)[C@@H]([NH3+])C1(O)C[NH2+]C1  958528467  ABAAMP.smi  5.947813    [6, 7, 8]
13184         [NH3+][C@@H](C1CC1)C1(O)C[NH2+]C1  958528617  ABAAMP.smi  6.006393    [7, 6, 8]
13185            C[N@H+]1CC[C@H](C(N)=[NH2+])C1   84298158  ABAAMP.smi  6.201146       [6, 7]
13186                 CC(C)([NH3+])C1=[NH+]CCN1  104653403  ABAAMP.smi  4.794812       [6, 7]

[13187 rows x 5 columns]
```

Further reducing the space to ~13k molecules.

Notice the molecules don’t have to contain all 3 atoms I specified, but they cannot contain any others.

#### Filter by substructure

This filter takes advantage of RDKit’s HasSubstructMatch function, so it isn’t perfect, but it’s a pretty useful tool nonetheless. Let’s say I want to further filter ‘Test02’ and create a new collection with only the molecules that contain the carbonyl and propyl substructures, I would need to pass a dataframe containing a dictionary-like structure of n molecule-smile pairs:

```python
import mongo_api as mapi

mapi.subsearch('Zinc15', 'Test02', substructs = pd.DataFrame(
        {
            'name': ['carbonyl', 'propyl'],
            'smiles': ['C=O','CCC']
        }),
        newname = 'Test03')
collection = mapi.visualize_contents('Zinc15', 'Test03')
print(collection)
```

#### Output

```text
12104 molecules match the filters

Total number of documents in (Test03): 12104
         zinc_id                              smiles_chain  substructure match
0       39059935             N=C1C(=O)C(=O)[N-]C(=O)N1[O-]  [carbonyl, propyl]
1        5670124  O=c1[n-]c(=O)c2c(=O)[n-]c(=O)[nH]c2[nH]1          [carbonyl]
2        6119146             O=c1[n-]c(=O)n([O-])c2nccnc12          [carbonyl]
3      687728024            Nc1c(NC=O)c(=O)[n-]c(=O)n1[O-]          [carbonyl]
4       17130245                Nn1[n-]c([O-])c2ncnc-2c1=O          [carbonyl]
...          ...                                       ...                 ...
12099  958528614         CC(C)[C@@H]([NH3+])C1(O)C[NH2+]C1            [propyl]
12100  958528467      CC(C)(C)[C@@H]([NH3+])C1(O)C[NH2+]C1            [propyl]
12101  958528617         [NH3+][C@@H](C1CC1)C1(O)C[NH2+]C1            [propyl]
12102   84298158            C[N@H+]1CC[C@H](C(N)=[NH2+])C1            [propyl]
12103  104653403                 CC(C)([NH3+])C1=[NH+]CCN1            [propyl]

[12104 rows x 3 columns]
```

Thus creating an even more reduced subspace. This substructure filter should have no trouble filtering by any amount of substructures at once, although it’s not recommended to apply these kinds of filters to super large datasets, as it can get really time consuming.

## Deleting databases and/or collections

If we visualize the contents of the ‘Zinc15’ database at this point, it would look something like this:

```python
import mongo_api as mapi

mapi.visualize_contents('Zinc15')
```

#### Output

```text
Collections in Zinc15:

Test02
Test01
Test03
Zinc15 Universe
```

If I wanted to delete the Test01 collection, I could use the following function:

```python
import mongo_api as mapi

mapi.rm_collection('Zinc15', 'Test03')
mapi.visualize_contents('Zinc15')
```

#### Output

```text
Collections in Zinc15:

Test02
Test01
Zinc15 Universe
```

However, if I prefer a nuclear option instead (deleting an entire database), I could use:

```python
import mongo_api as mapi

mapi.rm_db('Zinc15')
```

#### Output

```text
'name': 'admin', 'sizeOnDisk': 102400, 'empty': False
'name': 'config', 'sizeOnDisk': 110592, 'empty': False
'name': 'local', 'sizeOnDisk': 73728, 'empty': False
```

[Gaussian]: https://gaussian.com/
[Docker]: https://docs.docker.com/get-docker/