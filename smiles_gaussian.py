import numpy as np
import pandas as pd
import glob
import itertools
import time
from datetime import date, datetime
import os
import molfetcher as mf
import sys
import copy
from scipy.spatial.transform import Rotation as Rot
import math
import shutil
#import mysql_api as mysql_api
import mongo_api as mongo_api

DB_NAME = 'dna'

globalcounter = 0
coords = ['x', 'y', 'z']
vector = ['i', 'j', 'k']

atom46 = pd.Series(data=[-5.6144100, 1.7503100, -0.9838330], index=coords)
atom46b = pd.Series(data=[-5.0909910, 0.8250507, -1.1056020], index=coords)
atom128 = pd.Series(data=[5.1424560, 1.9134080, 0.3529960], index=coords)
atom128b = pd.Series(data=[4.6285837, 0.9950496, 0.5465239], index=coords)

def pause():
    input('\nPress enter to continue @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ')

def getAngle(a, b, c):
    '''angle between 3 points, b should be the anchor'''
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return ang + 360 if ang < 0 else ang

def fix_shape(carbons, input, index, inputvec, refvec, refcoord, refcoordpent, inputt, ref):
    #print('\ninputvectors: \n', inputvec, '\n\n', 'input: \n', input)
    #this function is a mess full of redundancies
    #this assumes that no 2 atoms will ever have the same x value in the same molecule
    #bonds is for the conditional final rotation - saves the 2 aromatic atoms bonded to current target atom
    #all 3 points make up the plane that dictates rotation
    bonds = carbons.iloc[0, index]
    globalindex = inputt.loc[inputt['x'] == bonds].index.tolist()[0]
    bonds = inputt.loc[inputt['x'] == bonds]
    bonds = bonds['Bonded To'].tolist()[0]
    bonds = inputt.iloc[[bonds[0]-1, bonds[1]-1], 0:3]
    bonds = bonds.transpose()
    


    #for printing a whole bunch of information
    testing = False



    angle = 120
    if testing:
        print('\nOriginal input: \n',inputt)
        print('\nInput Vectors C-H: \n', inputvec)
        print('\nBonds: \n', bonds)
        print('\nCarbons: \n', carbons)
        print('\nReference Coordinate: \n', refcoord)
        print('\nIndex: ', index)

        angle = getAngle(bonds.iloc[:, 0].to_list(), carbons.iloc[0:3, index].to_list(), bonds.iloc[:, 1].to_list())
        print('\nC-C-H angle pre-transformation(s) (should be multiple of 120): \n', angle)

    global globalcounter

    if refcoord[0] == -5.09099097:
        globalcounter+=1
    
    toobig = False
    if angle > 140 or angle < 110:
        print('Test #: ', globalcounter)
        print('\nAngle: ', angle)
        print('Reference Coordinate: \n', refcoord)
        toobig = True
    
    #vector to translate to 0, relative to current aromatic C being looked at
    originvector = 0 - carbons.iloc[0:3, index]
    input = input.iloc[:, 0:3].transpose()
    
    #first translation
    input.loc[len(input)] = 1
    bonds.loc[len(input)] = 1

    final = translation(originvector, input)
    bonds = translation(originvector, bonds)
    carbons = translation(originvector, carbons)

    if testing:
        print('\nOrigin Vector: \n', originvector)
        print('\nObject - Translated to origin: \n', final.transpose())
        print('\nBonds - Translated to origin: \n', bonds.transpose())
        print('\nCarbons - Translated to origin: \n', carbons.transpose())

    #### first rotation ####
    vec1 = inputvec.iloc[index]

    if testing:
        print('\nVec1 Ctarget-H: \n', vec1)

    vec1 = vec1 / np.linalg.norm(vec1)

    if testing:
        print('\nNormalized Vec1 Ctarget-H: \n', vec1)

    refvec = -refvec

    if testing:
        print('\nReference vector (first rotation): \n', refvec)

    refvec = refvec / np.linalg.norm(refvec)

    if testing:
        print('\nNormalized Reference Vector (first rotation): \n', refvec)

    theta = np.arccos(np.dot(vec1, refvec))
    rot_axis = pd.Series(np.cross(vec1, refvec), index=vector, dtype=np.float64)
    rot_axis = rot_axis / np.linalg.norm(rot_axis)

    if testing:
        print('\nAngle 1: \n', theta)
        print('\nRot Axis: \n', rot_axis)

    final = quat_rotation(theta, rot_axis, final)
    final.loc[len(input)] = 1
    
    bonds = quat_rotation(theta, rot_axis, bonds)
    bonds.loc[len(input)] = 1

    carbons = quat_rotation(theta, rot_axis, carbons)
    carbons.loc[len(carbons)] = 1
    carbons = carbons.astype(float)
    
    if testing:
        print('\nObject post-first-rotation: \n', final.transpose())
        print('\nBonds post-first-rotation: \n', bonds.transpose())
        print('\nCarbons post-first-rotation: \n', carbons.transpose())

    #### second rotation ####
    bonds = bonds.drop(4)
    testing1 = True
    if testing1:

        vec1 = bonds.iloc[:, 0] - carbons.iloc[0:3, index]
        vec1_2 = bonds.iloc[:, 1] - carbons.iloc[0:3, index]

        if testing:
            print('Second rotation vec1: \n', vec1)
            print('Second rotation vec1_2: \n', vec1_2)
        
        vec1 = vec1 / np.linalg.norm(vec1)
        vec1_2 = vec1_2 / np.linalg.norm(vec1_2)

        if testing:
            print('Second rotation vec1 (normalized): \n', vec1)
            print('Second rotation vec1_2 (normalized): \n', vec1_2)

        vec1 = np.cross(vec1.astype(float),vec1_2.astype(float))

        if testing:
            print('Vector for rotation: \n', vec1)

        vec1 = vec1 / np.linalg.norm(vec1)

        if testing:
            print('Vector for rotation (normalized): \n', vec1)

        #precalculated refvecs to 'sort of' vertically align undesirably rotated molecules
        if refcoord[0] == -5.09099097:
            refvec1 = pd.Series([-0.3863647478, -0.3318246062, 0.8605897468])
        else:
            refvec1 = pd.Series([0.1215815585,0.1627277551,0.9791514706], index=vector)

        if testing:
            print('Refvec for final rotation: \n', refvec1)

        theta = np.arccos(np.dot(vec1, refvec1))
        
        if testing:
            print('Final rotation angle: ',theta)

        rot_axis = pd.Series(np.cross(vec1, refvec1), index=vector, dtype=np.float64)
        
        if testing:
            print('Rotation axis: \n', rot_axis)

        rot_axis = rot_axis / np.linalg.norm(rot_axis)

        if testing:
            print('Rotation axis (normalized: \n', rot_axis)

        bonds.loc[len(bonds)] = 1
        final = quat_rotation(theta, rot_axis, final)
        bonds = quat_rotation(theta, rot_axis, bonds)
        carbons = quat_rotation(theta, rot_axis, carbons)
        final.loc[len(input)] = 1
        carbons.loc[len(carbons)] = 1

        if testing:
            print('All post-first-rotation objects: \n')
            print('1) Final: \n', final)
            print('2) Bonds: \n', bonds)
            print('3) Carbons: \n', carbons)

    #### final translation ####
    testing1 = True
    if testing1:
        bonds.loc[len(bonds)] = 1
        bonds = translation(refcoord,bonds)
        carbons = translation(refcoord,carbons)
        bonds = bonds.drop(3)
        carbons = carbons.drop(3)

    final = translation(refcoord,final)
    final = final.drop(3)
    
    if testing:
        print('\nFinal post-translation: \n',final.transpose())
        
        angle1 = getAngle(bonds.iloc[:, 0].to_list(), carbons.iloc[0:3, index].to_list(), refcoordpent.to_list())
        print('\nBond 0 section: \n', bonds.iloc[:, 0].to_list())
        print('\nCarbons section: \n', carbons.iloc[0:3, index].to_list())
        print('\nBond 1 section: \n', bonds.iloc[:, 1].to_list())
        
        print('\nC-C-H angle post-transformation(s) (should be multiple of 120): \n', angle1)
        print('\nAngle comparison post-pre: ', angle, 'vs', angle1, '\nDifference', angle1 - angle)
        print('\nTest #: ', globalcounter)
        print('\nIndex: ', index)

        if angle1 > 130 or angle1 < 110:
            toobig = True

    if toobig:
        #print('a')
        #pause()
        aaa = 1

    if testing:
        if globalcounter >= 0:
            #print('b')
            pause()

    return final.transpose()


def find_midpoints(origin1, origin2, target1, target2):
 
    origin_vector = (origin2 - origin1) / 2
    target_vector = (target2 - target1) / 2

    origin_point = origin1 + origin_vector
    target_point = target1 + target_vector

    return origin_point, target_point

def translation(tvec, shape):
    tx = tvec[0]
    ty = tvec[1]
    tz = tvec[2]
    
    trans_matrix = np.array([
        [1, 0, 0, tx],
        [0, 1, 0, ty],
        [0, 0, 1, tz],
        [0, 0, 0, 1]
    ])
    
    return np.matmul(trans_matrix, shape)

def find_vectors_angle(vector1, vector2):
    
    #obtaining angle between vectors
    unit_vector_1 = vector1 / np.linalg.norm(vector1)
    unit_vector_2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    theta = np.arccos(dot_product)
    
    #cross product
    crossp = np.cross(vector1, vector2)
    
    return theta, crossp

def rotation(theta, rotation_axis, shape):
    #implement to combine translation & rotation in one step? probably better
    rx = rotation_axis[0]
    ry = rotation_axis[1]
    rz = rotation_axis[2]
    
    rot_matrix = np.array([
        [np.cos(theta) + np.square(rx) * (1 - np.cos(theta)), #1,1
         rx * ry * (1 - np.cos(theta)) - rz * np.sin(theta), #1,2
         rx * rz * (1 - np.cos(theta)) + ry * np.sin(theta), #1,3
         0], #1,4
        
        [ry * rx * (1 - np.cos(theta)) + rz * np.sin(theta), #2,1
         np.cos(theta) + np.square(ry) * (1 - np.cos(theta)), #2,2
         ry * rz * (1 - np.cos(theta)) - rx * np.sin(theta), #2,3
         0], #2,4
        
        [rz * rx * (1 - np.cos(theta)) - ry * np.sin(theta), #3,1
         rz * ry * (1 - np.cos(theta)) + rx * np.sin(theta), #3,2
         np.cos(theta) + np.square(rz) * (1 - np.cos(theta)), #3,3
         0], #3,4
        
        [0, 0, 0, 1] #4,n
    ])
    
    #print('rotation matrix: \n', rot_matrix, '\n\nshape matrix (pre rotation): \n', shape, '\n')

    return np.matmul(rot_matrix, shape)

def quat_rotation(theta, rot_axis, input):


    Q = [rot_axis[0], rot_axis[1], rot_axis[2]]


    r = Rot.from_quat(
        [Q[0] * np.sin(theta/2), Q[1] * np.sin(theta/2), Q[2] * np.sin(theta/2), np.cos(theta/2)]
        )

    return np.matmul(r.as_matrix(), input[:-1])


def ref_shape_vectors(origin_list, target_list):

    xs = target_list[0] - origin_list[0]
    ys = target_list[1] - origin_list[1]
    zs = target_list[2] - origin_list[2]

    return [xs, ys, zs]


def fetch_files(dir_input=glob.glob('./gauss_files/input/*'), dir_ref='./gauss_files/reference/*'):
    smilist = []

    with open(dir_input[0], 'r') as file:
        for i in file:
            smilist.append(i)

    return smilist, glob.glob(dir_ref)

def open_file(filename):

    file = []
    raw = []

    with open(filename) as f:
        for i in f:
            file.append(i.split())
            raw.append(i)

    return file, raw

def grab_sects(file):
    '''Accepts entire .gjf file and returns coordinate and indicator sections'''

    atom_indicators = []
    coord_sect = []
    full = []

    for i in range(8, len(file)):

        if not file[i]:
            return (
                pd.DataFrame(coord_sect, columns = coords), 
                i, 
                pd.DataFrame(atom_indicators, columns = ['Atom']), 
                pd.DataFrame(full, columns = ['Atom', 'x', 'y', 'z'])
            )

        atom_indicators.append(file[i][0])
        coord_sect.append(file[i][1:4])
        full.append(file[i][0:4])
        
def ref_vectors(ref, atom_pair_1 = [128,170], atom_pair_2 = [46, 169]):
    
    coord_pair_1 = np.array([

        ref.loc[atom_pair_1[0] - 1].astype('float'),
        ref.loc[atom_pair_1[1] - 1].astype('float')

    ])

    
    coord_pair_2 = np.array([
    
        ref.loc[atom_pair_2[0] - 1].astype('float'),
        ref.loc[atom_pair_2[1] - 1].astype('float')
    
    ])

    
    ref_vectors = np.array([

        coord_pair_1[1] - coord_pair_1[0],
        coord_pair_2[1] - coord_pair_2[0]

    ])

    return ref_vectors, coord_pair_1, coord_pair_2

    #coord_mx = fix_shape(ref[0])
    #coord_mx = translation( [0,0,0], coord_mx)
    #print(np.full((1,42),1))
    
    #a = tr.translation( [0,0,0], fixed.astype('float') )
    #print(coord_mx)
    #print(ref[0].shape)

def input_vectors(input):
    veclist = []
    hydrogens = input.loc[input['atom'] == 'H']
    counter = 0
    for i in hydrogens['Bonded To']:
        
        if input.iloc[i[0]-1,4] == True:
            atomcoord = (input.loc[ [i[0]-1] , coords ]).to_numpy()
            hcoord = (hydrogens.iloc[counter][0:3]).to_numpy()
            veclist.append( hcoord - atomcoord[0] )
            #print(type(veclist)

        counter += 1

    return veclist

def all_combinations(list1, list2):

    return list(itertools.product(list1, list2))


def fix_bondsect(file, position=168):
    counter = 0
    replace = 1
    bondsect = []
    bondsectf = []

    for i in range(8, len(file)-2):
        #print(i)
        #print(file[i])
        if file[i] == '\n':
            counter+=1
            replace = replace + i
        #print(counter, replace)
        if counter == 1 and (i != replace - 1):
            bondsect.append(file[i+1][:-1].split())

    for i in bondsect:
        temp = []
        for s in i:
            if '.' in s:
                temp.append(float(s))
            else:
                temp.append(int(s) + position)
        bondsectf.append(temp)

    return bondsectf

def target_atoms():

    
    print('a')


def zmatrix(xyzfile):

    zmat = cc.Cartesian.read_xyz(xyzfile, start_index=1).get_zmat()

    return zmat

def xyzfile(molecule):

    with open('temp.xyz', 'w')  as f:

        f.write(str(len(molecule)) + '\n')
        f.write('tempfile\n')

        for line in molecule:
            line = ' '.join([str(item) for item in line])
            f.write(line + '\n')

def get_breakpoints(list):

    breakpoints = []
    counter = 0

    for i in list:
        if i != counter:
            counter = i
            breakpoints.append(i - 7)
            counter += 1
        else:
            counter += 1

    return breakpoints

def zipall(zipname):

    tmp_path = './gauss_files/output'
    save_path = './gauss_files/output'
    timestr = zipname + time.strftime("%Y%m%d-%H%M%S")
    
    zipname = os.path.join(save_path, timestr)
    
    shutil.make_archive(zipname, 'zip', tmp_path)

    return zipname

def write_file(filename, reference, final1, input1, lengths1, h1, c1, final2, input2, lengths2, h2, c2, timestr,
visual = 'default'):
    '''
    print('\ninput1: \n', input1)
    print('\nh1: \n', h1)
    print('\ninput2: \n', input2)
    print('\nh2: \n', h2)
    pause()
    '''
    save_path = './gauss_files/output/' + timestr
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    filename = os.path.join(save_path, filename)
    
    #totalatoms = 170 + input1.shape[0] + input2.shape[0]

    #change in connectivity on these specific lines is necessary, specline = 'special line'
    specline = [46, 50, 1.0, 169 + len(input1) + h2, 1.0]
    
    specline = [str(x) for x in specline]
    specline = ' '.join(specline)

    specline2 = [128, 132, 1.0, 169 + h1, 1.0]
    specline2 = [str(x) for x in specline2]
    specline2 = ' '.join(specline2)
    #######

    if visual == 'default':
        rangedef = range(0,176)
        rangecount = 169
    
    #visual = 'pentose'
    if visual == 'pentose':
        
        Ro = range(8)
        Ra = range(40 + 7, 47 + 7)
        Rb = range(49 + 7, 53 + 7)
        Rc = range(122 + 7, 129 + 7)
        Rd = range(131 + 7, 136 + 7)

        ranged = itertools.chain(Ro, Ra, Rb, Rc, Rd)
        rangedef = copy.deepcopy(ranged)
        rangedef = list(rangedef)
        rangecount = len(rangedef)

        


    with open(filename, 'w')  as f:

        #top section, specs for simulation & trimer coordinates
        for i in rangedef:
            f.write(' '.join(reference[i]) + '\n')
        
        #input1 coordinates
        for i in range(input1.shape[0]):
            f.write(
                str(input1['atom'][i]) + ' ' + 
                str(final1['x'][i]) + ' ' +
                str(final1['y'][i]) + ' ' +
                str(final1['z'][i]) + ' ' +
                '\n'
            )
        
        #input2 coordinates
        for i in range(input2.shape[0]+1):
            try:
                f.write(
                    str(input2['atom'][i]) + ' ' + 
                    str(final2['x'][i]) + ' ' +
                    str(final2['y'][i]) + ' ' +
                    str(final2['z'][i]) + ' ' +
                    '\n'
                )
            except:
                fff = 0
        
        if visual == 'default':
            #connectivity for trimer up to atom 169
            #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            if visual == 'default':
                rangedef = rangedef[8:]
            if visual == 'pentose':
                rangedef = rangedef[8:]
                breakpoints = get_breakpoints(rangedef[1:])
                print('here',rangedef, breakpoints)

            f.write('\n')

            sequencer = 0
            for i in rangedef:
                
                sequencer += 1


                #exception for trimer, atom 46
                if (i - 7)  == 46:
                    f.write(specline)
                    f.write('\n')
                
                #exception for trimer, atom 128
                elif (i - 7) == 128:
                    f.write(specline2)
                    f.write('\n')

                #everything else
                else:
                    f.write(' '.join(reference[i + 171]) + '\n')

                            
                            
            
            #this loop writes the input1 connectivity section
            #the try/except is necessary since one element was dropped from input1, yet the original indexing was retained
            for i in range(rangecount, rangecount + input1.shape[0] - len(input1.loc[input1['atom'] == 'H'])):
                f.write(
                    str(i) + ' '
                    )

                position1 = i - rangecount
                
                if (c1 == position1):
                    rangeb = range(len(input1.loc[i - rangecount,'Bonded To']) - 1)
                else:
                    rangeb = range(len(input1.loc[i - rangecount,'Bonded To']))

                try:
                    for b in rangeb:
                        if (input1.loc[i - rangecount,'Bonded To'][b] > lengths1):
                            #print('hello!')
                            displace = rangecount - 2
                        else:
                            displace = rangecount - 1
                        

                        f.write(
                            str(input1.loc[i - rangecount,'Bonded To'][b] + displace) + ' ' +
                            str(round(float(input1.loc[i - rangecount,'Bond Orders'][b]),1)) + ' '
                        )
                except:
                    fff = 1

                f.write('\n')

            #hydrogens don't need bonding information - this minimizes redundancy            
            for i in range(rangecount + input1.shape[0] - len(input1.loc[input1['atom'] == 'H']), rangecount + input1.shape[0]):
                f.write(
                    str(i) + '\n'
                )

            
            #this loop writes the input2 connectivity section
            #the try/except is necessary since one element was dropped from input2, yet the original indexing was retained
            for i in range(169 + input1.shape[0], 169 + input1.shape[0] + input2.shape[0] - len(input2.loc[input2['atom'] == 'H'])):

                f.write(
                    str(i) + ' '
                    )
                
                position2 = i - 169 - input1.shape[0]

                if (c2 == position2):
                    rangeb = range(len(input2.loc[i - 169 - input1.shape[0],'Bonded To']) - 1)
                else:
                    rangeb = range(len(input2.loc[i - 169 - input1.shape[0],'Bonded To']))

                try:

                    for b in rangeb:
                        
                        #print('\n')
                        if (input2.loc[i - 169 - input1.shape[0],'Bonded To'][b] > lengths2):
                            #print('hello!')
                            displace = 168 - 1
                        else:
                            displace = 168

                        #print('argument: ', input2.loc[i-169-input1.shape[0],'Bonded To'])
                        #print('write: ', input2.loc[i - 169 - input1.shape[0],'Bonded To'][b] + displace + len(input1))
                        #print(input2.loc[i - 169 - input1.shape[0],'Bonded To'][b])
                        #print('h2: ',h2)
                        #print('index2: ', lengths2)
                        
                        
                        f.write(
                            str(input2.loc[i - 169 - input1.shape[0],'Bonded To'][b] + displace + len(input1)) + ' ' +
                            str(round(float(input2.loc[i - 169 - input1.shape[0],'Bond Orders'][b]),1)) + ' '
                        )
                except:
                    #print('AAAAA')
                    fff = 2

                f.write('\n')


            #hydrogens don't need bonding information - this minimizes redundancy            
            for i in range(169 + input1.shape[0] + input2.shape[0] - len(input2.loc[input2['atom'] == 'H']), 169 + input1.shape[0] + input2.shape[0]):
                f.write(
                    str(i) + '\n'
                )
        '''
        '''

    #########
    '''
        f.write(str(specline))

        for i in range(179 + 129, len(reference) - 2):
            f.write(
                ' '.join(reference[i]) + '\n'
            )

        #here
        for i in range(1, len(molecule)):
            f.write( 
                ' ' + str(atomlist[i]) + '       ' + str(molecule[i][0]) + '   ' +
                str(molecule[i][1]) + '   ' + str(molecule[i][2]) + '\n'
            )

        for i in range(1, len(molecule2)):

            f.write( 
                ' ' + str(atomlist2[i]) + '       ' + str(molecule2[i][0]) + '   ' +
                str(molecule2[i][1]) + '   ' + str(molecule2[i][2]) + '\n'
            )

        for i in range(coordref, len(fullref) - 2):
            f.write(fullref[i])

        for i in bondlist:
            for s in i:
                f.write(' ' + str(s) + ' ')
            f.write('\n')

        for i in bondlist2:
            for s in i:
                f.write(' ' + str(s) + ' ')
            f.write('\n')
        '''

def get_aromatics(input):
    counter = 0
    hydrogens = input.loc[input['atom'] == 'H'].reset_index().drop(['index'], axis=1)

    #empty dataframe
    aromatics = pd.DataFrame()

    for i in range(len(hydrogens)):
        goodcarbon = True
        bond = hydrogens['Bonded To'].iloc[i - counter]
        aromatic = input['Aromaticity'].iloc[bond[0]-1]

        if not aromatic:
            hydrogens = hydrogens.drop(i)
            counter += 1
            goodcarbon = False
        
        if goodcarbon:
            aromatics = pd.concat([aromatics, input.iloc[bond[0]-1]], axis=1)
        
        
    return aromatics.transpose(), hydrogens


def getsmiles(file):

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "mol2")

    directory = '../gauss_files/input/SICS_ccpvtz_del-Copy.mol2'

    #with open(directory, 'r')  as f:
        #print(f.read())

    mol = openbabel.OBMol()
    obConversion.ReadString(mol, 'C1=C(C=c2c(=C1)cccc2)OC')

    smiles = obConversion.WriteString(mol)
    obConversion.WriteFile(mol, '1abc.gjf')

    return smiles

def graphtest(reference, input, final, filename):

    save_path = '../gauss_files/output'
    filename = os.path.join(save_path, filename)
    rangedef = range(0,8)
    
    with open(filename, 'w')  as f:

        #top section, specs for simulation & trimer coordinates
        for i in rangedef:
            f.write(' '.join(reference[i]) + '\n')

        if type(final) != str:
            print('AAA')
            final['atom'] = input['atom']

            for row in final.iterrows():
                f.write(row[1]['atom'] + ' ' +
                str(row[1][0]) + ' ' + str(row[1][1]) + ' ' + str(row[1][2]) 
                + '\n')
        else:
            for row in input.iterrows():
                f.write(row[1]['atom'] + ' ' +
                str(row[1]['x']) + ' ' + str(row[1]['y']) + ' ' + str(row[1]['z']) 
                + '\n')


def to_db_mysql(zip_directory, smiles_pair, TABLE_NAME):

    #convert tuple to str to avoid mysql conversion errors
    smiles_pair = str(smiles_pair)

    TABLES = {}
    TABLES[TABLE_NAME] = (

    "CREATE TABLE `" + TABLE_NAME + "` (" +
    "  `id` int(11) NOT NULL AUTO_INCREMENT," +
    "  `path_to_file` varchar(50) NOT NULL," +
    "  `test_date` date NOT NULL," +
    "  PRIMARY KEY (`id`)" + 
    ") ENGINE=InnoDB"
    
    )

    add_entry = ("INSERT INTO `" + TABLE_NAME + "` " +
               "(path_to_file, test_date) " +
               "VALUES (%s, %s)")

    data_entry = (zip_directory, date.today())

    global DB_NAME
    mysql_api.add_value(DB_NAME, TABLES, add_entry, data_entry)
    
def to_db_mongo(COLLECTION, batch, DB_NAME = 'dna'):
    mongo_api.insert(DB_NAME, COLLECTION, batch)


def gauss_gen(zipname='test_', name='test'):
    
    name = name + time.strftime("%Y%m%d_%H%M%S")
    #empty batch for the db
    batch = list()
    #empty .~/tmp before starting each test 
    dir = './gauss_files/tmp/*'
    files = glob.glob(dir)
    for f in files:
        os.remove(f)

    #all file directories, [0] > input files, [1] > reference file
    files = fetch_files()

    #entire reference file 
    ref, raw = open_file(files[1][0])

    #refsects - [0] > coordinate sect, [1] > end pos of coord sect, [2] > atomlist
    refsects = grab_sects(ref)

    #reference vectors dna molecule
    refvectors = ref_vectors(refsects[0])

    refvec1 = pd.Series(refvectors[0][0], index=coords)
    refvec2 = pd.Series(refvectors[0][1], index=coords)

    refcoord1 = pd.Series(refvectors[1][1], index=coords)
    refcoord2 = pd.Series(refvectors[2][1], index=coords)

    #all possible file combinations inside input folder
    combinations = itertools.combinations(files[0], 2)

    print('\n Converting SMILES to Gaussian output...\n')
    counter = 0
    print('\nCurrent pair: \n\n')
    for pair in combinations:
        print(' {} \n'.format(pair))
        #print(pair)
        #separating a split list from a raw source file is useful when rewriting to new files
        #input1, raw_input1 = open_file(pair[0])

        #saving smiles codes for several uses
        smiles1 = pair[0]
        smiles2 = pair[1]

        input1 = mf.getdf(smiles1, name)
        input2 = mf.getdf(smiles2, name)

        aromatics1, hydrogens1 = get_aromatics(input1)
        h1list = hydrogens1['Bonded To'].tolist()
        aromatics2, hydrogens2 = get_aromatics(input2)
        h2list = hydrogens2['Bonded To'].tolist()


        #NEVER FORGETTI WHEN I CASUALLY HAD DOUBLE BONDED HYDROGENS


        #columns we don't need for transformations
        todrop = ['index', 'atom', 'Aromaticity', 'Bonded To', 'Bond Orders']

        hydrogens1 = hydrogens1.reset_index().drop(todrop, axis=1).transpose()
        hydrogens2 = hydrogens2.reset_index().drop(todrop, axis=1).transpose()
        hydrogens1.loc[len(hydrogens1)] = 1
        hydrogens2.loc[len(hydrogens2)] = 1

        aromatics1 = aromatics1.reset_index().drop(todrop, axis=1).transpose()
        aromatics2 = aromatics2.reset_index().drop(todrop, axis=1).transpose()
        aromatics1.loc[len(aromatics1)] = 1
        aromatics2.loc[len(aromatics2)] = 1

        #hn is position of first hydrogen for molecule n
        
        input1vectors = pd.DataFrame(input_vectors(input1), columns=vector)
        input2vectors = pd.DataFrame(input_vectors(input2), columns=vector)

        #print('\n hydrogens 1: \n', hydrogens1)
        #print('\n vectors 1: \n', input1vectors)
        #print('\n hydrogens 2: \n', hydrogens2)
        #print('\n vectors 2: \n', input2vectors)
        
        counterh1 = 0
        counterh1_2 = 0
        input1t = input1
        input2t = input2
        hydrogens1t = hydrogens1
        hydrogens2t = hydrogens2
        #the innermost for loop is inefficient as fuck, this entire logic should probably be optimized later
        for index in range(len(aromatics1.columns)):
            #lengths = input1.shape[0] - aromatics1.transpose().shape[0] + index
            input1 = input1t
            hydrogens1 = hydrogens1t

        
            #this next part is dedicated to getting counterh1 to skip non-aromatic hydrogens
            #print('\nhydrogens1: \n', hydrogens1)
            #hydrogens1 = hydrogens1.drop([index], axis=1).transpose().reset_index().drop(['index'], axis=1).transpose()
            #print('\nhydrogens1: \n', hydrogens1)
            #print('\nthe first: \n', hydrogens1.iloc[0:3,index])
            #print('\ninput1: \n', input1)
            #print('\nthe second: \n', input1.iloc[len(input1) - len(input1.loc[input1['atom'] == 'H']) + index + counterh1, 0:3])
            #pause()

            aa = (hydrogens1.iloc[0:3,index] ==
            input1.iloc[len(input1) - len(input1.loc[input1['atom'] == 'H']) + index + counterh1, 0:3])

            if not aa.all():
                    counterh1_2 += 1

            while not aa.all():
                aa = (hydrogens1.iloc[0:3,index] ==
                input1.iloc[len(input1) - len(input1.loc[input1['atom'] == 'H']) 
                + index + counterh1, 0:3])
               
                if not aa.all():
                    counterh1 += 1

            ###
            lengths1 = input1.shape[0] - len(input1.loc[input1['atom'] == 'H']) + index + counterh1



            

            #print('\npre-change input1: \n', input1)
            #print('\nlengths1: ', lengths1)
            c1 = input1.iloc[lengths1]['Bonded To'][0]-1
            input1 = input1.drop([lengths1]).reset_index().drop(['index'], axis=1)
            #print('\npost-change input1: \n', input1)


            final1 = fix_shape(aromatics1, input1, index, input1vectors, refvec1, refcoord1, atom128, input1t, ref)
            final1.columns = coords

            h1 = h1list[index][0] - 1

            counterh2 = 0
            counterh2_2 = 0
            for index2 in range(len(aromatics2.columns)):
                #print('\n INDEX2: ', index2, '\n\n\n')
    
                input2 = input2t
                hydrogens2 = hydrogens2t

                #print('\n2: ', theta2, rot_axis2, final2)

                

                #this next part is dedicated to getting h2 to skip non-aromatic hydrogens
                #print('\npre-change hydrogens2: \n', hydrogens2)
                #if not hydrogens2.shape[1]-1 == index2:
                #    hydrogens2 = hydrogens2.drop([index2], axis=1).transpose().reset_index().drop(['index'], axis=1).transpose()
                #print('\npost hydrogens2: \n', hydrogens2)
                #print('\nthe first: \n', hydrogens2.iloc[0:3,index2])
                #print(input2)
                #print('\n1st: \n', len(input2) - len(input2.loc[input2['atom'] == 'H']))
                #print('\n others: ', index2, counterh2)
                #print('\nthe second: \n', input2.iloc[len(input2) - len(input2.loc[input2['atom'] == 'H']) + index2 + counterh2, 0:3])
                
                
                #pause()
                bb = (hydrogens2.iloc[0:3,index2] ==
                input2.iloc[len(input2) - len(input2.loc[input2['atom'] == 'H']) + index2 + counterh2, 0:3])

                if not bb.all():
                    counterh2_2 += 1

                while not bb.all():
                    bb = (hydrogens2.iloc[0:3,index2] ==
                    input2.iloc[len(input2) - len(input2.loc[input2['atom'] == 'H']) 
                    + index2 + counterh2, 0:3])

                    if not bb.all():
                        counterh2 += 1
                

                lengths2 = input2.shape[0] - len(input2.loc[input2['atom'] == 'H']) + index2 + counterh2
                #print(len(input2.loc[input2['atom'] == 'H']))
                #print(input2.shape)
                #print('index2: ', index2, counterh2)
                #print('\npre-change input2: \n', input2)
                #print('\nlengths2: ', lengths2)
                c2 = input2.iloc[lengths2]['Bonded To'][0]-1


                input2 = input2.drop([lengths2])

                
                #print('\npost-change input2: \n', input2)

                '''
                print('all things  being sent to fix_shape: \n')
                print('aromatics2: \n', aromatics2)
                print('input2: \n', input2)
                print('index2: \n', index2)
                print('input2vectors: \n',input2vectors)
                print('refvec2: \n',refvec2)
                print('refcoord2: \n', refcoord2)
                print('atom46: \n', atom46)
                print('input2t: \n', input2t)
                print('ref: \n', ref)
                '''


                final2 = fix_shape(aromatics2, input2, index2, input2vectors, refvec2, refcoord2, atom46, input2t, ref)
                final2.columns = coords
                
                ###

                h2 = h2list[index2][0] - 1

                filename = 'test' +  str(counter + 1) + '_' + str(index + 1) + str(index2 + 1) + '.gjf'
                
                batch.append(
                    {
                        '_id' : counter,
                        'pair' : pair,
                        'directory' : './gauss_files/output/' + filename
                    }
                )
                counter += 1
                write_file(filename, ref, final1, input1, lengths1, h1, c1, final2, input2, lengths2, h2, c2, name)

    #zipname = zipall(zipname)
    to_db_mongo(name, batch)

#gauss_gen(zipname='extra rotation v3 ')

'''

#benzene coordinates taken from:
#https://github.com/choderalab/openmoltools/blob/master/openmoltools/chemicals/benzene/benzene.mol2
benzene =  np.array([
    [-1.2130, -0.6880, 0.0000],
    [-1.2030, 0.7060, 0.0000],
    [-0.0100, -1.3950, 0.0000],
    [0.0100, 1.3950, -0.0000],
    [ 1.2030, -0.7060, 0.0000],
    [1.2130, 0.6880, 0.0000],
    [ -2.1580, -1.2240, 0.0000],
    [ -2.1390, 1.2560, 0.0000],
    [-0.0180, -2.4810, -0.0000],
    [0.0180, 2.4810, 0.0000],
    [2.1390, -1.2560, 0.0000],
    [2.1580, 1.2240, 0.0000]
])

#2 opposing carbons to be used as reference for transformation
origin_point_1 = np.array([ benzene[0,0],benzene[0,1], benzene[0,2] ]).reshape(3,1)
origin_point_2 = np.array([ benzene[5,0],benzene[5,1], benzene[5,2] ]).reshape(3,1)
dist_origin = np.linalg.norm(origin_point_1 - origin_point_2)

#calculating sides of cube so that cross distance between opposite vertices = dist_origin
cross = np.sqrt((np.square(dist_origin)/3))

#cube modified to perfectly fit benzene reference across opposing vertices
cube = np.array([
    [-2,-2,-2], [-2, -2 + cross, -2], [-2,-2, -2 + cross], [-2,-2 + cross, -2 + cross],
    [-2 + cross,-2, -2 + cross], [-2 + cross,-2,-2], [-2 + cross,-2 + cross,-2],
    [-2 + cross,-2 + cross,-2 + cross]
])

#target vertices on cube to be alligned
target_1 = np.array([cube[0,0], cube[0,1], cube[0,2]]).reshape(3,1)
target_2 = np.array([cube[7,0], cube[7,1], cube[7,2]]).reshape(3,1)

#separating all x, y, z coordinates for both shapes
cubex, cubey, cubez = cube[:, 0], cube[:, 1], cube[:, 2]
benx, beny, benz = benzene[:, 0], benzene[:, 1], benzene[:, 2]

#redefining benzene coordinates to ignore hydrogens
benx, beny, benz = benzene[:-6, 0], benzene[:-6, 1], benzene[:-6, 2]

#origin and target coordinates
org, tar = find_midpoints(origin_point_1, origin_point_2, target_1, target_2)

#transformation vector calculation (this should be inside of a function on the final version)
transf_vector_shapes = tar - org

#transform to origin pre-rotation (move to function later)
#doesn't do anything in this particular case lolol midpoint is already 0,0,0
transf_vector_origin = 0 - org
mod_mx = fix_shape(benzene[:-6].transpose())
mod_mx = translation(transf_vector_origin, mod_mx)

#finding rotation angle and rotation axis, find more general solution to angle correction
angle, rot_axis = find_vectors_angle(origin_point_1, origin_point_2, target_1, target_2)
angle = 2*np.pi - angle

#applying rotation
#print('angle: ', angle, end='\n\n')
#print('rot axis: ', rot_axis, end='\n\n')
rot_axis_norm = rot_axis / np.linalg.norm(rot_axis)
mod_mx_1 = rotation(angle, rot_axis_norm, mod_mx)

#final translation
mod_mx_2 = translation(tar, mod_mx_1)
#print('final shape matrix: \n', mod_mx_2)
'''



'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#cube
ax.scatter(cubex, cubey, cubez)

#benzene
#pre
ax.scatter(benx, beny, benz)
#post
ax.scatter(mod_mx_2[0], mod_mx_2[1], mod_mx_2[2])

plt.savefig('f')
plt.show()
'''