from __future__ import print_function

import mysql.connector
from mysql.connector import errorcode
from datetime import date, datetime
import os

def pause():
    input('@@@ Press ENTER to continue')

'''
config = {
  'user': 'root',
  'password': 'dna',
  'host': 'mysql',
  'port': '3306',
  'database': 'test',
  'raise_on_warnings': True
}
'''

config = {
    'user':'root',
    'password':'dna',
    'host':'mysql',
    'port':'3306'
}

def getcursor():
    global config

    try:
        cnx = mysql.connector.connect(**config)
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)

    return cnx.cursor(), cnx

def closecursor(cursor, cnx):
    cursor.close()
    cnx.close()


def create_database(DB_NAME, cursor):
    
    #visual confirmation of database creation
    try:
        cursor.execute(
            "CREATE DATABASE {} DEFAULT CHARACTER SET 'utf8'".format(DB_NAME))
    except mysql.connector.Error as err:
        print("Failed creating database: {}".format(err))
        exit(1)

def use_db(DB_NAME, cursor, cnx):

    try:
        cursor.execute("USE {}".format(DB_NAME))
    except mysql.connector.Error as err:
        print("Database {} does not exist.".format(DB_NAME))
        if err.errno == errorcode.ER_BAD_DB_ERROR:
            create_database(DB_NAME, cursor)
            print("Database {} created successfully.".format(DB_NAME))
            cnx.database = DB_NAME
        else:
            print(err)
            exit(1)

def new_table(DB_NAME, TABLES, returncursor = False):
    
    cursor, cnx = getcursor()
    use_db(DB_NAME, cursor, cnx)

    #visual confirmation of writing tables
    for table_name in TABLES:
        table_description = TABLES[table_name]
        try:
            cursor.execute(table_description)            
        except mysql.connector.Error as err:
            if not err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
                print(err.msg)
    if returncursor:
        return cursor, cnx
    closecursor(cursor, cnx)

def add_value(DB_NAME, TABLES, command, data):

    cursor, cnx = new_table(DB_NAME, TABLES, returncursor = True)
    cursor.execute(command, data)
    cnx.commit()
    closecursor(cursor, cnx)

def rm_database(DB_NAME):

    cursor,cnx = getcursor()
    command = ('DROP DATABASE ' + DB_NAME)
    cursor.execute(command)
    closecursor(cursor,cnx)


def rm_table(DB_NAME, TABLE_NAME):

    cursor,cnx = getcursor()
    use_db(DB_NAME, cursor, cnx)
    command = ('DROP TABLE ' + TABLE_NAME)
    cursor.execute(command)
    closecursor(cursor,cnx)

############# NEEDS FIXING DOESN'T WORK
def rm_value(DB_NAME, TABLE_NAME, id = None):

    cursor,cnx = getcursor()
    use_db(DB_NAME, cursor, cnx)
    command = ('DELETE FROM ' + TABLE_NAME)
    
    if id != None:
        command = command + 'WHERE id = ' + id

    cursor.execute(command)
    closecursor(cursor,cnx)

def grab_smiles(DB_NAME, TABLE_NAME = 'smiles_' + str(date.today()).replace('-', '')):

    #first we fetch all files in the download directory
    allfile = set()
    dir = './tranches/'

    for file in os.listdir(dir):
        if file.endswith(".smi"):
            allfile.add(os.path.join(dir, file))
    
    #now we send it to the db
    TABLES = {}
    TABLES[TABLE_NAME] = ("CREATE TABLE `" + TABLE_NAME + "` (" +
    "  `id` int(11) NOT NULL AUTO_INCREMENT," +
    "  `smiles_chain` varchar(100) NOT NULL," +
    "  `zinc_id` varchar(20) NOT NULL," +
    "  `tranche` varchar(25) NOT NULL," +
    "  PRIMARY KEY (`id`)" + 
    ") ENGINE=InnoDB")

    add_entry = ("INSERT INTO `" + TABLE_NAME + "` " +
               "(smiles_chain, zinc_id, tranche) " +
               "VALUES (%s, %s, %s)")

    #each file is a tranche, so we gotta iterate through them
    for i in allfile:
        tranche = i[-10:]
        with open(i, 'r') as file:
            for e in file:
                if not e.startswith('smiles'):
                    smiles_chain = e.split()[0]
                    zinc_id = e.split()[1]
                    data_entry = (smiles_chain, zinc_id, tranche)
                    add_value(DB_NAME, TABLES, add_entry, data_entry)
                    
