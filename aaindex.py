################################################################################
#################                    AAIndex                   #################
################################################################################

"""
The AAindex is a database of numerical indices representing various physicochemical
and biochemical properties of amino acids and pairs of amino acids. The focus
on this module is on the AAIndex 1 database which stores the amino acid index of
20 numerical values for the 20 amino acids: http://www.genome.jp/aaindex/

Data format of AAI1:
************************************************************************
*                                                                      *
* Each entry has the following format.                                 *
*                                                                      *
* H Accession number                                                   *
* D Data description                                                   *
* R PMID                                                               *
* A Author(s)                                                          *
* T Title of the article                                               *
* J Journal reference                                                  *
* * Comment or missing                                                 *
* C Accession numbers of similar entries with the correlation          *
*   coefficients of 0.8 (-0.8) or more (less).                         *
*   Notice: The correlation coefficient is calculated with zeros       *
*   filled for missing values.                                         *
* I Amino acid index data in the following order                       *
*   Ala    Arg    Asn    Asp    Cys    Gln    Glu    Gly    His    Ile *
*   Leu    Lys    Met    Phe    Pro    Ser    Thr    Trp    Tyr    Val *
* //                                                                   *
************************************************************************

e.g ANDN920101
--------------------------------------------------------------------------------
H ANDN920101
D alpha-CH chemical shifts (Andersen et al., 1992)
R PMID:1575719
A Andersen, N.H., Cao, B. and Chen, C.
T Peptide/protein structure analysis using the chemical shift index method:
    upfield alpha-CH values reveal dynamic helices and aL sites
J Biochem. and Biophys. Res. Comm. 184, 1008-1014 (1992)
C BUNA790102    0.949
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    4.35    4.38    4.75    4.76    4.65    4.37    4.29    3.97    4.63    3.95
    4.17    4.36    4.52    4.66    4.44    4.50    4.35    4.70    4.60    3.95
--------------------------------------------------------------------------------

References 
----------
[1]: Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., and Kanehisa, 
     M.; AAindex: amino acid index database, progress report 2008. Nucleic Acids Res. 36, 
     D202-D205 (2008). [PMID:17998252]
[2]: https://github.com/harmslab/hops/blob/master/hops/features/data/util/aaindex2json.py
"""

#importing required modules and dependencies
import json
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
import sys, copy, re
import requests
import shutil
import csv
import urllib.request as request
from contextlib import closing
from collections import defaultdict
import itertools
from difflib import get_close_matches

aaindex_json = {}
aaindex_ftp_url = "ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1"
aaindex_http_url = "https://www.genome.jp/ftp/db/community/aaindex/aaindex1"
aaindex_filename = "aaindex1"
DATA_DIR = 'data'

# def parse_categories(self, aaindex_category_file='aaindex_to_category.txt'):
#     """
#     Parse category file which maps each AAI record in the databse into 1 of 8 categories.
#     Category file and parsing code inspired from: https://github.com/harmslab/hops. 

#     Parameters
#     ----------
#     : aaindex_category_file : str
#         Name of category file to parse (default is "aaindex_to_category.txt"). If parameter
#         is a file name then it will be imported from the 'data' directory by default, else
#         it will be imported from its full path.

#     Returns
#     -------
#     : aaindex_category : dict
#         Dictionary that maps each AAI record into 1 of 8 categories.
#     """
#     #if input parameter is a full path, read it else read from default DATA_DIR
#     if os.path.isfile(aaindex_category_file):
#         f = open(aaindex_category_file,'r')
#     else:
#         try:
#             f = open((os.path.join('aaindex', DATA_DIR, aaindex_category_file)),'r')
#         except IOError:
#             print('Error opening AAIndex1 category file, check {} in {} directory'
#                 .format(aaindex_category_file, DATA_DIR))

#     #get total number of lines in file
#     total_lines = len(f.readlines(  ))

#     #open new file in data directory to store parsed category file
#     f.seek(0)
#     category_output_file = "aaindex_categories.txt"
#     f_out = open((os.path.join('aaindex', DATA_DIR, category_output_file)), "w")

#     #lines starting with '#' are file metadata so don't write these to parsed output file
#     for line in f.readlines():
#         if not (line.startswith('#')):
#             f_out.write(line)
#     #close both files
#     f.close()
#     f_out.close()

#     #open parsed category file for reading
#     with open(os.path.join('aaindex', DATA_DIR, category_output_file)) as f_out:
#         reader = csv.reader(f_out, delimiter="\t")
#         d = list(reader)

#     aaindex_category = {}

#     #iterate through all lines in parsed file and store AAI indices as keys
#     #  in the output dict and their respective categories as the values of dict
#     for i in range(0,len(d)):
#         for j in range(0,len(d[i])):
#             category_substring = d[i][1].strip()
#             category_substring = category_substring.split(" ",1)
#             aaindex_category[d[i][0]] = category_substring[0]

#     categories = aaindex_category

def get_all_categories(aaindex_categories_file="aaindex_categories.txt"):
    """
    Return list of all indices, their description and associated categories from 
    parsed aaindex_categories.txt file created from parse_categories() function.

    Parameters
    ----------
    :aaindex_categories_file : str (default = aaindex_categories.txt)
        name of parsed categories file, created from parse_categories() function. 
    
    Returns
    -------
    aaindex_category : dict
        Dictionary that maps each AAI record into 1 of 8 categories.
    """
    if not (os.path.isfile(os.path.join('aaindex',DATA_DIR,aaindex_categories_file))):
        # parse_categories()
        print('here')

    try:
        with open(os.path.join('aaindex', DATA_DIR,'aaindex_categories.txt')) as f_out:
            reader = csv.reader(f_out, delimiter="\t")
            d = list(reader)
    except IOError:
        print('Error opening AAIndex1 category file.')

    aaindex_category = {}

    #iterate through all lines in parsed file and store AAI indices as keys
    #   in the output dict and their respective categories as the values of dict
    for i in range(0,len(d)):
        for j in range(0,len(d[i])):
            category_substring = d[i][1].strip()
            category_substring = category_substring.split(" ",1)
            aaindex_category[d[i][0]] = category_substring[0]

    return aaindex_category

def parse_aaindex_to_json():
    """
    Parse AAI database into JSON format. Each AAI record will be indexed by
    its feature code/index code, and will be in the format as shown in the
    example above. The file will be stored in a json file called according
    to (aaindex_filename+'.json') variable.

    Parameters
    ----------
    None

    Returns
    -------
    aaindex_json : dict
        parsed AAI database in dict form.
    """
    #initialise keys of AAI database
    template_dict = {"H":[],
                        "D":[],
                        "R":[],
                        "A":[],
                        "*":[],
                        "T":[],
                        "J":[],
                        "C":[],
                        "I":[]}

    #open AAI file for reading and parsing, by default it should be stored in DATA_DIR
    try:
        tmp_filepath = os.path.join('aaindex',DATA_DIR,aaindex_filename)
        f = open(tmp_filepath,'r')
    except IOError:
        print('Error opening file, check filename = {} and is stored in \
                {} directory.'.format(aaindex_filename,DATA_DIR))

    #read lines of file
    lines = f.readlines()
    f.close()

    clean_up_pattern = re.compile("\"")

    #initilaise parsed AAI database dictionary
    aaindex_json = {}
    current_dict = copy.deepcopy(template_dict)

    '''
    iterate through each line in the AAI database, parsing each record and its
    amino acid values into its own entry in the dictionary. Each index/record
    is seperated by a '//'. Remove any duplicate records, set any missing
    ('-') or NA amino acid values to 0. Store resulting dict into aaindex_json
    instance variable.
    '''
    for l in lines:

        if l.startswith("//"):

            # deal with meta data
            name = " ".join(current_dict["H"])
            name = clean_up_pattern.sub("'",name)

            description = " ".join(current_dict["D"])
            description = clean_up_pattern.sub("'",description)

            citation = "{} '{}' {}".format(" ".join(current_dict["A"]),
                                            " ".join(current_dict["T"]),
                                            " ".join(current_dict["J"]))

            citation = citation + "; Kawashima, S. and Kanehisa, M. \
                'AAindex: amino acid index database.'  Nucleic Acids Res. 28, 374 (2000)."

            citation = clean_up_pattern.sub("'",citation)

            notes = " ".join(current_dict["*"])
            notes = clean_up_pattern.sub("'",notes)

            # parse amino acid data
            aa_lines = current_dict["I"]

            aa_names = aa_lines[0].split()
            row_0_names = [aa.split("/")[0] for aa in aa_names]
            row_1_names = [aa.split("/")[1] for aa in aa_names]

            row_0_values = aa_lines[1].split()
            row_1_values = aa_lines[2].split()

            values = {}
            for i in range(len(row_0_values)):
                try:
                    values[row_0_names[i]] = float(row_0_values[i])
                except ValueError:
                    values[row_0_names[i]] = "NA"

                try:
                    values[row_1_names[i]] = float(row_1_values[i])
                except ValueError:
                    values[row_1_names[i]] = "NA"

            # look for duplicate name entries
            try:
                aaindex_json[name]
                err = "duplicate value name ({})".format(name)
                raise ValueError('Duplicate AAI Record name found')
            except KeyError:
                pass

            aaindex_json[name] = {"description":description,
                                "refs":citation,
                                "notes":notes,
                                "values":values}

            current_dict = copy.deepcopy(template_dict)
            continue

        this_entry = l[0]
        if l[0] != " ":
            current_entry = this_entry

        current_dict[current_entry].append(l[1:].strip())

    categories = get_all_categories()

    #append '-' to each aa index entry to account for missing AA in a protein sequence
    for index in aaindex_json:

        #add category to aaindex file
        aaindex_json[index]['category'] = categories[index]
        aaindex_json[index]['values']['-'] = 0

        #set any NA amino acid values to 0
        for val in aaindex_json[index]['values']:
            if aaindex_json[index]['values'][val] == 'NA':
                aaindex_json[index]['values'][val] = 0

    # for index in
    #save parsed dictionary into JSON format to DATA_DIR
    with open((os.path.join('aaindex', DATA_DIR, aaindex_filename+'.json')),'w') as output_F:
        json.dump(aaindex_json, output_F, indent=4, sort_keys=True)

    return aaindex_json

#if parsed json of AAIndex already in file then read it and return  <- Not working at the moment
if (os.path.isfile(os.path.join('aaindex', DATA_DIR, aaindex_filename+'.json'))):
    with open(os.path.join('aaindex', DATA_DIR, aaindex_filename+'.json')) as aai_json:
        aaindex_json = json.load(aai_json)
else:
    #parse AAIndex database file into JSON format
    aaindex_json = parse_aaindex_to_json()

def parse_aaindex_to_category(aaindex_category_file='aaindex_to_category.txt'):
    """
    Parse category file which maps each AAI record in the databse into 1 of 8 categories.
    Category file and parsing code inspired from:
    https://github.com/harmslab/hops
    Parameters
    ----------
    :aaindex_category_file : str
        Name of category file to parse (default is "aaindex-to-category.txt").
    Returns
    -------
    :aaindex_category : dict
        Dictionary that maps each AAI record into 1 of 8 categories.
    """
    #if input parameter is a full path, read it else read from default DATA_DIR
    if os.path.isfile(aaindex_category_file):
        f = open(aaindex_category_file,'r')
    else:
        try:
            f = open((os.path.join('aaindex',DATA_DIR, aaindex_category_file)),'r')
        except IOError:
            print('Error opening AAIndex1 category file, check {} in {} directory'
                .format(aaindex_category_file, DATA_DIR))

    #get total number of lines in file
    total_lines = len(f.readlines(  ))

    #open new file in data directory to store parsed category file
    f.seek(0)
    category_output_file = "aaindex-to-category-parsed.txt"
    f_out = open((os.path.join('aaindex',DATA_DIR, category_output_file)), "w")

    #lines starting with '#' are file metadata so don't write these to parsed output file
    for line in f.readlines():
        if not (line.startswith('#')):
            f_out.write(line)
    f.close()
    f_out.close()

    #open parsed category file for reading
    with open(os.path.join('aaindex',DATA_DIR, category_output_file)) as f_out:
        reader = csv.reader(f_out, delimiter="\t")
        d = list(reader)
    f_out.close()

    aaindex_category = {}

    #iterate through all lines in parsed file and store AAI indices as keys
    #   in the output dict and their respective categories as the values of dict
    for i in range(0,len(d)):
        for j in range(0,len(d[i])):
            category_substring = d[i][1].strip()
            category_substring = category_substring.split(" ",1)
            aaindex_category[d[i][0]] = category_substring[0]

    return aaindex_category

def download_aaindex(download_using="ftp"):
    """
    If AAI database not found in DATA directory, then it will be downloaded from
    the dedicated FTP or HTTPS server from https://www.genome.jp/aaindex/.
    FTP is the default method used for downloading the database with the URL:
    "ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1". If you want to
    download by https, set the 'download_using' instance variable to 'https'.

    Parameters
    ----------
    None
    """
    url = ""

    if download_using=='ftp':
        url = aaindex_ftp_url
    elif download_using=='http' or download_using=='https':
        url = aaindex_http_url
    else:
        url = aaindex_ftp_url

    #fetch AAI database from URL if not present in DATA_DIR
    try:
        if not(os.path.isfile(os.path.join('aaindex',DATA_DIR, aaindex_filename))):
            try:
                with closing(request.urlopen(url)) as r:
                    with open((os.path.join('aaindex',DATA_DIR,aaindex_filename)), 'wb') as f:
                        shutil.copyfileobj(r, f)
                print('AAIndex1 successfully downloaded.')
            except requests.exceptions.RequestException:
                print('Error downloading or exporting AAIndex from the url {}.'.format(url))
        else:
            pass    #AAIndex already in folder
    except:
        raise OSError('Save Directory does not exist: {}.'.format('data'))
        
def get_amino_acids():
    """
    Get all canonical amino acid letters. The '-' value will also be included
    in the list from this function as it accounts for the abcense of any
    amino acid or gaps in a AAI record.

    Returns
    -------
    amino_acids : list
        List of all 20 canoniocal amino acid letters as found in each record
        of the AAI database.
    """
    amino_acids = list(aaindex_json[list(aaindex_json.keys())[0]]["values"].keys())
    amino_acids.sort()

    return amino_acids

def get_amino_acids_encoding():
    """
    Get one-hot encoding of amino acids.

    Returns
    -------
    onehot_encoded : np.ndarray
        one hot encoded array of the 20 canonical amino acids.
    """
    all_amino_acids = get_amino_acids()
    values = np.array(all_amino_acids[1:])    #convert amino acids to np array

    #encode amino acids with value between 0 and n_classes-1.
    label_encoder = LabelEncoder()
    integer_encoded = label_encoder.fit_transform(values)

    #encode amino acids as a one-hot numeric array.
    onehot_encoder = OneHotEncoder(sparse=False)
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.fit_transform(integer_encoded)

    return onehot_encoded
    
def get_record_codes():
    """
    Get list of all AAI index codes/Accession numbers for each record in the database.

    Returns
    -------
    records : list
        list of record names/index codes for all records in the AAI database.
    """
    records = list(aaindex_json.keys())
    records.sort()       #sort into alphabetical order

    return records

def get_num_records():
    """
    Calculate total number of records/indices in the AAI database.

    Returns
    -------
    len(get_record_codes()) : int
        number of indices/records found in the AAI database.
    """
    return len(get_record_codes())

def get_record_names():
    """
    Return a list of all index descriptions for all records in the AAI database.

    Returns
    -------
    desc : list
        list of descriptions for all records in the AAI database.
    """
    desc = []

    #iterate through database, appending all descriptions to desc list
    for name in list(aaindex_json.values()):
        desc.append(name['description'])

    return desc

def get_record_from_name(name):
    """
    Return full AAI database record details from its description. Search
    through the descriptions/names of all records in the database, returning
    record that matches name input parameter.

    Parameters
    ----------
    :name : str
        AAI database record feature name/description.

    Returns
    -------
    :record : dict/None
        dict of full AAI database record or None if not found.
    """
    correct_index = 0
    for index, value in aaindex_json.items():
        if value['description'] != name:
            continue
        else:
            correct_index = index
            return aaindex_json[correct_index]

    return None

def get_values_from_record(index_code):
    """
    Return amino acid values from database record of index in the AAI from
    its feature/index code.

    Parameters
    ----------
    index_code : str
        AAI database record feature code/index code.

    Returns
    -------
    values : dict
    amino acid values for specified AAI database record.
    """
    #stripping input of whitespace
    try:
        index_code = index_code.strip()
    except:
        raise TypeError('Input parameter {} is not of correct datatype string, got {}' \
            .format(index_code, type(index_code)))

    #check that inputted index_code/Accession num does exist in the AAI database
    if index_code not in (get_record_codes()):
        raise ValueError('Record index {} not found in AAI'.format(index_code))

    #get amino acid values for specified index
    values = (aaindex_json[index_code]['values'])

    return values

def get_ref_from_record(index_code):
    """
    Return reference details of a AAI database record from its feaure/index code.

    Parameters
    ----------
    index_code : str
    AAI database record feature code/index code.

    Returns
    -------
    refs : list
    reference details for specifed AAI database record.
    """
    #stripping input of whitespace
    try:
        index_code = index_code.strip()
    except:
        raise TypeError('Input parameter {} is not of correct datatype string, got {}' \
            .format(index_code, type(index_code)))

    #check that inputted index_code does exist in the AAI database
    if index_code not in (get_record_codes()):
        raise ValueError('Record index {} not found in AAI'.format(index_code))

    refs = (aaindex_json[index_code]['refs'])    #get refs from record

    return refs

def get_category_from_record(index_code):
    """
    Return category of a AAI database record from its feaure/index code.

    Parameters
    ----------
    index_code : str
        AAI database record feature code/index code.

    Returns
    -------
    cat : str
    category of AAI database record according to parsed aaindex_to_category.txt file.
    """
    #stripping input of whitespace
    try:
        index_code = index_code.strip()
    except:
        raise TypeError('Input parameter {} is not of correct datatype string, got {}' \
            .format(index_code, type(index_code)))

    #check that inputted index_code does exist in the AAI database
    if index_code not in (get_record_codes()):
        raise ValueError('Record index {} not in AAI database'.format(index_code))

    #get category from dict
    cat = categories[index_code]

    return cat

def __str__(self):
    return "AAIndex1 Database, with {} records - stored in {} ".format(
    self.get_num_records(),self.aaindex_filename)

def __repr__(self):
    return (json.dumps(self.aaindex_json, indent=4, sort_keys=True))

def __sizeof__(self):
    """ Return size of AAI database file """
    return os.path.getsize(os.path.isfile(os.path.join('aaindex',DATA_DIR,self.aaindex_filename)))

def __getitem__(self, record_code):
    """
    Return full AAI database record details from its feature/index code/Accession number.
    Parameters
    ----------
    : record_code : str
        AAI database record feature index/Accession number.
    Returns
    -------
    : record : dict
        dict of AAI database record.
    Usage
    -----
        aaindex = AAIndex()
        full_record = aaindex['ZIMJ680105']
    """
    #stripping input of whitespace
    try:
        record_code = record_code.strip()
    except:
        raise TypeError('Input parameter {} is not of correct datatype string, got {}' \
            .format(record_code, type(record_code)))

    #check that inputted record_code/Accession num does exist in the AAI database
    if record_code not in (self.record_codes()):
        raise ValueError('Record Index ({}) not found in AAIndex'.format(record_code))
    
    #get full record from parsed JSON
    record = self.aaindex_json[record_code]

    return 

categories = parse_aaindex_to_category()
