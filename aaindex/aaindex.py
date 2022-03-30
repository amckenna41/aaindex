################################################################################
#################                    AAIndex                   #################
################################################################################

#importing required modules and dependencies
import json
import numpy as np
import os
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
import sys, copy, re
import requests
import shutil
from difflib import get_close_matches
import csv
import urllib.request as request
from contextlib import closing
from difflib import get_close_matches

#global vars
DATA_DIR = 'data'
AAINDEX_FILENAME = 'aaindex1'
DOWNLOAD_USING = 'ftp'

class AAIndex():
    """
            Python parser for AAindex1: Amino Acid Index Database
                      (**abbreviated to AAI onwards**)
    The AAindex is a database of numerical indices representing various physiochemical
    and biochemical properties of amino acids and pairs of amino acids. This class 
    stores the amino acid index of 20 numerical values for the 20 amino acids 
    (http://www.genome.jp/aaindex/).
    Some of the parsing functionality is patrially inspired off the AAIndex1 parser
    by harmslab: https://github.com/harmslab/hops/blob/master/hops/features/data/util/aaindex2json.py

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

    Methods
    -------
    parse_aaindex():
        parse AAI database into JSON format.
    parse_categories(aaindex_category_file):
        parse AAI indices into their respective categories.
    download_aaindex():
        download AAI database from its FTP server.
    get_amino_acids():
        get list of 20 amino acid letters.
    get_amino_acids_encoding():
        get one-hot encoding of amino acids.
    record_codes():
        get list of AAI record codes/Accession numbers.
    get_num_records():
        return total number of records in AAI database.
    get_record_names():
        return list of all descriptions for all records in AAI database.
    get_record_from_desc():
        return AAI record from its description.
    get_ref_from_record():
        return references for AAI record.
    get_category_from_record():
        return category of AAI index record.
    __getitem__():
        access full aaindex record using its record code.

    References
    ----------
    [1]: Nakai, K., Kidera, A., and Kanehisa, M.;  Cluster analysis of
        amino acid indices for prediction of protein structure and
        function.  Protein Eng. 2, 93-100 (1988)
    [2]: Tomii, K. and Kanehisa, M.;  Analysis of amino acid indices and
        mutation matrices for sequence comparison and structure
        prediction of proteins.  Protein Eng. 9, 27-36 (1996).
    [3]: Kawashima, S., Ogata, H., and Kanehisa, M.;  AAindex: amino acid
        index database.  Nucleic Acids Res. 27, 368-369 (1999).
    [4]: Kawashima, S. and Kanehisa, M.;  AAindex: amino acid index
        database.  Nucleic Acids Res. 28, 374 (2000).
    """
    def __init__(self):

        self.aaindex_module_path = os.path.dirname(os.path.abspath(sys.modules[self.__module__].__file__))

        #download AAI database using ftp or https
        if DOWNLOAD_USING =='ftp':
            self.url = "ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1"
        elif DOWNLOAD_USING =='http' or DOWNLOAD_USING =='https':
            self.url = "https://www.genome.jp/ftp/db/community/aaindex/aaindex1"
        else:
            self.url = "ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1"

        #if AAIndex database not found in data directory then call download method
        if not (os.path.isfile(os.path.join(self.aaindex_module_path, DATA_DIR, AAINDEX_FILENAME))):
            self.download_aaindex()

        #if parsed json of AAIndex already in file then read it and return 
        if (os.path.isfile(os.path.join(self.aaindex_module_path, DATA_DIR, AAINDEX_FILENAME + '.json'))):
            with open(os.path.join(self.aaindex_module_path, DATA_DIR, AAINDEX_FILENAME +'.json')) as aai_json:
                self.aaindex_json = json.load(aai_json)
        else:
            #parse AAIndex database file into JSON format
            self.aaindex_json = self.parse_aaindex()

        #get dict of categories
        self.categories = self.get_all_categories()

    def parse_aaindex(self):
        """
        Parse AAI database into JSON format. Each AAI record will be indexed by
        its feature code/index code, and will be in the format as shown in the
        docstring above. The file will be stored in a json file called according
        to (AAINDEX_FILENAME) variable.

        Returns
        -------
        :aaindex_json : dict
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
            tmp_filepath = os.path.join(self.aaindex_module_path, DATA_DIR, AAINDEX_FILENAME)
            f = open(tmp_filepath,'r')
        except IOError:
            print('Error opening file, check filename = {} and is stored in {} directory.'.format(AAINDEX_FILENAME, DATA_DIR))

        #read lines of file
        lines = f.readlines()
        f.close()

        #pattern that seperates each AAI record
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

                #get values for each row of amino acids in record
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

                #look for duplicate name entries
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
        
        #get categories for each record
        categories = self.get_all_categories()

        #append '-' to each aa index entry to account for missing AA in a protein sequence
        for index in aaindex_json:

            #add category to aaindex file
            aaindex_json[index]['category'] = categories[index]
            aaindex_json[index]['values']['-'] = 0

            #set any NA amino acid values to 0
            for val in aaindex_json[index]['values']:
                if aaindex_json[index]['values'][val] == 'NA':
                    aaindex_json[index]['values'][val] = 0

        #save parsed dictionary into JSON format to DATA_DIR
        with open((os.path.join(self.aaindex_module_path, DATA_DIR, AAINDEX_FILENAME + '.json')),'w') as output_F:
          json.dump(aaindex_json, output_F, indent=4, sort_keys=True)

        return aaindex_json

    def download_aaindex(self, save_dir=DATA_DIR):
        """
        If AAI database not found in DATA directory, then it will be downloaded from
        the dedicated FTP or HTTPS server from https://www.genome.jp/aaindex/.
        FTP is the default method used for downloading the database with the URL:
        "ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1". If you want to
        download by https, set the 'DOWNLOAD_USING' variable to 'https'. The class
        looks for the AAI database and its parsed JSON form in the 'data' directory, an
        error will throw if not found.

        Parameters
        ----------
        :save_dir : str (default = DATA_DIR)
            Directory to save the AAI database to. By default it wil be stored in
            the global var DATA_DIR = 'data'. Parameter has to point to a directory
            within the aaindex directory or there will be an error on instantiation 
            of the class. 
        """
        #fetch AAI database from URL if not present in DATA_DIR
        try:
            if not(os.path.isfile(os.path.join(self.aaindex_module_path, save_dir, AAINDEX_FILENAME))):
                try:
                    with closing(request.urlopen(self.url)) as r:
                        with open((os.path.join(self.aaindex_module_path, save_dir, AAINDEX_FILENAME)), 'wb') as f:
                            shutil.copyfileobj(r, f)
                    print('AAIndex1 successfully downloaded.')
                except requests.exceptions.RequestException:
                    print('Error downloading or exporting AAIndex from the url {}.'.format(self.url))
            else:
                pass    #AAIndex already in folder
        except:
            raise OSError('Save Directory does not exist: {}.'.format(save_dir))

    def parse_categories(self, aaindex_category_file='aaindex_to_category.txt'):
        """
        Parse category file which maps each AAI record in the databse into 1 of 8 categories.
        Category file and parsing code inspired from: https://github.com/harmslab/hops. 

        Parameters
        ----------
        :aaindex_category_file : str
            Name of category file to parse (default is "aaindex_to_category.txt"). If parameter
            is a file name then it will be imported from the 'data' directory by default, else
            it will be imported from its full path.

        Returns
        -------
        :aaindex_category : dict
            Dictionary that maps each AAI record into 1 of 8 categories.
        """
        #if input parameter is a full path, read it else read from default DATA_DIR
        if os.path.isfile(aaindex_category_file):
            try:
                f = open(aaindex_category_file,'r')
            except IOError:
                print('Error opening AAIndex1 category file {} .'
                    .format(aaindex_category_file))
        else:
            try:
                f = open((os.path.join(self.aaindex_module_path, DATA_DIR, aaindex_category_file)),'r')
            except IOError:
                print('Error opening AAIndex1 category file, check {} in {} directory'
                    .format(aaindex_category_file, DATA_DIR))

        #get total number of lines in file
        # total_lines = len(f.readlines(  ))

        #open new file in data directory to store parsed category file
        f.seek(0)
        category_output_file = "aaindex_categories.txt"
        f_out = open((os.path.join(self.aaindex_module_path, DATA_DIR, category_output_file)), "w")

        #lines starting with '#' are file metadata so don't write these to parsed output file
        for line in f.readlines():
            if not (line.startswith('#')):
                f_out.write(line)
        #close both files
        f.close()
        f_out.close()

        #open parsed category file for reading
        with open(os.path.join(self.aaindex_module_path, DATA_DIR, category_output_file)) as f_out:
            reader = csv.reader(f_out, delimiter="\t")
            d = list(reader)

        aaindex_category = {}

        #iterate through all lines in parsed file and store AAI indices as keys
        #  in the output dict and their respective categories as the values of dict
        for i in range(0, len(d)):
          for j in range(0, len(d[i])):
            category_substring = d[i][1].strip()
            category_substring = category_substring.split(" ", 1)
            aaindex_category[d[i][0]] = category_substring[0]

        self.categories = aaindex_category

    def get_all_categories(self, category_file="aaindex_categories.txt"):
        """
        Return list of all indices, their description and associated categories from 
        parsed aaindex_categories.txt file created from parse_categories() function.

        Parameters
        ----------
        :category_file : str
            path to categories mapping file.
        
        Returns
        -------
        :aaindex_category : dict
            Dictionary that maps each AAI record into 1 of 8 categories.
        """
        #if parsed categories file doesn't exist in DATA_DIR then call function
        if not (os.path.isfile(os.path.join(self.aaindex_module_path, DATA_DIR, category_file))):
            self.parse_categories()

        #read categories file and its content
        try:
            with open(os.path.join(self.aaindex_module_path, DATA_DIR, category_file)) as f_out:
                reader = csv.reader(f_out, delimiter="\t")
                d = list(reader)
        except IOError:
            print('Error opening AAIndex1 category file {}.'.format(category_file))

        aaindex_category = {}

        #iterate through all lines in parsed file and store AAI indices as keys
        #   in the output dict and their respective categories as the values of dict
        for i in range(0, len(d)):
          for j in range(0, len(d[i])):
            category_substring = d[i][1].strip()
            category_substring = category_substring.split(" ", 1)
            aaindex_category[d[i][0]] = category_substring[0]

        return aaindex_category

    def get_category_from_record(self, record_code):
        """
        Return category of a AAI database record from its 
        feaure/index code/accession number.

        Parameters
        ----------
        :record_code : str
            AAI database record feature code/index code.
        Returns
        -------
        :cat : str
            category of AAI database record according to parsed aaindex_to_category.txt file.
        """
        #stripping input of whitespace
        try:
            record_code = record_code.strip()
        except:
            raise TypeError('Input parameter {} is not of correct datatype string, got {}' \
                .format(record_code, type(record_code)))

        #check that inputted record_code does exist in the AAI database
        if record_code not in (self.record_codes()):
            raise ValueError('Record index {} not in AAI database'.format(record_code))

        #get category from dict
        cat = self.categories[record_code]

        return cat

    def get_record_from_desc(self, description, all_matches=False):
        """
        Return full AAI database record details from its description. Search
        through the descriptions/names of all records in the database, returning
        record that matches name input parameter. All matches can be returned or
        just the most relevant one, by default just the first match is found, this
        can be changed using the all_matches input param.

        Parameters
        ----------
        :description : str
            AAI database record feature name/description.
        :all_matches : bool (default=false)
            if set to true, all found matches in database returned else if false only the
            first found match is returned.

        Returns
        -------
        :matches : list(dict)
            closest found AAIndex record that matches user input description parameter,
            according to closeness function OR all found relevant AAIndex records.
        """
        matches = []
        all_desc = []

        #iterate through all records in AAI finding matching records with description
        for index, value in self.aaindex_json.items():
            all_desc.append(value['description'].lower())

        #use closeness function to calculate whether description mostly matches record description
        record_matches = get_close_matches(description.lower(), all_desc, cutoff=0.5)

        #if no matches found, return empty list
        if (record_matches==[]):
            return []
        #else iterate through all records in AAI, finding one with most matching description to
        #the input parameter, using the closeness function. Return most found record or all 
        #records, depending on all_matches parameter.
        else:
            for index, value in aaindex.aaindex_json.items():
                if not all_matches:
                    if (record_matches[0] == value['description']):
                        matches.append(self.aaindex_json[index])
                else:
                    for i in range(0,len(record_matches)):
                        if (record_matches[i] == value['description']):
                            matches.append(self.aaindex_json[index])    

        return matches


    def get_amino_acids(self):
        """
        Get all canonical amino acid letters. The '-' value will also be included
        in the list from this function as it accounts for the abcense of any
        amino acid or gaps in a AAI record.

        Returns
        -------
        :amino_acids : list
            List of all 20 canonical amino acid letters as found in each record
            of the AAI database.
        """
        amino_acids = list(self.aaindex_json[list(self.aaindex_json.keys())[0]]["values"].keys())
        amino_acids.sort()  #sort into alphabetical order
        
        return amino_acids

    def get_amino_acids_encoding(self):
        """
        Get one-hot encoding of amino acids.

        Returns
        -------
        :onehot_encoded : np.ndarray
            one hot encoded array of the 20 canonical amino acids.
        """
        all_amino_acids = self.get_amino_acids()
        #convert amino acid letters to np array
        values = np.array(all_amino_acids[1:])    

        #encode amino acids with value between 0 and n_classes-1.
        label_encoder = LabelEncoder()
        integer_encoded = label_encoder.fit_transform(values)

        #encode amino acids as a one-hot numeric array.
        onehot_encoder = OneHotEncoder(sparse=False)
        integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
        onehot_encoded = onehot_encoder.fit_transform(integer_encoded)

        return onehot_encoded

    def record_codes(self):
        """
        Get list of all AAI index codes/Accession numbers for each record in the database.

        Returns
        -------
        :records : list
            list of record names/index codes for all records in the AAI database.
        """
        records = list(self.aaindex_json.keys())
        records.sort()  #sort into alphabetical order

        return records

    def get_num_records(self):
        """
        Calculate total number of records/indices in the AAI database.

        Returns
        -------
        :len(aaindex.record_codes()) : int
            number of indices/records found in the AAI database.
        """
        return len(self.record_codes())

    def get_record_names(self):
        """
        Return a list of all index descriptions for all records in the AAI database.

        Returns
        -------
        :desc : list
            list of descriptions for all records in the AAI database.
        """
        desc = []

        #iterate through database, appending all descriptions to desc list
        for name in list(self.aaindex_json.values()):
            desc.append(name['description'])

        return desc

######################          Getters & Setters          ######################

    @property
    def url(self):
        return self._url

    @url.setter
    def url(self, value):
        self._url = value

    @property
    def categories(self):
        return self._categories

    @categories.setter
    def categories(self, value):
        self._categories = value

################################################################################

    def __repr__(self):
        return "<AAIndex: {}>".format(self)

    def __sizeof__(self):
        """ Return size of AAI database file """
        return os.path.getsize(os.path.isfile(os.path.join(DATA_DIR, AAINDEX_FILENAME)))

    def __getitem__(self, record_code):
        """
        Return full AAI database record details from its feature/index code/Accession number by
        making the whole database subscriptable.

        Parameters
        ----------
        :record_code : str
            AAI database record feature index/Accession number.

        Returns
        -------
        :record : dict
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

        return record

aaindex = AAIndex()

# get_ref_from_record()
# get_values_from_record