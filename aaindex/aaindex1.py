################################################################################
#################                    AAIndex                   #################
################################################################################

#importing required modules and dependencies
import json
import numpy as np
import os
import sys, copy, re
import requests
import shutil
from difflib import get_close_matches
import csv
import urllib.request as request
from contextlib import closing

#global vars
DOWNLOAD_USING = 'ftp'

class AAIndex1():
    """
            Python parser for AAindex1: Amino Acid Index Database (AAI)
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
    * R Pub med article ID (PMID)                                          *
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
    amino_acids():
        get list of 20 amino acid letters.
    record_codes():
        get list of AAI record codes/Accession numbers.
    num_records():
        return total number of records in AAI database.
    record_names():
        return list of all descriptions for all records in AAI database.
    search():
        return 1 or more AAI records from its description.
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
        self.data_dir = "data"
        self.aaindex_filename = "aaindex1"
        self.download_using = DOWNLOAD_USING

        #download AAI database using ftp or https
        if self.download_using =='ftp':
            self.url = "ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1"
        elif self.download_using =='http' or self.download_using =='https':
            self.url = "https://www.genome.jp/ftp/db/community/aaindex/aaindex1"
        else:
            self.url = "ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1"

        #if AAIndex database not found in data directory then call download method
        if not (os.path.isfile(os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename))):
            self.download_aaindex()

        #get dict of categories
        self.categories = self.get_all_categories()

        #if parsed json of AAIndex already in file then read it and return 
        if (os.path.isfile(os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename + '.json'))):
            with open(os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename + '.json')) as aai_json:
                self.aaindex_json = json.load(aai_json)
        else:
            #parse AAIndex database file into JSON format
            self.aaindex_json = self.parse_aaindex()

        #get last updated - as shown on homepage (https://www.genome.jp/aaindex/)
        self.last_updated = "February 13, 2017" 

    def parse_aaindex(self):
        """
        Parse AAI database into JSON format. Each AAI record will be indexed by
        its accession number/index code, and will be in the format as shown in the
        docstring above, with the addition of a category key which shows the type/
        category of the record in question. The file will be stored in a json file 
        determined by AAINDEX_FILENAME variable.

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

        #open AAI file for reading and parsing, by default it should be stored in self.data_dir
        try:
            tmp_filepath = os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename)
            f = open(tmp_filepath,'r')
        except IOError:
            print('Error opening file, check filename = {} and is stored in {} directory.'.format(
                    self.aaindex_filename, self.data_dir))

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

                # handle meta data of each record - name, desc, citation, notes
                name = " ".join(current_dict["H"])
                name = clean_up_pattern.sub("'", name)

                description = " ".join(current_dict["D"])
                description = clean_up_pattern.sub("'", description)

                #append author, title and journal name to reference
                references = "{} '{}' {}".format(" ".join(current_dict["A"]),
                                               " ".join(current_dict["T"]),
                                               " ".join(current_dict["J"]))
                references = clean_up_pattern.sub("'", references)

                # references = references + "; Kawashima, S. and Kanehisa, M."
                #     "AAindex: amino acid index database.'  Nucleic Acids Res. 28, 374 (2000)."

                #parse pub med article ID
                pmid = " ".join(current_dict["R"])
                pmid = clean_up_pattern.sub("'", pmid)
                pmid = pmid.replace("PMID:", "")

                #appending correlation coefficient of current record, convert into dict
                correlation_coefficient = " ".join(current_dict["C"])
                correlation_coefficient = clean_up_pattern.sub("'", correlation_coefficient)
                correlation_coefficient_ = {}
                correlation_coefficient_list = [correlation_coefficient.split()[n:n+2] for n in range(0, len(correlation_coefficient.split()), 2)]

                for correlation in correlation_coefficient_list:
                    correlation_coefficient_[correlation[0]] = correlation[1]
            
                #parse notes from record
                notes = " ".join(current_dict["*"])
                notes = clean_up_pattern.sub("'", notes)

                # parse individual amino acid values
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

                #append all aaindex record data to object
                aaindex_json[name] = {"description": description,
                                  "references": references,
                                  "pmid": pmid,
                                  "correlation_coefficients": correlation_coefficient_,
                                  "notes": notes,
                                  "values": values
                            }
                
                #create copy of dict 
                current_dict = copy.deepcopy(template_dict)
                continue

            this_entry = l[0]
            if l[0] != " ":
                current_entry = this_entry

            current_dict[current_entry].append(l[1:].strip())
        
        #append '-' to each aa index entry to account for missing AA in a protein sequence
        for index in aaindex_json:

            #add category to aaindex file
            aaindex_json[index]['category'] = self.categories[index]
            aaindex_json[index]['values']['-'] = 0

            #set any NA amino acid values to 0
            for val in aaindex_json[index]['values']:
                if aaindex_json[index]['values'][val] == 'NA':
                    aaindex_json[index]['values'][val] = 0

        #save parsed dictionary into JSON format to self.data_dir
        with open((os.path.join(self.aaindex_module_path, self.data_dir, 
            self.aaindex_filename + '.json')),'w') as output_F:
          json.dump(aaindex_json, output_F, indent=4, sort_keys=True)

        return aaindex_json

    def download_aaindex(self, save_dir="data"):
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
        :save_dir : str (default = 'data')
            Directory to save the AAI database to. By default it wil be stored in
            the 'data' dir. Parameter has to point to a directory within the aaindex 
            directory or there will be an error on instantiation of the class. 
        """
        #fetch AAI database from URL if not present in 'data' dir
        try:
            if not(os.path.isfile(os.path.join(self.aaindex_module_path, save_dir, self.aaindex_filename))):
                try:
                    with closing(request.urlopen(self.url)) as r:
                        with open((os.path.join(self.aaindex_module_path, save_dir, self.aaindex_filename)), 'wb') as f:
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
        #if input parameter is a full path, read it else read from default 'data'
        if os.path.isfile(aaindex_category_file):
            try:
                f = open(aaindex_category_file,'r')
            except IOError:
                print('Error opening AAIndex1 category file {}.'
                    .format(aaindex_category_file))
        else:
            try:
                f = open((os.path.join(self.aaindex_module_path, self.data_dir, aaindex_category_file)),'r')
            except IOError:
                print('Error opening AAIndex1 category file, check {} is in {} directory.'
                    .format(aaindex_category_file, self.data_dir))

        #get total number of lines in file
        # total_lines = len(f.readlines(  ))

        #open new file in data directory to store parsed category file
        f.seek(0)
        category_output_file = "aaindex_categories.txt"
        f_out = open((os.path.join(self.aaindex_module_path, self.data_dir, category_output_file)), "w")

        #lines starting with '#' are file metadata so don't write these to parsed output file
        for line in f.readlines():
            if not (line.startswith('#')):
                f_out.write(line)

        #close both files
        f.close()
        f_out.close()

        #open parsed category file for reading
        with open(os.path.join(self.aaindex_module_path, self.data_dir, category_output_file)) as f_out:
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
        
        for code, cat in aaindex_category.items():
            self.aaindex_json[code]["category"] = cat

        return aaindex_category

    def get_all_categories(self, category_file="aaindex_categories.txt"):
        """
        Return list of all indices, their description and associated categories from 
        parsed aaindex_categories.txt file created from parse_categories() function.

        Parameters
        ----------
        :category_file : str (default = "aaindex_categories.txt")
            path to categories mapping file.
        
        Returns
        -------
        :aaindex_category : dict
            Dictionary that maps each AAI record into 1 of 8 categories.
        """
        #if parsed categories file doesn't exist in 'data' then call function to get it
        if not (os.path.isfile(os.path.join(self.aaindex_module_path, self.data_dir, category_file))):
            self.parse_categories()

        #read categories file and its content
        try:
            with open(os.path.join(self.aaindex_module_path, self.data_dir, category_file)) as f_out:
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

    def search(self, description):
        """
        Return 1 or more full AAI database records from their description. Search
        through the descriptions/names of all records in the database, returning
        record(s) that match description input parameter.

        Parameters
        ----------
        :description : str
            keywords to search for in all AAI1 records.

        Returns
        -------
        :all_indices : dict
            AAIndex record(s) that contain the keywords in the description input parameter.
        """
        all_indices = {}

        #if input param not a list or string then raisse type error
        if not (isinstance(description, list)) and not (isinstance(description, str)):
            raise TypeError("Input Description parameter must be a list or str, got {}.".format(type(description)))

        #convert description parameter to list to make iterable
        if not (isinstance(description, list)):
            description = [description]

        #iterate over description list, if keywords are in an aaindex record's description
        #   then add record to all_indices object
        for desc in description:
            for index, value in self.aaindex_json.items():
                if (desc.lower() in self.aaindex_json[index]['description'].lower()):
                    all_indices[index] = self.aaindex_json[index]

        return all_indices

    def amino_acids(self):
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

    def num_records(self):
        """
        Calculate total number of records/indices in the AAI database.

        Returns
        -------
        :len(aaindex.record_codes()) : int
            number of indices/records found in the AAI database.
        """
        return len(self.record_codes())

    def record_names(self):
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

    @property
    def last_updated(self):
        return self._last_updated

    @last_updated.setter
    def last_updated(self, value):
        self._last_updated = value

################################################################################

    # def __repr__(self):
    #     return "<AAIndex: {}>".format(self)

    def __sizeof__(self):
        """ Return size of AAI database file """
        return os.path.getsize(os.path.isfile(os.path.join(self.data_dir, self.aaindex_filename)))

    def __getitem__(self, record_code):
        """
        Return full AAI database record details from its index code/Accession number by
        making the whole database subscriptable.

        Parameters
        ----------
        :record_code : str
            AAI database record Accession number.

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

        #create instance of map class which allows you to access dict using dot notation 
        record = Map(record)

        return record

class Map(dict):
    """
    Instantiating this class will convert a dict such that it can be accessed using 
    dot notation which makes it easier for accessing the individual elements of the 
    aaindex records. It also works for nested dicts.

    Parameters 
    ----------
    :dict : dict 
        input dictionary to be mapped into dot notation.

    Usage
    -----
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
   
    Reference
    ---------
    [1]: https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
    """
    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]

#create instance of aaindex1 class
aaindex1 = AAIndex1()