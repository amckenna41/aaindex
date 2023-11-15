################################################################################
################                    AAindex1                   #################
################################################################################

#importing required modules and dependencies
import json
import os
import sys, copy, re
import csv

class AAIndex1():
    """
    Python parser for AAindex1: Amino Acid Index Database.
    
    The AAindex is a database of numerical indices representing various physiochemical
    and biochemical properties of amino acids and pairs of amino acids. This class 
    stores the amino acid index of 20 numerical values for the 20 amino acids - 
    AAindex1 (http://www.genome.jp/aaindex/).

    Data format of AAi1:
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

    Parameters
    ==========
    None
    
    Methods
    =======
    parse_aaindex():
        parse AAi database into JSON format.
    parse_categories(aaindex_category_file):
        parse and map AAi indices to their respective categories.
    amino_acids():
        get list of 20 canonical amino acid letters.
    record_codes():
        get list of all AAi record codes/Accession numbers in AA1 database.
    num_records():
        return total number of records currently in AAi1 database.
    record_names():
        return list of all descriptions for all records in AAi1 database.
    search():
        return 1 or more AAi records from its description.
    __getitem__():
        access full AAi record using its record code/accession number.
    
    References
    ==========
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
    [5]: https://github.com/harmslab/hops
    """
    def __init__(self):
        
        #initialise class variables/attributes
        self.aaindex_module_path = os.path.dirname(os.path.abspath(sys.modules[self.__module__].__file__))
        self.data_dir = "data"
        self.aaindex_filename = "aaindex1"

        #get dict of categories
        self.categories = self.get_all_categories()

        #if parsed json of AAindex (aaindex1.json) already in file then read it and return 
        if (os.path.isfile(os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename + '.json'))):
            with open(os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename + '.json')) as aai_json:
                self.aaindex_json = json.load(aai_json)
        else:
            #parse AAindex database file into JSON format
            self.aaindex_json = self.parse_aaindex()

        #get last updated - as shown on homepage (https://www.genome.jp/aaindex/)
        self.last_updated = "February 13, 2017" 

    def parse_aaindex(self):
        """
        Parse AAi database into JSON format. Each AAi record will be indexed by
        its accession number/index code, and will be in the format as shown in the
        docstring above, with the addition of a category key which shows the type/
        category of the record in question. The file will be stored in the data 
        folder in the aaindex module path in json format.
        
        Parameters
        ==========
        None
    
        Returns
        =======
        :aaindex_json: dict
          parsed AAi database in dict form.
        """
        #initialise keys of AAi database
        template_dict = {"H":[],
                         "D":[],
                         "R":[],
                         "A":[],
                         "*":[],
                         "T":[],
                         "J":[],
                         "C":[],
                         "I":[]}

        #open AAi file for reading and parsing, by default it should be stored in self.data_dir
        try:
            tmp_filepath = os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename)
            f = open(tmp_filepath,'r')
        except IOError:
            print('Error opening AAindex1 file, check file is in filepath: {}.'.format(
                    os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename)))

        #read lines of file
        lines = f.readlines()
        f.close()

        #pattern that seperates each AAi record
        clean_up_pattern = re.compile("\"")

        #initilaise parsed AAi database dictionary
        aaindex_json = {}
        current_dict = copy.deepcopy(template_dict)

        '''
        iterate through each line in the AAi database, parsing each record and its
        amino acid values into its own entry in the dictionary. Each index/record
        is seperated by a '//'. Remove any duplicate records, set any missing
        ('-') or NA amino acid values to 0. Store resulting dict into aaindex_json
        instance variable.
        '''
        for l in lines:
            if l.startswith("//"):

                #handle meta data of each record - name, desc, citation, pmid, notes
                name = " ".join(current_dict["H"])
                name = clean_up_pattern.sub("'", name)

                description = " ".join(current_dict["D"])
                description = clean_up_pattern.sub("'", description)

                #append author, title and journal name to reference
                references = "{} '{}' {}".format(" ".join(current_dict["A"]),
                                               " ".join(current_dict["T"]),
                                               " ".join(current_dict["J"]))
                references = clean_up_pattern.sub("'", references)

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

                #parse individual amino acid values
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

                #look for duplicate name entries - raise value error if already present
                if (name in aaindex_json):
                    raise ValueError('Duplicate AAi Record found: {}.'.format(name))

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

            #set any NA amino acid values to 0
            for val in aaindex_json[index]['values']:
                if aaindex_json[index]['values'][val] == 'NA':
                    aaindex_json[index]['values'][val] = 0

            #add category to aaindex file
            aaindex_json[index]['category'] = self.categories[index]
            aaindex_json[index]['values']['-'] = 0

        #save parsed dictionary into JSON format to self.data_dir
        with open((os.path.join(self.aaindex_module_path, self.data_dir, 
            self.aaindex_filename + '.json')),'w') as output_F:
          json.dump(aaindex_json, output_F, indent=4, sort_keys=True)

        return aaindex_json

    def parse_categories(self, aaindex_category_file='aaindex_to_category.txt'):
        """
        Parse category file which maps each AAi record in the databse into 1 of 8 categories.
        Category file and parsing code inspired from: https://github.com/harmslab/hops. 

        Parameters
        ==========
        :aaindex_category_file: str
            Name of category file to parse (default is "aaindex_to_category.txt"). If parameter
            is a filename then it will be imported from the 'data' directory by default, else
            it will be imported from its full path.

        Returns
        =======
        :aaindex_category: dict
            Dictionary that maps each AAi record into 1 of 8 categories.
        """
        #if input parameter is a full path, read it else read from default 'data' dir
        if os.path.isfile(aaindex_category_file):
            try:
                f = open(aaindex_category_file,'r')
            except IOError:
                print('Error opening AAindex1 category file: {}.'.format(aaindex_category_file))
        else:
            try:
                f = open((os.path.join(self.aaindex_module_path, self.data_dir, aaindex_category_file)),'r')
            except IOError:
                print('Error opening AAindex1 category file: {}.'.format(os.path.join(self.aaindex_module_path, self.data_dir, aaindex_category_file)))

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
            lines = list(reader)

        aaindex_category = {}

        #iterate through all lines in parsed file and store AAi indices as keys
        #  in the output dict and their respective categories as the values of dict
        for i in range(0, len(lines)):
          for j in range(0, len(lines[i])):
            category_substring = lines[i][1].strip()
            category_substring = category_substring.split(" ", 1)
            aaindex_category[lines[i][0]] = category_substring[0]
        
        #add cateogry to main aaindex json dict
        for code, cat in aaindex_category.items():
            self.aaindex_json[code]["category"] = cat

        return aaindex_category

    def get_all_categories(self, category_file="aaindex_categories.txt"):
        """
        Return list of all indices, their description and associated categories from 
        parsed aaindex_categories.txt file created from parse_categories() function.

        Parameters
        ==========
        :category_file: str (default="aaindex_categories.txt")
            path to categories mapping file.
        
        Returns
        =======
        :aaindex_category: dict
            dictionary that maps each AAi record into 1 of 8 categories.
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
            print('Error opening AAindex1 category file: {}.'.format(os.path.join(self.aaindex_module_path, self.data_dir, category_file)))

        aaindex_category = {}

        #iterate through all lines in parsed file and store AAi indices as keys
        # in the output dict and their respective categories as the values of dict
        for i in range(0, len(d)):
          for j in range(0, len(d[i])):
            category_substring = d[i][1].strip()
            category_substring = category_substring.split(" ", 1)
            aaindex_category[d[i][0]] = category_substring[0]

        return aaindex_category

    def search(self, description):
        """
        Return 1 or more full AAi database record(s) from their description. Search
        through the descriptions/names of all records in the database, returning
        record(s) that match description input parameter. If no matching records
        found then an empty object will be returned.

        Parameters
        ==========
        :description: list/str
            keywords to search for in all AAi1 records. Can be a single string
            or list of keywords to search for.

        Returns
        =======
        :all_indices: dict
            AAindex record(s) that contain the keywords in the description input parameter.
            If no matching records found then an empty object will be returned.
        """
        all_indices = {}

        #if input param not a list or string then raisse type error
        if not (isinstance(description, list)) and not (isinstance(description, str)):
            raise TypeError("Input description parameter must be a list or str, got {}.".format(type(description)))

        #convert description parameter to list to make iterable
        if not (isinstance(description, list)):
            description = [description]

        #iterate over description list, if keywords are in an aaindex record's description
        # then add record to all_indices object
        for desc in description:
            for index, value in self.aaindex_json.items():
                if (desc.lower() in self.aaindex_json[index]['description'].lower()):
                    all_indices[index] = self.aaindex_json[index]

        return all_indices

    def amino_acids(self):
        """
        Get list of all canonical amino acid letters. The '-' value will also be included
        in the list from this function as it accounts for the abcense of any amino acid 
        or gaps in a AAi record.

        Parameters
        ==========
        None

        Returns
        =======
        :amino_acids: list
            List of all 20 canonical amino acid letters as found in each record
            of the AAi database.
        """
        amino_acids = list(self.aaindex_json[list(self.aaindex_json.keys())[0]]["values"].keys())
        amino_acids.sort()  #sort into alphabetical order
        
        return amino_acids

    def record_codes(self):
        """
        Get list of all AAi index codes/Accession numbers for each record in the database.

        Parameters
        ==========

        Returns
        =======
        :records: list
            list of record names/index codes for all records in the AAi database.
        """
        records = list(self.aaindex_json.keys())
        records.sort()  #sort into alphabetical order

        return records

    def num_records(self):
        """
        Calculate total number of records/indices in the AAi database.

        Parameters
        ==========
        None

        Returns
        =======
        :len(aaindex.record_codes()): int
            number of indices/records found in the AAi database.
        """
        return len(self.record_codes())

    def record_names(self):
        """
        Return a list of all index descriptions for all records in the AAi database.

        Parameters
        ==========
        None

        Returns
        =======
        :desc: list
            list of descriptions for all records in the AAi database.
        """
        desc = []

        #iterate through database, appending all descriptions to desc list
        for name in list(self.aaindex_json.values()):
            desc.append(name['description'])

        return desc

    def __getitem__(self, record_code):
        """
        Return full AAi database record details from its index code/Accession number by
        making the whole database subscriptable. Raise error if invalid accession number
        input.

        Parameters
        ==========
        :record_code: str
            AAi database record Accession number.

        Returns
        =======
        :record: dict
            dict of AAi database record.

        Usage
        =====
        aaindex = AAIndex()
        full_record = aaindex['ZIMJ680105']
        """
        #stripping input of whitespace, uppercase
        try:
            record_code = record_code.strip().upper()
        except:
            raise TypeError('Input parameter {} is not of correct datatype string, got {}.' \
                .format(record_code, type(record_code)))

        #check that inputted record_code/Accession num does exist in the AAi database
        if record_code not in (self.record_codes()):
            raise ValueError('Record Index ({}) not found in AAindex1.'.format(record_code))
        
        #get full record from parsed JSON
        record = self.aaindex_json[record_code]

        #create instance of map class which allows you to access dict using dot notation 
        record = Map(record)

        return record

    def __sizeof__(self):
        """ Return size of AAi database file """
        return os.path.getsize(os.path.isfile(os.path.join(self.data_dir, self.aaindex_filename)))
    
######################          Getters & Setters          ######################

    @property
    def categories(self):
        return self._categories

    @categories.setter
    def categories(self, value):
        self._categories = value

    @property
    def data_dir(self):
        return self._data_dir

    @data_dir.setter
    def data_dir(self, value):
        self._data_dir = value

    @property
    def aaindex_filename(self):
        return self._aaindex_filename

    @aaindex_filename.setter
    def aaindex_filename(self, value):
        self._aaindex_filename = value

    @property
    def last_updated(self):
        return self._last_updated

    @last_updated.setter
    def last_updated(self, value):
        self._last_updated = value

################################################################################
class Map(dict):
    """
    Instantiating this class will convert a dict such that it can be accessed using 
    dot notation which makes it easier for accessing the individual elements of the 
    aaindex records. It also works for nested dicts.

    Parameters 
    ==========
    :dict: dict 
        input dictionary to be mapped into dot notation.

    Usage
    =====
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    # Add new key
    m.new_key = 'Hello world!'
    # Or
    m['new_key'] = 'Hello world!'
    print m.new_key
    print m['new_key']
    # Update values
    m.new_key = 'Yay!'
    # Or
    m['new_key'] = 'Yay!'
    # Delete key
    del m.new_key
    # Or
    del m['new_key']
    
    References
    ==========
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