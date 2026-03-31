################################################################################
################                    AAindex1                   #################
################################################################################

#importing required modules and dependencies
import json
import os
import sys
import copy
import re
import csv
from typing import Dict, Iterator, List, Union

from ._aaindex_matrix import Map

__all__: List[str] = ['AAIndex1', 'aaindex1']


class AAIndex1():
    """Python parser for AAindex1: Amino Acid Index Database.

    The AAindex is a database of numerical indices representing various
    physicochemical and biochemical properties of amino acids. This class
    stores the amino acid index of 20 numerical values for the 20 amino
    acids — AAindex1 (http://www.genome.jp/aaindex/).

    Attributes:
        aaindex_module_path: Absolute path to the aaindex package directory.
        data_dir: Subdirectory name containing raw and cached data files.
        aaindex_filename: Base filename for this database (no extension).
        aaindex_json: Parsed database keyed by accession number.
        categories: Dict mapping each record code to its category.
        last_updated: Date string of the last published database update.
    """
    def __init__(self) -> None:
        #resolve the package directory for data file lookups
        self.aaindex_module_path = os.path.dirname(os.path.abspath(sys.modules[self.__module__].__file__))
        self.data_dir = "data"
        self.aaindex_filename = "aaindex1"

        #get dict of categories
        self.categories = self.get_all_categories()

        #load from cached JSON if available, otherwise parse the raw file
        json_path = os.path.join(self.aaindex_module_path, self.data_dir, f"{self.aaindex_filename}.json")
        if os.path.isfile(json_path):
            with open(json_path) as aai_json:
                self.aaindex_json = json.load(aai_json)
        else:
            self.aaindex_json = self.parse_aaindex()

        #date as shown on https://www.genome.jp/aaindex/
        self.last_updated = "February 13, 2017"

        #cache amino acid list once at init to avoid re-sorting on every call
        self._amino_acids_cache: List[str] = sorted(
            self.aaindex_json[next(iter(self.aaindex_json))]["values"].keys()
        )

    def parse_aaindex(self) -> Dict:
        """Parse the raw AAindex1 database file into a nested dict and cache as JSON.

        Each record is keyed by its accession number and stores metadata, amino
        acid values, and category. The result is written to a .json file in the
        data directory for fast subsequent loads.

        Returns:
            Parsed database keyed by accession number.

        Raises:
            IOError: If the raw database file cannot be opened.
            ValueError: If a duplicate accession number is encountered.
        """
        #initialise keys of AAi database
        template_dict = {
            "H": [], "D": [], "R": [], "A": [],
            "*": [], "T": [], "J": [], "C": [], "I": [],
        }

        #open AAi file for reading and parsing
        tmp_filepath = os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename)
        try:
            with open(tmp_filepath, 'r') as f:
                lines = f.readlines()
        except IOError as e:
            raise IOError(f"Error opening AAindex1 file, check file is in filepath: {tmp_filepath}.") from e

        #regex to normalise double-quote characters in field values
        clean_up_pattern = re.compile("\"")

        aaindex_json: Dict = {}
        current_dict = copy.deepcopy(template_dict)
        current_entry: str = "H"  # first non-space line in any block is always H

        #iterate through each line, parsing records delimited by '//'
        for line in lines:
            if line.startswith("//"):

                #handle meta data of each record
                name = clean_up_pattern.sub("'", " ".join(current_dict["H"]))
                description = clean_up_pattern.sub("'", " ".join(current_dict["D"]))

                #append author, title and journal name to reference
                a = " ".join(current_dict["A"])
                t = " ".join(current_dict["T"])
                j = " ".join(current_dict["J"])
                references = clean_up_pattern.sub("'", f"{a} '{t}' {j}")

                #parse pub med article ID
                pmid = clean_up_pattern.sub("'", " ".join(current_dict["R"]))
                pmid = pmid.replace("PMID:", "")

                #parse correlation coefficients into a dict
                correlation_coefficient = clean_up_pattern.sub("'", " ".join(current_dict["C"]))
                correlation_coefficient_ = {}
                correlation_coefficient_list = [
                    correlation_coefficient.split()[n:n + 2]
                    for n in range(0, len(correlation_coefficient.split()), 2)
                ]
                for correlation in correlation_coefficient_list:
                    correlation_coefficient_[correlation[0]] = correlation[1]

                #parse notes from record
                notes = clean_up_pattern.sub("'", " ".join(current_dict["*"]))

                #parse individual amino acid values from I-lines
                aa_lines = current_dict["I"]
                aa_names = aa_lines[0].split()
                row_0_names = [aa.split("/")[0] for aa in aa_names]
                row_1_names = [aa.split("/")[1] for aa in aa_names]
                row_0_values = aa_lines[1].split()
                row_1_values = aa_lines[2].split()

                values: Dict = {}
                for i in range(len(row_0_values)):
                    try:
                        values[row_0_names[i]] = float(row_0_values[i])
                    except ValueError:
                        values[row_0_names[i]] = "NA"
                    try:
                        values[row_1_names[i]] = float(row_1_values[i])
                    except ValueError:
                        values[row_1_names[i]] = "NA"

                #guard against duplicate accession numbers
                if name in aaindex_json:
                    raise ValueError(f"Duplicate AAi Record found: {name}.")

                aaindex_json[name] = {
                    "description": description,
                    "references": references,
                    "pmid": pmid,
                    "correlation_coefficients": correlation_coefficient_,
                    "notes": notes,
                    "values": values,
                }

                current_dict = copy.deepcopy(template_dict)
                continue

            if line[0] != " ":
                current_entry = line[0]

            current_dict[current_entry].append(line[1:].strip())

        #post-process: set NA values to 0, add category and '-' gap placeholder
        for index in aaindex_json:
            for val in aaindex_json[index]['values']:
                if aaindex_json[index]['values'][val] == 'NA':
                    aaindex_json[index]['values'][val] = 0
            aaindex_json[index]['category'] = self.categories[index]
            aaindex_json[index]['values']['-'] = 0

        #cache parsed database as JSON for fast subsequent loads
        json_out_path = os.path.join(
            self.aaindex_module_path, self.data_dir, f"{self.aaindex_filename}.json"
        )
        with open(json_out_path, 'w') as output_f:
            json.dump(aaindex_json, output_f, indent=4, sort_keys=True)

        return aaindex_json

    def parse_categories(self, aaindex_category_file: str = 'aaindex_to_category.txt') -> Dict:
        """Parse category file mapping each AAi record to one of 8 categories.

        Category file and parsing code inspired from https://github.com/harmslab/hops.

        Args:
            aaindex_category_file: Filename or full path of the category mapping
                file. Defaults to ``aaindex_to_category.txt`` in the data directory.

        Returns:
            Dict mapping each record code to its category string.

        Raises:
            IOError: If the category file cannot be opened.
        """
        #if input parameter is a full path, use it directly, else read from default 'data' dir
        if os.path.isfile(aaindex_category_file):
            category_filepath = aaindex_category_file
        else:
            category_filepath = os.path.join(self.aaindex_module_path, self.data_dir, aaindex_category_file)

        try:
            with open(category_filepath, 'r') as f:
                category_lines = f.readlines()
        except IOError as e:
            raise IOError(f"Error opening AAindex1 category file: {category_filepath}.") from e

        #write non-comment lines to parsed output file (strip '#' metadata lines)
        category_output_file = "aaindex_categories.txt"
        with open(os.path.join(self.aaindex_module_path, self.data_dir, category_output_file), "w") as f_out:
            for line in category_lines:
                if not line.startswith('#'):
                    f_out.write(line)

        #open parsed category file for reading
        with open(os.path.join(self.aaindex_module_path, self.data_dir, category_output_file)) as f_out:
            reader = csv.reader(f_out, delimiter="\t")
            lines = list(reader)

        aaindex_category: Dict = {}

        #iterate through all lines, map each record code to its category
        for row in lines:
            if len(row) >= 2:
                category_substring = row[1].strip().split(" ", 1)
                aaindex_category[row[0]] = category_substring[0]

        return aaindex_category

    def get_all_categories(self, category_file: str = "aaindex_categories.txt") -> Dict:
        """Return dict mapping every record code to its category.

        Reads from the parsed ``aaindex_categories.txt`` file produced by
        :meth:`parse_categories`. If the file does not yet exist, it is
        generated first.

        Args:
            category_file: Filename of the pre-parsed categories file
                inside the data directory.

        Returns:
            Dict mapping each record code to its category string.

        Raises:
            IOError: If the categories file cannot be opened.
        """
        #if parsed categories file doesn't exist in 'data' then call function to get it
        if not os.path.isfile(os.path.join(self.aaindex_module_path, self.data_dir, category_file)):
            self.parse_categories()

        #read categories file and its content
        cat_filepath = os.path.join(self.aaindex_module_path, self.data_dir, category_file)
        try:
            with open(cat_filepath) as f_out:
                reader = csv.reader(f_out, delimiter="\t")
                d = list(reader)
        except IOError as e:
            raise IOError(f"Error opening AAindex1 category file: {cat_filepath}.") from e

        aaindex_category: Dict = {}

        #iterate through all lines, map each record code to its category
        for row in d:
            if len(row) >= 2:
                category_substring = row[1].strip().split(" ", 1)
                aaindex_category[row[0]] = category_substring[0]

        return aaindex_category

    def search(self, description: Union[str, List[str]]) -> Dict:
        """Search records by keyword(s) present in their description field.

        Args:
            description: Keyword string or list of keyword strings.
                         Matching is case-insensitive.

        Returns:
            Dict of matching records keyed by accession number.
            Returns an empty dict if no records match.

        Raises:
            TypeError: If description is not a str or list.
        """
        all_indices: Dict = {}

        if not isinstance(description, (list, str)):
            raise TypeError(f"Input description parameter must be a list or str, got {type(description)}.")

        #convert description parameter to list to make iterable
        if not isinstance(description, list):
            description = [description]

        #iterate over description list, matching keywords case-insensitively
        for desc in description:
            for index, value in self.aaindex_json.items():
                if desc.lower() in self.aaindex_json[index]['description'].lower():
                    all_indices[index] = self.aaindex_json[index]

        return all_indices

    def amino_acids(self) -> List[str]:
        """Return sorted list of amino acid single-letter codes.

        Includes the ``-`` placeholder for absent/gap amino acids.

        Returns:
            Sorted list of amino acid codes including ``-``.
        """
        #return pre-computed cache from __init__
        return self._amino_acids_cache

    def record_codes(self) -> List[str]:
        """Return sorted list of all accession numbers in the database.

        Returns:
            Sorted list of accession number strings.
        """
        records = list(self.aaindex_json.keys())
        records.sort()
        return records

    def num_records(self) -> int:
        """Return the total number of records in the database.

        Returns:
            Number of records as int.
        """
        return len(self.aaindex_json)

    def record_names(self) -> List[str]:
        """Return a list of description strings for all records.

        Returns:
            List of description strings in database insertion order.
        """
        return [v['description'] for v in self.aaindex_json.values()]

    def values(self, record_code: str) -> Dict:
        """Return the amino acid values dict for a given record.

        Shortcut to avoid accessing the full record when only
        the values are needed.

        Args:
            record_code: AAindex accession number.

        Returns:
            Dict of amino acid values for the specified record.

        Raises:
            ValueError: If record_code is not found in the database.
        """
        return self[record_code]['values']

    def get_record_by_category(self, category: str) -> Dict:
        """Return all records belonging to a given category.

        Args:
            category: Category name to filter records by (case-insensitive).

        Returns:
            Dict of matching records keyed by accession number.

        Raises:
            TypeError: If category is not a string.
        """
        if not isinstance(category, str):
            raise TypeError(f"Input category parameter must be a str, got {type(category)}.")

        #filter records by matching category field
        category_records = {
            code: record for code, record in self.aaindex_json.items()
            if record.get('category', '').lower() == category.lower()
        }
        return category_records

    def plot(self, record_code: str) -> None:
        """Display a bar chart of the 20 amino acid values for a given record.

        Requires matplotlib to be installed. The ``-`` gap placeholder is
        excluded from the plot.

        Args:
            record_code: AAindex accession number.

        Raises:
            ImportError: If matplotlib is not installed.
            ValueError: If record_code is not found in the database.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as e:
            raise ImportError(
                "matplotlib is required for plot(). "
                "Install it with: pip install matplotlib"
            ) from e

        record = self[record_code]
        vals = record.values
        #exclude the '-' gap placeholder from the plot
        aa_list = sorted(k for k in vals if k != '-')
        scores = [vals[aa] for aa in aa_list]

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.bar(aa_list, scores)
        ax.set_xlabel("Amino Acid")
        ax.set_ylabel("Value")
        ax.set_title(f"{record_code}: {record.description}")
        plt.tight_layout()
        plt.show()

    def __getitem__(self, record_code: str) -> "Map":
        """Return a record by accession number wrapped in a Map (dot-notation dict).

        Args:
            record_code: AAindex accession number (case-insensitive,
                         leading/trailing whitespace is stripped).

        Returns:
            Record data as a Map, accessible via dict or dot notation.

        Raises:
            TypeError: If record_code is not a string.
            ValueError: If record_code is not found in the database.
        """
        try:
            record_code = record_code.strip().upper()
        except AttributeError:
            raise TypeError(
                f"Input parameter {record_code} is not of correct datatype string, got {type(record_code)}."
            )

        if record_code not in self.aaindex_json:
            raise ValueError(f"Record Index ({record_code}) not found in AAindex1.")

        return Map(self.aaindex_json[record_code])

    def __sizeof__(self) -> int:
        """Return the on-disk size of the raw AAindex data file in bytes."""
        return os.path.getsize(
            os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename)
        )

    def __len__(self) -> int:
        """Return total number of records in the database."""
        return len(self.aaindex_json)

    def __contains__(self, record_code: object) -> bool:
        """Return True if record_code exists in the database."""
        return record_code in self.aaindex_json

    def __iter__(self) -> Iterator[str]:
        """Iterate over all record codes in the database."""
        return iter(self.aaindex_json)

    def __repr__(self) -> str:
        """Return a canonical string representation of this instance."""
        return f"AAIndex1(records={len(self.aaindex_json)}, last_updated='{self.last_updated}')"

######################          Getters & Setters          ######################

    @property
    def categories(self) -> Dict:
        return self._categories

    @categories.setter
    def categories(self, value: Dict) -> None:
        self._categories = value

    @property
    def data_dir(self) -> str:
        return self._data_dir

    @data_dir.setter
    def data_dir(self, value: str) -> None:
        self._data_dir = value

    @property
    def aaindex_filename(self) -> str:
        return self._aaindex_filename

    @aaindex_filename.setter
    def aaindex_filename(self, value: str) -> None:
        self._aaindex_filename = value

    @property
    def last_updated(self) -> str:
        return self._last_updated

    @last_updated.setter
    def last_updated(self, value: str) -> None:
        self._last_updated = value


#create instance of AAIndex1 class
aaindex1 = AAIndex1()