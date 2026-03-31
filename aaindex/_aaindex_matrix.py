################################################################################
################         AAindex Matrix Base Class             #################
################################################################################

#importing required modules and dependencies
import json
import os
import sys
import copy
import re
from typing import Dict, Iterator, List, Optional, Union


class Map(dict):
    """A dict subclass that enables attribute-style (dot notation) access to keys.

    Works for nested dicts. Each AAindex record returned by __getitem__ is
    wrapped in this class so fields can be read as record.description as well
    as record['description'].

    References:
        https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v
        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError(f"'Map' object has no attribute '{attr}'")

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super().__delitem__(key)
        del self.__dict__[key]

    def __repr__(self) -> str:
        return f"Map({dict.__repr__(self)})"


class _AAIndexMatrix:
    """Base class for AAindex2 and AAindex3 matrix database parsers.

    Provides shared parsing, lookup, search, and protocol methods for the
    lower-triangular 20x20 matrix databases. Subclasses call
    super().__init__(filename) with the appropriate base filename so the
    correct data file is loaded.

    Attributes:
        aaindex_module_path: Absolute path to the aaindex package directory.
        data_dir: Subdirectory name containing raw and cached data files.
        aaindex_filename: Base filename for this database (no extension).
        aaindex_json: Parsed database keyed by accession number.
        last_updated: Date string of the last published database update.
    """

    def __init__(self, filename: str) -> None:
        #resolve the package directory for data file lookups
        self.aaindex_module_path = os.path.dirname(
            os.path.abspath(sys.modules[self.__module__].__file__)
        )
        self.data_dir = "data"
        self.aaindex_filename = filename

        #load from cached JSON if available, otherwise parse the raw file
        json_path = os.path.join(
            self.aaindex_module_path, self.data_dir, self.aaindex_filename + ".json"
        )
        if os.path.isfile(json_path):
            with open(json_path) as aai_json:
                self.aaindex_json = json.load(aai_json)
        else:
            self.aaindex_json = self.parse_aaindex()

        #date as shown on https://www.genome.jp/aaindex/
        self.last_updated = "February 13, 2017"

    def parse_aaindex(self) -> Dict:
        """Parse the raw AAindex database file into a nested dict and cache as JSON.

        Each record is keyed by its accession number and stores metadata
        alongside the full symmetric 20x20 matrix reconstructed from the
        lower-triangular source data. The result is written to a .json file in
        the data directory for fast subsequent loads.

        Returns:
            dict: Parsed database keyed by accession number.

        Raises:
            IOError: If the raw database file cannot be opened.
            ValueError: If a duplicate accession number is encountered.
        """
        #template for each record's metadata fields
        template_dict = {
            "H": [], "D": [], "R": [], "A": [],
            "*": [], "T": [], "J": [], "C": [], "M": [],
        }

        tmp_filepath = os.path.join(
            self.aaindex_module_path, self.data_dir, self.aaindex_filename
        )
        try:
            with open(tmp_filepath) as f:
                lines = f.readlines()
        except OSError as e:
            raise OSError(
                f"Error opening {self.aaindex_filename} file, "
                f"check it is present at: {tmp_filepath}."
            ) from e

        #regex to normalise double-quote characters in field values
        clean_up_pattern = re.compile("\"")

        aaindex_json: Dict = {}
        current_dict = copy.deepcopy(template_dict)
        current_entry: str = "H"  # first non-space line in any block is always H

        for line in lines:
            if line.startswith("//"):
                #parse plain metadata fields
                name = clean_up_pattern.sub("'", " ".join(current_dict["H"]))
                description = clean_up_pattern.sub("'", " ".join(current_dict["D"]))
                a = " ".join(current_dict["A"])
                t = " ".join(current_dict["T"])
                j = " ".join(current_dict["J"])
                references = clean_up_pattern.sub("'", f"{a} '{t}' {j}")
                pmid = (
                    clean_up_pattern.sub("'", " ".join(current_dict["R"]))
                    .replace("PMID:", "")
                    .replace("LIT:", "")
                    .strip()
                )
                notes = clean_up_pattern.sub("'", " ".join(current_dict["*"]))

                #parse correlation coefficients into a dict
                corr_str = clean_up_pattern.sub("'", " ".join(current_dict["C"]))
                correlation_coefficients: Dict = {}
                corr_pairs = [
                    corr_str.split()[n:n + 2]
                    for n in range(0, len(corr_str.split()), 2)
                ]
                for pair in corr_pairs:
                    correlation_coefficients[pair[0]] = pair[1]

                #parse the M block into a full symmetric matrix
                matrix: Dict = {}
                row_order: List[str] = []
                col_order: List[str] = []
                matrix_rows: List[List] = []

                for m_line in current_dict["M"]:
                    stripped = m_line.strip()
                    if stripped.startswith("rows"):
                        #format: rows = <AA_STRING> cols = <AA_STRING>
                        parts = stripped.replace(",", "").split()
                        row_order = list(parts[2])
                        col_order = list(parts[5])
                    else:
                        row_vals: List = []
                        for token in stripped.split():
                            if token in ("NA", "-"):
                                row_vals.append(None)
                            else:
                                try:
                                    row_vals.append(float(token))
                                except ValueError:
                                    row_vals.append(None)
                        if row_vals:
                            matrix_rows.append(row_vals)

                #reconstruct full symmetric matrix from lower-triangular data
                for i, row_aa in enumerate(row_order):
                    if i >= len(matrix_rows):
                        break
                    matrix.setdefault(row_aa, {})
                    for j, val in enumerate(matrix_rows[i]):
                        if j >= len(col_order):
                            break
                        col_aa = col_order[j]
                        matrix[row_aa][col_aa] = val
                        #fill the symmetric counterpart
                        matrix.setdefault(col_aa, {})
                        matrix[col_aa][row_aa] = val

                if name in aaindex_json:
                    raise ValueError(f"Duplicate accession number found: {name}.")

                aaindex_json[name] = {
                    "description": description,
                    "references": references,
                    "pmid": pmid,
                    "correlation_coefficients": correlation_coefficients,
                    "notes": notes,
                    "matrix": matrix,
                    "row_order": row_order,
                    "col_order": col_order,
                }
                current_dict = copy.deepcopy(template_dict)
                continue

            #route non-separator lines to their field bucket
            if line[0] != " ":
                current_entry = line[0]
            current_dict[current_entry].append(line[1:].strip())

        #cache parsed database as JSON for fast subsequent loads
        json_out_path = os.path.join(
            self.aaindex_module_path, self.data_dir, self.aaindex_filename + ".json"
        )
        with open(json_out_path, "w") as output_f:
            json.dump(aaindex_json, output_f, indent=4, sort_keys=True)

        return aaindex_json

    def get(self, record_code: str, aa1: str, aa2: str) -> Optional[float]:
        """Return the pairwise matrix score for two amino acids from a given record.

        The matrix is symmetric so get(code, aa1, aa2) == get(code, aa2, aa1).
        Returns None when either amino acid carries an NA value in the source data
        or when the amino acid letter is not present in this record's matrix.

        Args:
            record_code: AAindex accession number.
            aa1: Single-letter code for the first amino acid.
            aa2: Single-letter code for the second amino acid.

        Returns:
            Pairwise score as float, or None if data is not available.

        Raises:
            TypeError: If aa1 or aa2 are not strings.
            ValueError: If record_code is not found in the database.
        """
        record = self[record_code]
        try:
            aa1 = aa1.strip().upper()
            aa2 = aa2.strip().upper()
        except AttributeError:
            raise TypeError("aa1 and aa2 must be single-letter string amino acid codes.")
        matrix = record.matrix
        if aa1 in matrix and aa2 in matrix[aa1]:
            return matrix[aa1][aa2]
        return None

    def values(self, record_code: str) -> Dict:
        """Return the full 20x20 matrix dict for a given record.

        Shortcut to avoid accessing the whole record when only the matrix is
        needed. Consistent with AAIndex1.values() which returns amino acid values.

        Args:
            record_code: AAindex accession number.

        Returns:
            Nested dict of pairwise scores keyed by single-letter amino acid codes.

        Raises:
            ValueError: If record_code is not found in the database.
        """
        return self[record_code].matrix

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
            raise TypeError(
                f"description must be a list or str, got {type(description)}."
            )
        if not isinstance(description, list):
            description = [description]
        for desc in description:
            for index, value in self.aaindex_json.items():
                if desc.lower() in value["description"].lower():
                    all_indices[index] = value
        return all_indices

    def amino_acids(self) -> List[str]:
        """Return sorted list of the 20 canonical amino acid single-letter codes.

        Derived from the row_order field of the first record in the database.

        Returns:
            Sorted list of single-letter amino acid codes.
        """
        first_record = self.aaindex_json[next(iter(self.aaindex_json))]
        return sorted(first_record["row_order"])

    def record_codes(self) -> List[str]:
        """Return sorted list of all accession numbers in the database.

        Returns:
            Sorted list of accession number strings.
        """
        return sorted(self.aaindex_json.keys())

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
        return [v["description"] for v in self.aaindex_json.values()]

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
                f"record_code must be a string, got {type(record_code)}."
            )
        if record_code not in self.aaindex_json:
            raise ValueError(
                f"Record ({record_code}) not found in {self.__class__.__name__}."
            )
        return Map(self.aaindex_json[record_code])

    def __len__(self) -> int:
        """Return total number of records in the database."""
        return len(self.aaindex_json)

    def __contains__(self, record_code: object) -> bool:
        """Return True if record_code exists in the database."""
        return record_code in self.aaindex_json

    def __iter__(self) -> Iterator[str]:
        """Iterate over all accession numbers in the database."""
        return iter(self.aaindex_json)

    def __repr__(self) -> str:
        """Return a canonical string representation of this instance."""
        return (
            f"{self.__class__.__name__}("
            f"records={len(self.aaindex_json)}, "
            f"last_updated='{self.last_updated}')"
        )

    def __sizeof__(self) -> int:
        """Return the on-disk size of the raw AAindex data file in bytes."""
        return os.path.getsize(
            os.path.join(self.aaindex_module_path, self.data_dir, self.aaindex_filename)
        )

######################          Getters & Setters          ######################

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
