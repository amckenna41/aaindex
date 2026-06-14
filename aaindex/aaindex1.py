################################################################################
################                    AAindex1                   #################
################################################################################

#importing required modules and dependencies
import difflib
import json
import math
import os
import sys
import copy
import re
import csv
import warnings
from importlib.metadata import version as _pkg_version, PackageNotFoundError
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

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
        #date as shown on https://www.genome.jp/aaindex/
        self.last_updated = "February 13, 2017"
        #database and categories are loaded lazily on first access
        self._aaindex_json: Optional[Dict] = None
        self._categories: Optional[Dict] = None
        self._amino_acids_cache: Optional[List[str]] = None
        try:
            self.version = _pkg_version("aaindex")
        except PackageNotFoundError:
            self.version = "unknown"

    def _load_data(self) -> None:
        """Load the database from the JSON cache, or parse the raw file if the cache is absent."""
        json_path = os.path.join(self.aaindex_module_path, self.data_dir, f"{self.aaindex_filename}.json")
        if os.path.isfile(json_path):
            with open(json_path) as aai_json:
                self._aaindex_json = json.load(aai_json)
            #normalise types that may differ in older cached JSON
            for record in self._aaindex_json.values():
                cc = record.get('correlation_coefficients', {})
                for k, v in list(cc.items()):
                    if isinstance(v, str):
                        try:
                            cc[k] = float(v)
                        except ValueError:
                            pass
                vals = record.get('values', {})
                if '-' in vals and isinstance(vals['-'], int):
                    vals['-'] = float(vals['-'])
        else:
            self._aaindex_json = self._parse_aaindex()
        self._amino_acids_cache = sorted(
            self._aaindex_json[next(iter(self._aaindex_json))]["values"].keys()
        )

    def _parse_aaindex(self) -> Dict:
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
            with open(tmp_filepath) as f:
                lines = f.readlines()
        except OSError as e:
            raise OSError(f"Error opening AAindex1 file, check file is in filepath: {tmp_filepath}.") from e

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

                #parse correlation coefficients into a dict (values stored as float)
                correlation_coefficient = clean_up_pattern.sub("'", " ".join(current_dict["C"]))
                correlation_coefficient_ = {}
                correlation_coefficient_list = [
                    correlation_coefficient.split()[n:n + 2]
                    for n in range(0, len(correlation_coefficient.split()), 2)
                ]
                for correlation in correlation_coefficient_list:
                    if len(correlation) == 2:
                        try:
                            correlation_coefficient_[correlation[0]] = float(correlation[1])
                        except ValueError:
                            pass

                #parse notes from record
                notes = clean_up_pattern.sub("'", " ".join(current_dict["*"]))

                #parse individual amino acid values from I-lines
                aa_lines = current_dict["I"]
                if len(aa_lines) < 3:
                    warnings.warn(
                        f"Record {name!r} has a malformed I-block ({len(aa_lines)} line(s)); skipping.",
                        UserWarning, stacklevel=2,
                    )
                    current_dict = copy.deepcopy(template_dict)
                    continue
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
                        values[row_0_names[i]] = None
                    try:
                        values[row_1_names[i]] = float(row_1_values[i])
                    except ValueError:
                        values[row_1_names[i]] = None

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

        #post-process: assign category and '-' gap placeholder (NA stays as None)
        for index in aaindex_json:
            aaindex_json[index]['category'] = self.categories[index]
            aaindex_json[index]['values']['-'] = 0.0

        #cache parsed database as JSON for fast subsequent loads; non-fatal if the data dir is read-only
        json_out_path = os.path.join(
            self.aaindex_module_path, self.data_dir, f"{self.aaindex_filename}.json"
        )
        try:
            with open(json_out_path, 'w') as output_f:
                json.dump(aaindex_json, output_f, indent=4, sort_keys=True)
        except OSError:
            pass

        return aaindex_json

    def parse_aaindex(self) -> Dict:
        """Deprecated: use _parse_aaindex() instead.

        .. deprecated::
            This method is an internal implementation detail and will be
            removed in a future version.
        """
        warnings.warn(
            "parse_aaindex() is deprecated and will be removed in a future version. "
            "It is an internal implementation detail.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._parse_aaindex()

    def _parse_categories(self, aaindex_category_file: str = 'aaindex_to_category.txt') -> Dict:
        """Parse category file mapping each AAi record to one of 8 categories.

        Category file and parsing code inspired from https://github.com/harmslab/hops.

        Accepts only a bare filename (resolved against the data directory) or an
        absolute path. Relative paths containing directory components are rejected
        to prevent path traversal.

        Args:
            aaindex_category_file: Bare filename or absolute path of the category
                mapping file. Defaults to ``aaindex_to_category.txt`` in the data
                directory.

        Returns:
            Dict mapping each record code to its category string.

        Raises:
            ValueError: If *aaindex_category_file* contains path separators and is
                not an absolute path.
            IOError: If the category file cannot be opened.
        """
        #S1 fix: reject relative paths with directory components to prevent traversal
        if os.path.isabs(aaindex_category_file):
            category_filepath = aaindex_category_file
        elif os.path.basename(aaindex_category_file) != aaindex_category_file:
            raise ValueError(
                f"aaindex_category_file must be a bare filename or an absolute path; "
                f"relative paths with directory components are not allowed, "
                f"got: {aaindex_category_file!r}"
            )
        else:
            category_filepath = os.path.join(self.aaindex_module_path, self.data_dir, aaindex_category_file)

        try:
            with open(category_filepath) as f:
                category_lines = f.readlines()
        except OSError as e:
            raise OSError(f"Error opening AAindex1 category file: {category_filepath}.") from e

        #B7 fix: filter comment lines in memory — no intermediate file write needed
        aaindex_category: Dict = {}
        reader = csv.reader(
            (line for line in category_lines if not line.startswith('#')),
            delimiter="\t",
        )
        for row in reader:
            if len(row) >= 2:
                category_substring = row[1].strip().split(" ", 1)
                aaindex_category[row[0]] = category_substring[0]

        return aaindex_category

    def parse_categories(self, aaindex_category_file: str = 'aaindex_to_category.txt') -> Dict:
        """Deprecated: use _parse_categories() instead.

        .. deprecated::
            This method is an internal implementation detail and will be
            removed in a future version.
        """
        warnings.warn(
            "parse_categories() is deprecated and will be removed in a future version. "
            "It is an internal implementation detail.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._parse_categories(aaindex_category_file)

    def get_all_categories(self, category_file: str = "aaindex_categories.txt") -> Dict:
        """Return dict mapping every record code to its category.

        Delegates to :meth:`_parse_categories` which parses the source category
        file in memory. The *category_file* parameter is retained for API
        compatibility but is no longer used.

        Returns:
            Dict mapping each record code to its category string.
        """
        return self._parse_categories()

    def search(self, query: Union[str, List[str]]) -> Dict:
        """Search records by keyword(s) across all text fields.

        Searches the description, accession code, PMID, references, and notes
        fields of each record. Matching is case-insensitive.

        Args:
            query: Keyword string or list of keyword strings.

        Returns:
            Dict of matching records keyed by accession number.
            Returns an empty dict if no records match.

        Raises:
            TypeError: If query is not a str or list.
        """
        all_indices: Dict = {}

        if not isinstance(query, (list, str)):
            raise TypeError(f"Input query parameter must be a list or str, got {type(query)}.")

        #convert query parameter to list to make iterable
        if not isinstance(query, list):
            query = [query]

        #iterate over query list, matching keywords case-insensitively across all text fields
        for term in query:
            term_lower = term.lower()
            for index, record in self.aaindex_json.items():
                if (
                    term_lower in record['description'].lower()
                    or term_lower in index.lower()
                    or term_lower in record['pmid'].lower()
                    or term_lower in record['references'].lower()
                    or term_lower in record['notes'].lower()
                ):
                    all_indices[index] = record

        return all_indices

    def search_fuzzy(
        self,
        query: str,
        n: int = 10,
        cutoff: float = 0.0,
    ) -> Dict:
        """Search records using fuzzy matching with ranked results.

        Scores each record by how closely *query* matches any of its text
        fields (description, accession code, PMID, references, notes, category)
        using ``difflib.SequenceMatcher``. Results are returned in descending
        score order (most relevant first) with no external dependencies.

        Args:
            query: Search string.
            n: Maximum number of results to return. Defaults to ``10``.
            cutoff: Minimum similarity score in [0, 1] to include a record.
                    Defaults to ``0.0`` (all records scored, top *n* returned).

        Returns:
            Dict of matching records keyed by accession number, ordered by
            descending similarity score.

        Raises:
            TypeError: If *query* is not a str.
            ValueError: If *n* < 1 or *cutoff* is outside [0, 1].
        """
        if not isinstance(query, str):
            raise TypeError(f"query must be a str, got {type(query)}.")
        if not isinstance(n, int) or n < 1:
            raise ValueError(f"n must be a positive int, got {n}.")
        if not (0.0 <= cutoff <= 1.0):
            raise ValueError(f"cutoff must be in [0, 1], got {cutoff}.")
        query_lower = query.lower()
        scored: List[Tuple[float, str]] = []
        for code, record in self.aaindex_json.items():
            candidate_fields = [
                record.get("description", ""),
                code,
                record.get("pmid", ""),
                record.get("references", ""),
                record.get("notes", ""),
                record.get("category", ""),
            ]
            score = max(
                difflib.SequenceMatcher(None, query_lower, field.lower()).ratio()
                for field in candidate_fields
                if field
            )
            if score >= cutoff:
                scored.append((score, code))
        scored.sort(key=lambda x: x[0], reverse=True)
        return {code: self.aaindex_json[code] for _, code in scored[:n]}

    def amino_acids(self) -> List[str]:
        """Return sorted list of amino acid single-letter codes.

        Includes the ``-`` placeholder for absent/gap amino acids.

        Returns:
            Sorted list of amino acid codes including ``-``.
        """
        #trigger lazy load if necessary (also sets _amino_acids_cache)
        _ = self.aaindex_json
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

    def valid_categories(self) -> List[str]:
        """Return a sorted list of all distinct category strings in the database.

        Returns:
            Sorted list of valid category name strings.
        """
        return sorted({record['category'] for record in self.aaindex_json.values()})

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

        #validate the category against known values; raise with a helpful hint if unknown
        valid = self.valid_categories()
        valid_lower = [c.lower() for c in valid]
        if category.lower() not in valid_lower:
            close = difflib.get_close_matches(category.lower(), valid_lower, n=3, cutoff=0.6)
            if close:
                close_actual = [valid[valid_lower.index(c)] for c in close]
                hint = f" Did you mean: {close_actual}?"
            else:
                hint = f" Valid categories are: {valid}."
            raise ValueError(f"Category '{category}' not found in AAindex1.{hint}")

        #filter records by matching category field
        category_records = {
            code: record for code, record in self.aaindex_json.items()
            if record.get('category', '').lower() == category.lower()
        }
        return category_records

    def encode(self, sequence: str, record_code: str, gap_value: float = 0.0) -> List[Optional[float]]:
        """Encode a protein sequence as a numeric feature vector using one index.

        Each amino acid in *sequence* is replaced by its corresponding value
        from the named AAindex1 record. Gaps (``-``) are replaced by
        *gap_value*. Amino acids with an NA value in the source data are
        represented as ``None``.

        Args:
            sequence: Protein sequence string of single-letter amino acid codes.
                      Case-insensitive; non-standard characters are allowed
                      and map to ``None``.
            record_code: AAindex1 accession number.
            gap_value: Numeric value substituted for gap characters (``-``).
                       Defaults to ``0.0``.

        Returns:
            List of floats (or ``None``) with the same length as *sequence*.

        Raises:
            TypeError: If *sequence* is not a string.
            ValueError: If *record_code* is not found in the database.

        Example:
            >>> aaindex1.encode('ACDE', 'KYTJ820101')
            [1.8, 2.5, -3.5, -3.5]
        """
        if not isinstance(sequence, str):
            raise TypeError(f"sequence must be a str, got {type(sequence)}.")
        values = self[record_code]['values']
        result: List[Optional[float]] = []
        for aa in sequence.upper():
            if aa == '-':
                result.append(gap_value)
            else:
                result.append(values.get(aa))
        return result

    def batch_encode(
        self,
        sequences: List[str],
        record_code: str,
        gap_value: float = 0.0,
    ) -> List[List[Optional[float]]]:
        """Encode multiple protein sequences as numeric feature vectors using one index.

        Convenience wrapper around :meth:`encode` for processing a collection of
        sequences with a single AAindex1 record in one call. The record is looked
        up once and reused for all sequences, making this more efficient than
        calling :meth:`encode` in a loop.

        Args:
            sequences: List of protein sequence strings of single-letter amino
                       acid codes. Case-insensitive; non-standard characters map
                       to ``None``.
            record_code: AAindex1 accession number.
            gap_value: Numeric value substituted for gap characters (``-``).
                       Defaults to ``0.0``.

        Returns:
            List of encoded sequences; each inner list has the same length as the
            corresponding input sequence and contains floats or ``None``.

        Raises:
            TypeError: If *sequences* is not a list, or any element is not a str.
            ValueError: If *record_code* is not found in the database.

        Example:
            >>> aaindex1.batch_encode(['ACDE', 'FGHI'], 'KYTJ820101')
            [[1.8, 2.5, -3.5, -3.5], [-3.2, -0.4, -3.2, -3.2]]
        """
        if not isinstance(sequences, list):
            raise TypeError(f"sequences must be a list, got {type(sequences)}.")
        values = self[record_code]['values']
        result: List[List[Optional[float]]] = []
        for i, seq in enumerate(sequences):
            if not isinstance(seq, str):
                raise TypeError(f"sequences[{i}] must be a str, got {type(seq)}.")
            encoded: List[Optional[float]] = []
            for aa in seq.upper():
                if aa == '-':
                    encoded.append(gap_value)
                else:
                    encoded.append(values.get(aa))
            result.append(encoded)
        return result

    def get_correlated_indices(
        self,
        record_code: str,
        min_correlation: float = 0.8,
        depth: int = 1,
    ) -> Dict[str, float]:
        """Return a dict of indices correlated with *record_code*.

        Traverses the ``correlation_coefficients`` links stored on each record
        (pre-computed by the AAindex authors for |r| ≥ 0.8 pairs) up to
        *depth* hops from the seed record. Only indices whose stored |r| meets
        *min_correlation* are included.

        Args:
            record_code: Seed AAindex1 accession number.
            min_correlation: Minimum absolute correlation coefficient to
                             include a neighbour. Must be in [0, 1].
                             Defaults to ``0.8``.
            depth: Number of hops to traverse from the seed. ``1`` returns
                   only direct neighbours; ``2`` also includes their
                   neighbours, etc. Defaults to ``1``.

        Returns:
            Dict mapping each correlated accession number to its stored
            correlation coefficient (float). The seed record itself is
            excluded from the result.

        Raises:
            TypeError: If *min_correlation* or *depth* are of the wrong type.
            ValueError: If *record_code* is not found, or if *min_correlation*
                        is outside [0, 1], or if *depth* is less than 1.

        Example:
            >>> aaindex1.get_correlated_indices('KYTJ820101', depth=1)
            {'EISD840101': '0.949', ...}
        """
        if not isinstance(min_correlation, (int, float)):
            raise TypeError(f"min_correlation must be a float, got {type(min_correlation)}.")
        if not isinstance(depth, int):
            raise TypeError(f"depth must be an int, got {type(depth)}.")
        if not (0.0 <= min_correlation <= 1.0):
            raise ValueError(f"min_correlation must be in [0, 1], got {min_correlation}.")
        if depth < 1:
            raise ValueError(f"depth must be >= 1, got {depth}.")

        #validate seed — __getitem__ raises ValueError if not found
        self[record_code]

        visited: Dict[str, float] = {}
        frontier = {record_code.strip().upper()}

        for _ in range(depth):
            next_frontier: set = set()
            for code in frontier:
                if code not in self.aaindex_json:
                    continue
                for neighbour, coeff_str in self.aaindex_json[code]['correlation_coefficients'].items():
                    try:
                        coeff = float(coeff_str)
                    except (ValueError, TypeError):
                        continue
                    if abs(coeff) >= min_correlation and neighbour not in visited:
                        visited[neighbour] = coeff
                        next_frontier.add(neighbour)
            frontier = next_frontier - {record_code.strip().upper()}

        #exclude the seed from output
        seed_upper = record_code.strip().upper()
        visited.pop(seed_upper, None)
        return visited

    def compare_indices(
        self,
        code1: str,
        code2: str,
    ) -> Optional[float]:
        """Compute the Pearson correlation coefficient between two AAindex1 records.

        Uses only the amino acids that have a non-``None`` value in *both*
        records (i.e. NA positions are excluded pairwise). Returns ``None``
        when fewer than two valid positions remain after masking.

        Args:
            code1: First AAindex1 accession number.
            code2: Second AAindex1 accession number.

        Returns:
            Pearson r as a float in [-1, 1], or ``None`` if there are
            insufficient non-NA values to compute a correlation.

        Raises:
            ValueError: If either accession number is not found in the database.

        Example:
            >>> aaindex1.compare_indices('KYTJ820101', 'EISD840101')
            0.949...
        """
        vals1 = self[code1]['values']
        vals2 = self[code2]['values']

        #collect paired values — both must be non-None and not the gap placeholder
        pairs = [
            (vals1[aa], vals2[aa])
            for aa in vals1
            if aa != '-' and vals1[aa] is not None and vals2.get(aa) is not None
        ]

        n = len(pairs)
        if n < 2:
            return None

        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]

        mean_x = sum(xs) / n
        mean_y = sum(ys) / n

        cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
        std_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
        std_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))

        if std_x == 0.0 or std_y == 0.0:
            return None

        return cov / (std_x * std_y)

    def to_dict(self, record_code: Optional[str] = None) -> Dict:
        """Export one record or the full database as a plain Python dict.

        Args:
            record_code: If given, export only that record keyed by its
                         accession number. If ``None`` (default), export
                         the entire database.

        Returns:
            Dict containing the requested records.

        Raises:
            ValueError: If *record_code* is not found in the database.
        """
        if record_code is not None:
            code = record_code.strip().upper()
            if code not in self.aaindex_json:
                raise ValueError(f"Record ({code}) not found in AAindex1.")
            return {code: dict(self.aaindex_json[code])}
        return dict(self.aaindex_json)

    def to_json(self, record_code: Optional[str] = None, indent: int = 4) -> str:
        """Serialise one record or the full database to a JSON string.

        Args:
            record_code: If given, serialise only that record. If ``None``
                         (default), serialise the entire database.
            indent: JSON indentation level. Defaults to ``4``.

        Returns:
            JSON-formatted string.

        Raises:
            ValueError: If *record_code* is not found in the database.
        """
        return json.dumps(self.to_dict(record_code), indent=indent)

    def to_dataframe(self, record_code: Optional[str] = None) -> Any:
        """Export amino acid values as a pandas DataFrame.

        Requires ``pandas`` to be installed. If pandas is not available a
        clear ``ImportError`` is raised rather than a silent failure.

        When *record_code* is ``None`` the returned DataFrame has one column
        per amino acid and one row per AAindex1 record (index = accession
        numbers). When a specific record code is given the DataFrame has a
        single row.

        Args:
            record_code: If given, export only that record's values.
                         If ``None`` (default), export all records.

        Returns:
            ``pandas.DataFrame`` with amino acids as columns and record
            accession numbers as the index.

        Raises:
            ImportError: If ``pandas`` is not installed.
            ValueError: If *record_code* is not found in the database.
        """
        try:
            import pandas as pd  # noqa: PLC0415
        except ImportError as exc:
            raise ImportError(
                "pandas is required for to_dataframe(). "
                "Install it with: pip install pandas"
            ) from exc

        if record_code is not None:
            code = record_code.strip().upper()
            records_subset = {code: self[code]['values']}
        else:
            records_subset = {
                code: record['values'] for code, record in self.aaindex_json.items()
            }

        return pd.DataFrame.from_dict(records_subset, orient='index')

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
            close = difflib.get_close_matches(
                record_code, self.aaindex_json.keys(), n=3, cutoff=0.6
            )
            hint = f" Did you mean: {close}?" if close else ""
            raise ValueError(f"Record Index ({record_code}) not found in AAindex1.{hint}")

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
    def aaindex_json(self) -> Dict:
        """Parsed database dict, loaded lazily on first access."""
        if self._aaindex_json is None:
            self._load_data()
        return self._aaindex_json

    @aaindex_json.setter
    def aaindex_json(self, value: Dict) -> None:
        self._aaindex_json = value

    @property
    def categories(self) -> Dict:
        if self._categories is None:
            self._categories = self.get_all_categories()
        return self._categories

    @categories.setter
    def categories(self, value: Dict) -> None:
        self._categories = value

    @property
    def data_dir(self) -> str:
        return self._data_dir

    @data_dir.setter
    def data_dir(self, value: str) -> None:
        if os.path.basename(value) != value:
            raise ValueError(
                f"data_dir must be a simple directory name with no path separators, got: {value!r}"
            )
        self._data_dir = value

    @property
    def aaindex_filename(self) -> str:
        return self._aaindex_filename

    @aaindex_filename.setter
    def aaindex_filename(self, value: str) -> None:
        if os.path.basename(value) != value:
            raise ValueError(
                f"aaindex_filename must be a simple filename with no path separators, got: {value!r}"
            )
        self._aaindex_filename = value

    @property
    def last_updated(self) -> str:
        return self._last_updated

    @last_updated.setter
    def last_updated(self, value: str) -> None:
        self._last_updated = value


#create instance of AAIndex1 class
aaindex1 = AAIndex1()