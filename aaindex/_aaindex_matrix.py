################################################################################
################         AAindex Matrix Base Class             #################
################################################################################

#importing required modules and dependencies
import difflib
import json
import os
import sys
import copy
import re
import warnings
from importlib.metadata import version as _pkg_version, PackageNotFoundError
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union


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
        self.__dict__.pop(key, None)

    def __repr__(self) -> str:
        """Return a structured, human-readable representation of the record."""
        summary_keys = ("description", "pmid", "notes", "category", "is_symmetric")
        lines = []
        for key in summary_keys:
            if key in self:
                val = self[key]
                if isinstance(val, str) and len(val) > 80:
                    val_repr = repr(val[:77] + "...")
                else:
                    val_repr = repr(val)
                lines.append(f"  {key}: {val_repr}")
        for key in ("values", "matrix", "correlation_coefficients"):
            if key in self:
                val = self[key]
                if isinstance(val, dict):
                    lines.append(f"  {key}: {{...}} ({len(val)} entries)")
                else:
                    lines.append(f"  {key}: {val!r}")
        if "references" in self:
            ref = self["references"]
            if isinstance(ref, str) and len(ref) > 80:
                lines.append(f"  references: {repr(ref[:77] + '...')}")
            else:
                lines.append(f"  references: {ref!r}")
        return "AAIndexRecord(\n" + "\n".join(lines) + "\n)"

    def __str__(self) -> str:
        """Return a human-readable string representation of the record."""
        return self.__repr__()


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
        #date as shown on https://www.genome.jp/aaindex/
        self.last_updated = "February 13, 2017"
        #database is loaded lazily on first access to avoid unnecessary I/O at import time
        self._aaindex_json: Optional[Dict] = None
        self._amino_acids_cache: Optional[List[str]] = None
        try:
            self.version = _pkg_version("aaindex")
        except PackageNotFoundError:
            self.version = "unknown"

    def _load_data(self) -> None:
        """Load the database from the JSON cache, or parse the raw file if the cache is absent."""
        json_path = os.path.join(
            self.aaindex_module_path, self.data_dir, self.aaindex_filename + ".json"
        )
        if os.path.isfile(json_path):
            with open(json_path) as aai_json:
                self._aaindex_json = json.load(aai_json)
            #normalise correlation_coefficients to float; older cached JSON stored them as strings
            for record in self._aaindex_json.values():
                cc = record.get('correlation_coefficients', {})
                for k, v in list(cc.items()):
                    if isinstance(v, str):
                        try:
                            cc[k] = float(v)
                        except ValueError:
                            pass
        else:
            self._aaindex_json = self._parse_aaindex()

    def _parse_aaindex(self) -> Dict:
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
                    if len(pair) == 2:
                        try:
                            correlation_coefficients[pair[0]] = float(pair[1])
                        except ValueError:
                            pass

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
                        if len(parts) >= 6:
                            row_order = list(parts[2])
                            col_order = list(parts[5])
                        else:
                            warnings.warn(
                                f"Malformed rows/cols header in record {name!r}: {stripped!r}",
                                UserWarning, stacklevel=2,
                            )
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
                #detect whether matrix is lower-triangular (symmetric) or full (non-symmetric)
                #a lower-triangular matrix has row i containing i+1 values;
                #a full matrix has every row containing len(col_order) values
                is_full = bool(
                    matrix_rows
                    and col_order
                    and len(matrix_rows[0]) == len(col_order)
                )

                for i, row_aa in enumerate(row_order):
                    if i >= len(matrix_rows):
                        warnings.warn(
                            f"Matrix for record {name!r} has fewer rows than expected "
                            f"({len(matrix_rows)} vs {len(row_order)}); matrix may be incomplete.",
                            UserWarning, stacklevel=2,
                        )
                        break
                    matrix.setdefault(row_aa, {})
                    for j, val in enumerate(matrix_rows[i]):
                        if j >= len(col_order):
                            warnings.warn(
                                f"Matrix row {row_aa!r} in record {name!r} has fewer columns "
                                f"than expected ({j} vs {len(col_order)}); matrix may be incomplete.",
                                UserWarning, stacklevel=2,
                            )
                            break
                        col_aa = col_order[j]
                        matrix[row_aa][col_aa] = val
                        if not is_full:
                            #lower-triangular: fill the symmetric counterpart
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
                    "is_symmetric": not is_full,
                }
                current_dict = copy.deepcopy(template_dict)
                continue

            #route non-separator lines to their field bucket
            if line[0] != " ":
                current_entry = line[0]
            current_dict[current_entry].append(line[1:].strip())

        #cache parsed database as JSON for fast subsequent loads; non-fatal if the data dir is read-only
        json_out_path = os.path.join(
            self.aaindex_module_path, self.data_dir, self.aaindex_filename + ".json"
        )
        try:
            with open(json_out_path, "w") as output_f:
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

    def get(self, record_code: str, aa1: str, aa2: str) -> Optional[float]:
        """Return the pairwise matrix score for two amino acids from a given record.

        For symmetric (lower-triangular) records, get(code, aa1, aa2) ==
        get(code, aa2, aa1). For non-symmetric records, order matters:
        get(code, aa1, aa2) may differ from get(code, aa2, aa1).
        Returns None when the amino acid pair has an NA value in the source data
        or when the amino acid letter is not present in this record's matrix.

        Args:
            record_code: AAindex accession number.
            aa1: Single-letter code for the first amino acid (row).
            aa2: Single-letter code for the second amino acid (column).

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

    def search(self, query: Union[str, List[str]]) -> Dict:
        """Search records by keyword(s) across all text fields.

        Searches description, accession code, PMID, references, and notes.
        Matching is case-insensitive.

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
            raise TypeError(
                f"query must be a list or str, got {type(query)}."
            )
        if not isinstance(query, list):
            query = [query]
        for term in query:
            term_lower = term.lower()
            for index, value in self.aaindex_json.items():
                if (
                    term_lower in value["description"].lower()
                    or term_lower in index.lower()
                    or term_lower in value["pmid"].lower()
                    or term_lower in value["references"].lower()
                    or term_lower in value["notes"].lower()
                ):
                    all_indices[index] = value
        return all_indices

    def search_fuzzy(
        self,
        query: str,
        n: int = 10,
        cutoff: float = 0.0,
    ) -> Dict:
        """Search records using fuzzy matching with ranked results.

        Scores each record by how closely *query* matches any of its text
        fields (description, accession code, PMID, references, notes) using
        ``difflib.SequenceMatcher``. Results are returned in descending score
        order (most relevant first) with no external dependencies.

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
        """Return sorted list of the 20 canonical amino acid single-letter codes.

        Derived from the row_order field of the first record in the database.
        Result is cached after the first call.

        Returns:
            Sorted list of single-letter amino acid codes.
        """
        if self._amino_acids_cache is None:
            first_record = self.aaindex_json[next(iter(self.aaindex_json))]
            self._amino_acids_cache = sorted(first_record["row_order"])
        return self._amino_acids_cache

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
                raise ValueError(
                    f"Record ({code}) not found in {self.__class__.__name__}."
                )
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
        """Export pairwise matrix scores as a pandas DataFrame.

        Requires ``pandas`` to be installed. If pandas is not available a
        clear ``ImportError`` is raised rather than a silent failure.

        When *record_code* is given the returned DataFrame is the full 20×20
        (or N×N) matrix for that record, with row and column labels from
        ``row_order`` and ``col_order``. When *record_code* is ``None`` a
        MultiIndex DataFrame is returned with index levels
        ``(accession, row_aa)`` and column labels from ``col_order``.

        Args:
            record_code: If given, export only that record's matrix.
                         If ``None`` (default), export all records.

        Returns:
            ``pandas.DataFrame`` representing the requested matrix data.

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
            record = self[code]
            return pd.DataFrame.from_dict(record['matrix'], orient='index')

        #build a MultiIndex DataFrame across all records
        rows = []
        index_tuples = []
        first_record = self.aaindex_json[next(iter(self.aaindex_json))]
        col_order = first_record.get('col_order', sorted(first_record.get('matrix', {}).keys()))
        for code, record in self.aaindex_json.items():
            matrix = record.get('matrix', {})
            for row_aa in record.get('row_order', sorted(matrix.keys())):
                row_data = matrix.get(row_aa, {})
                rows.append([row_data.get(col_aa) for col_aa in col_order])
                index_tuples.append((code, row_aa))
        multi_index = pd.MultiIndex.from_tuples(index_tuples, names=['record', 'aa'])
        return pd.DataFrame(rows, index=multi_index, columns=col_order)

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
            close = difflib.get_close_matches(
                record_code, self.aaindex_json.keys(), n=3, cutoff=0.6
            )
            hint = f" Did you mean: {close}?" if close else ""
            raise ValueError(
                f"Record ({record_code}) not found in {self.__class__.__name__}.{hint}"
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
    def aaindex_json(self) -> Dict:
        """Parsed database dict, loaded lazily on first access."""
        if self._aaindex_json is None:
            self._load_data()
        return self._aaindex_json

    @aaindex_json.setter
    def aaindex_json(self, value: Dict) -> None:
        self._aaindex_json = value

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
