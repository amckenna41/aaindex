################################################################################
################             AAindex2 Module Tests             #################
################################################################################

import unittest
from aaindex import aaindex2, __version__

class AAIndex2_Tests(unittest.TestCase):
    """
    Test suite for testing the aaindex2 module in the aaindex Python software package.

    Test Cases
    ==========
    test_num_records:
        testing to check that the correct number of records are present in the AAi2 database.
    test_records:
        testing that correct record data and substitution matrix are returned for given
        accession numbers using the __getitem__ function.
    test_records_dot_notation:
        testing record data is accessible via dot notation using the Map class.
    test_get:
        testing the pairwise amino acid lookup, including symmetry and invalid input handling.
    test_search:
        testing search functionality returning records whose description contains keywords.
    test_record_codes:
        testing that correct record codes are present in the parsed AAi2 database.
    test_amino_acids:
        testing that only valid amino acid single-letter codes are returned.
    test_record_names:
        testing that record descriptions are returned correctly for all records.
    test_last_updated:
        testing the last updated date attribute matches the known database version.
    test_dunder_methods:
        testing __len__, __contains__, __iter__, and __repr__ dunder methods.
    """
    def test_num_records(self):
        """ Test Case to check the correct number of records are present in the AAi2 database.
        To date, 94 records are present in the database. """
#1.)
        self.assertEqual(aaindex2.num_records(), 94,
            f'Expected 94 records in AAi2, got {aaindex2.num_records()}.')

    def test_records(self):
        """ Test Case to check that correct record data and substitution matrix are
        returned for given accession numbers using the __getitem__ function. """
        index_code1 = 'ALTS910101'
        index_code2 = 'AZAE970101'
        index_code3 = 'alts910101'  # lowercase input should be normalised
        index_code4 = 'ABCDEFGH'
        index_code5 = ''
        index_code6 = '123456'
#1.)
        record = aaindex2[index_code1]  # ALTS910101

        #verify all expected keys are present in the record
        self.assertIn('description', list(record.keys()), 'description key not found in AAi2 record.')
        self.assertIn('references', list(record.keys()), 'references key not found in AAi2 record.')
        self.assertIn('notes', list(record.keys()), 'notes key not found in AAi2 record.')
        self.assertIn('pmid', list(record.keys()), 'pmid key not found in AAi2 record.')
        self.assertIn('correlation_coefficients', list(record.keys()), 'correlation_coefficients key not found in AAi2 record.')
        self.assertIn('matrix', list(record.keys()), 'matrix key not found in AAi2 record.')
        self.assertIn('row_order', list(record.keys()), 'row_order key not found in AAi2 record.')
        self.assertIn('col_order', list(record.keys()), 'col_order key not found in AAi2 record.')
        self.assertIn('is_symmetric', list(record.keys()), 'is_symmetric key not found in AAi2 record.')

        self.assertEqual(record['description'], 'The PAM-120 matrix (Altschul, 1991)',
            f"Unexpected description for ALTS910101, got {record['description']}.")
        self.assertEqual(record['pmid'], '1713145 2051488',
            f"Unexpected pmid for ALTS910101, got {record['pmid']}.")
        self.assertEqual(record['notes'], '',
            f"Unexpected notes for ALTS910101, got {record['notes']}.")

        #verify matrix is a nested dict keyed by amino acid
        self.assertIsInstance(record['matrix'], dict,
            f"matrix should be a dict, got {type(record['matrix'])}.")
        self.assertIn('A', record['matrix'],
            'Expected amino acid A as a key in the matrix dict.')
        self.assertIsInstance(record['matrix']['A'], dict,
            'Inner matrix values should be dicts keyed by amino acid.')

        #verify row_order and col_order are lists of 20 amino acids
        self.assertIsInstance(record['row_order'], list,
            f"row_order should be a list, got {type(record['row_order'])}.")
        self.assertEqual(len(record['row_order']), 20,
            f"row_order should contain 20 amino acids, got {len(record['row_order'])}.")
        self.assertEqual(len(record['col_order']), 20,
            f"col_order should contain 20 amino acids, got {len(record['col_order'])}.")
#2.)
        record2 = aaindex2[index_code2]  # AZAE970101

        self.assertTrue(record2['description'].startswith('The single residue substitution matrix'),
            f"Unexpected description for AZAE970101, got {record2['description']}.")
        self.assertEqual(record2['pmid'], '9488136',
            f"Unexpected pmid for AZAE970101, got {record2['pmid']}.")
        self.assertIsInstance(record2['matrix'], dict,
            f"matrix should be a dict, got {type(record2['matrix'])}.")
#3.)
        #lowercase record code should produce the same result as uppercase
        record3 = aaindex2[index_code3]
        self.assertEqual(record3['description'], record['description'],
            'Lowercase record code should return the same record as uppercase.')
#4.)
        #invalid record codes should raise ValueError
        with self.assertRaises(ValueError):
            aaindex2[index_code4]
        with self.assertRaises(ValueError):
            aaindex2[index_code5]
        with self.assertRaises(ValueError):
            aaindex2[index_code6]

    def test_records_dot_notation(self):
        """ Test Case to check record fields are accessible via dot notation
        using the Map class. """
        index_code1 = 'ALTS910101'
        index_code2 = 'AZAE970101'
#1.)
        record = aaindex2[index_code1]

        self.assertEqual(record.description, 'The PAM-120 matrix (Altschul, 1991)',
            f'Unexpected description via dot notation, got {record.description}.')
        self.assertEqual(record.pmid, '1713145 2051488',
            f'Unexpected pmid via dot notation, got {record.pmid}.')
        self.assertEqual(record.notes, '',
            f'Unexpected notes via dot notation, got {record.notes}.')
        self.assertTrue(record.references.startswith('Altschul, S.F.'),
            f'Unexpected references via dot notation, got {record.references}.')
        self.assertIsInstance(record.matrix, dict,
            'matrix via dot notation should be a dict.')
        self.assertIsInstance(record.row_order, list,
            'row_order via dot notation should be a list.')
#2.)
        record2 = aaindex2[index_code2]

        self.assertTrue(record2.description.startswith('The single residue substitution matrix'),
            f'Unexpected description for AZAE970101, got {record2.description}.')
        self.assertEqual(record2.pmid, '9488136',
            f'Unexpected pmid via dot notation, got {record2.pmid}.')

    def test_get(self):
        """ Test Case for the get() method: pairwise lookup, symmetry,
        and invalid input handling. """
        index_code1 = 'ALTS910101'
        index_code2 = 'AZAE970101'
#1.)
        #test known scalar values from the substitution matrix
        self.assertEqual(aaindex2.get(index_code1, 'A', 'A'), 3.0,
            f"Expected A,A = 3.0 for ALTS910101, got {aaindex2.get(index_code1, 'A', 'A')}.")
        self.assertEqual(aaindex2.get(index_code1, 'A', 'R'), -3.0,
            f"Expected A,R = -3.0 for ALTS910101, got {aaindex2.get(index_code1, 'A', 'R')}.")
        self.assertEqual(aaindex2.get(index_code2, 'L', 'L'), 9.0,
            f"Expected L,L = 9.0 for AZAE970101, got {aaindex2.get(index_code2, 'L', 'L')}.")
        self.assertEqual(aaindex2.get(index_code2, 'V', 'I'), 9.0,
            f"Expected V,I = 9.0 for AZAE970101, got {aaindex2.get(index_code2, 'V', 'I')}.")
#2.)
        #symmetry: (aa1, aa2) must equal (aa2, aa1)
        self.assertEqual(aaindex2.get(index_code1, 'A', 'R'), aaindex2.get(index_code1, 'R', 'A'),
            'get() should be symmetric: (A,R) should equal (R,A) for ALTS910101.')
        self.assertEqual(aaindex2.get(index_code2, 'V', 'I'), aaindex2.get(index_code2, 'I', 'V'),
            'get() should be symmetric: (V,I) should equal (I,V) for AZAE970101.')
#3.)
        #lowercase amino acid codes should be normalised to uppercase
        self.assertEqual(aaindex2.get(index_code1, 'a', 'r'), aaindex2.get(index_code1, 'A', 'R'),
            'get() should accept lowercase amino acid codes.')
#4.)
        #invalid record code should raise ValueError
        with self.assertRaises(ValueError):
            aaindex2.get('BLAH999999', 'A', 'R')
#5.)
        #non-string amino acid input should raise TypeError
        with self.assertRaises(TypeError):
            aaindex2.get(index_code1, 1, 'R')

    def test_search(self):
        """ Test Case for search(), returning records whose description or other text
        fields contain the given keyword(s). """
        description1 = 'substitution'
        description2 = 'blahblahblah'
        description3 = 'not a description'
        description4 = 1234
        description5 = True
#1.)
        #search for 'substitution' - at least the original 39 records match via description;
        #searching all fields (references, notes, pmid, accession) may return more
        search1 = aaindex2.search(description1)
        self.assertGreaterEqual(len(search1), 39,
            f'Expected at least 39 records for substitution search, got {len(search1)}.')
        #records matched via description still have the keyword in their description
        desc_matched = {k: v for k, v in search1.items() if description1.lower() in v['description'].lower()}
        self.assertEqual(len(desc_matched), 39,
            f'Expected 39 records to match via description, got {len(desc_matched)}.')
#2.)
        #nonsense search terms should return empty dict
        search2 = aaindex2.search(description2)
        self.assertEqual(len(search2), 0,
            f'Expected 0 records for nonsense search, got {len(search2)}.')
#3.)
        search3 = aaindex2.search(description3)
        self.assertEqual(len(search3), 0,
            f'Expected 0 records for nonsense search, got {len(search3)}.')
#4.)
        #non-string/list input should raise TypeError
        with self.assertRaises(TypeError):
            aaindex2.search(description4)
        with self.assertRaises(TypeError):
            aaindex2.search(description5)
#5.) search by PMID should return the matching record
        pmid_search = aaindex2.search('1713145')
        self.assertIn('ALTS910101', pmid_search,
            'search() should find ALTS910101 when searching by its PMID 1713145.')
#6.) search by exact accession code should return the record itself
        accession_search = aaindex2.search('ALTS910101')
        self.assertIn('ALTS910101', accession_search,
            'search() should find ALTS910101 when searching by its accession code.')

    def test_record_codes(self):
        """ Test Case to check that correct record codes are present in the parsed
        AAi2 database, and that closely named erroneous codes are not. """
        index_code1 = 'ALTS910101'
        index_code2 = 'AZAE970101'
        index_code3 = 'WEIL970101'
        index_code4 = 'TUDE900101'
#1.)
        #known record codes should be present
        self.assertIn(index_code1, aaindex2.record_codes(),
            f'Expected {index_code1} in AAi2 record codes.')
        self.assertIn(index_code2, aaindex2.record_codes(),
            f'Expected {index_code2} in AAi2 record codes.')
        self.assertIn(index_code3, aaindex2.record_codes(),
            f'Expected {index_code3} in AAi2 record codes.')
        self.assertIn(index_code4, aaindex2.record_codes(),
            f'Expected {index_code4} in AAi2 record codes.')
#2.)
        #erroneous codes should not be present
        self.assertNotIn('BLAH999999', aaindex2.record_codes(),
            'Expected BLAH999999 to not be in AAi2 record codes.')
        self.assertNotIn('ALTS910199', aaindex2.record_codes(),
            'Expected ALTS910199 to not be in AAi2 record codes.')

    def test_amino_acids(self):
        """ Test Case to check that only valid amino acid single-letter codes are returned. """
        valid_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
            'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        amino_acids = aaindex2.amino_acids()
#1.)
        self.assertEqual(len(amino_acids), 20,
            f'Expected 20 amino acids, got {len(amino_acids)}.')
        for aa in amino_acids:
            self.assertIn(aa, valid_amino_acids,
                f'Amino acid {aa} not in list of valid amino acids.')

    def test_record_names(self):
        """ Test Case to check that all record descriptions are returned correctly. """
        record_names = aaindex2.record_names()
#1.)
        #number of names must equal number of records
        self.assertEqual(len(record_names), aaindex2.num_records(),
            f'Number of record names should equal num_records(), got {len(record_names)}.')
#2.)
        #known descriptions should be present
        self.assertIn('The PAM-120 matrix (Altschul, 1991)', record_names,
            'Expected ALTS910101 description in record names.')
        self.assertTrue(any('single residue substitution matrix' in n for n in record_names),
            'Expected AZAE970101 description in record names.')

    def test_last_updated(self):
        """ Testing the last updated class attribute matches the known database version. """
        self.assertEqual(aaindex2.last_updated, "February 13, 2017",
            f'Last updated value does not match expected, got {aaindex2.last_updated}.')

    def test_dunder_methods(self):
        """ Test Case for __len__, __contains__, __iter__, and __repr__ dunder methods. """
#1.) __len__
        self.assertEqual(len(aaindex2), 94,
            f'Expected __len__ to return 94, got {len(aaindex2)}.')
#2.) __contains__
        self.assertIn('ALTS910101', aaindex2,
            'Expected ALTS910101 in aaindex2 via __contains__.')
        self.assertNotIn('BLAH999999', aaindex2,
            'Expected BLAH999999 to not be in aaindex2 via __contains__.')
#3.) __iter__
        codes_from_iter = list(aaindex2)
        self.assertEqual(len(codes_from_iter), 94,
            f'Expected 94 codes from __iter__, got {len(codes_from_iter)}.')
        self.assertIn('ALTS910101', codes_from_iter,
            'Expected ALTS910101 in __iter__ output.')
#4.) __repr__
        self.assertEqual(repr(aaindex2), "AAIndex2(records=94, last_updated='February 13, 2017')",
            f'Unexpected __repr__ output: {repr(aaindex2)}.')

    def test_values(self):
        """ Test Case for values() which returns the full matrix dict for a record. """
#1.)
        matrix = aaindex2.values('ALTS910101')
        self.assertIsInstance(matrix, dict,
            f'values() should return a dict, got {type(matrix)}.')
        self.assertIn('A', matrix,
            'Expected amino acid A as key in values() result.')
        self.assertEqual(matrix['A']['A'], 3.0,
            'Expected A,A = 3.0 in values() result.')
#2.)
        #invalid record code should raise ValueError
        with self.assertRaises(ValueError):
            aaindex2.values('BLAH999999')

    def test_symmetry(self):
        """ Explicit symmetry test: for symmetric records, (aa1, aa2) == (aa2, aa1);
        for non-symmetric records, is_symmetric flag is False. """
        amino_acids = aaindex2.amino_acids()
        #test symmetry on two representative symmetric records
        for code in ['ALTS910101', 'AZAE970101']:
            record = aaindex2[code]
            self.assertTrue(record['is_symmetric'],
                f'Expected is_symmetric=True for {code}.')
            for aa1 in amino_acids:
                for aa2 in amino_acids:
                    val_forward = aaindex2.get(code, aa1, aa2)
                    val_reverse = aaindex2.get(code, aa2, aa1)
                    self.assertEqual(val_forward, val_reverse,
                        f'Symmetry violation for {code}: ({aa1},{aa2})={val_forward} != ({aa2},{aa1})={val_reverse}')

    def test_non_symmetric_matrices(self):
        """ Test that non-symmetric records store directional values correctly
        and expose is_symmetric=False. """
        #LINK010101 is a known non-symmetric substitution matrix
        record = aaindex2['LINK010101']
#1.) is_symmetric flag must be False
        self.assertFalse(record['is_symmetric'],
            'LINK010101 should be flagged as non-symmetric.')
#2.) A->R and R->A should differ in the raw data (0.030 vs 0.034)
        a_to_r = aaindex2.get('LINK010101', 'A', 'R')
        r_to_a = aaindex2.get('LINK010101', 'R', 'A')
        self.assertIsNotNone(a_to_r, 'A->R should not be None for LINK010101.')
        self.assertIsNotNone(r_to_a, 'R->A should not be None for LINK010101.')
        self.assertNotEqual(a_to_r, r_to_a,
            f'Expected A->R ({a_to_r}) != R->A ({r_to_a}) for non-symmetric LINK010101.')
        self.assertAlmostEqual(a_to_r, 0.030, places=3,
            msg=f'Expected A->R = 0.030 for LINK010101, got {a_to_r}.')
        self.assertAlmostEqual(r_to_a, 0.034, places=3,
            msg=f'Expected R->A = 0.034 for LINK010101, got {r_to_a}.')
#3.) symmetric records should retain is_symmetric=True
        sym_record = aaindex2['ALTS910101']
        self.assertTrue(sym_record['is_symmetric'],
            'ALTS910101 should be flagged as symmetric.')

    def test_edge_cases(self):
        """ Test edge cases: invalid AA in get(), empty search, non-string input. """
#1.) non-canonical amino acid "Z" should return None from get()
        result = aaindex2.get('ALTS910101', 'Z', 'A')
        self.assertIsNone(result,
            'get() with non-canonical amino acid Z should return None.')
#2.) empty string search should match all records
        empty_search = aaindex2.search("")
        self.assertEqual(len(empty_search), aaindex2.num_records(),
            'Empty string search should match all records.')
        self.assertIn('ALTS910101', empty_search,
            'Expected ALTS910101 in empty string search results.')
#3.) non-string record code for __getitem__ should raise TypeError
        with self.assertRaises(TypeError):
            aaindex2[12345]
        with self.assertRaises(TypeError):
            aaindex2[None]
#4.) whitespace-padded record code should resolve
        record = aaindex2['  ALTS910101  ']
        self.assertEqual(record.description, 'The PAM-120 matrix (Altschul, 1991)',
            'Whitespace-padded record code should resolve correctly.')

    def test_to_dict(self):
        """ Test Case for to_dict() exporting records as plain Python dicts. """
        code = 'ALTS910101'
#1.) single record export
        d = aaindex2.to_dict(code)
        self.assertIsInstance(d, dict,
            f'to_dict() should return a dict, got {type(d)}.')
        self.assertIn(code, d,
            f'Expected {code} as key in to_dict() result.')
        self.assertIn('matrix', d[code],
            'Expected matrix field in exported record.')
#2.) full database export
        full = aaindex2.to_dict()
        self.assertEqual(len(full), aaindex2.num_records(),
            f'Full to_dict() should contain all {aaindex2.num_records()} records.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex2.to_dict('BLAH999999')

    def test_to_json(self):
        """ Test Case for to_json() serialising records to a JSON string. """
        import json
        code = 'ALTS910101'
#1.) single record serialisation returns valid JSON
        j = aaindex2.to_json(code)
        self.assertIsInstance(j, str,
            f'to_json() should return a str, got {type(j)}.')
        parsed = json.loads(j)
        self.assertIn(code, parsed,
            f'Expected {code} as key in parsed to_json() output.')
#2.) full database serialisation is valid JSON containing all records
        full_parsed = json.loads(aaindex2.to_json())
        self.assertEqual(len(full_parsed), aaindex2.num_records(),
            f'Full to_json() should contain all {aaindex2.num_records()} records.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex2.to_json('BLAH999999')

    def test_to_dataframe(self):
        """ Test Case for to_dataframe(). Skipped if pandas is not installed. """
        try:
            import pandas as pd
        except ImportError:
            self.skipTest('pandas not installed — skipping test_to_dataframe.')
        code = 'ALTS910101'
#1.) single record returns a DataFrame with amino acids as both index and columns
        df = aaindex2.to_dataframe(code)
        self.assertIsInstance(df, pd.DataFrame,
            f'to_dataframe() should return a DataFrame, got {type(df)}.')
        self.assertGreater(len(df), 0,
            'Single-record to_dataframe() should have at least 1 row.')
        for aa in ['A', 'R']:
            self.assertIn(aa, df.index, f'Expected row {aa} in to_dataframe() result.')
            self.assertIn(aa, df.columns, f'Expected column {aa} in to_dataframe() result.')
#2.) full database returns a MultiIndex DataFrame
        full_df = aaindex2.to_dataframe()
        self.assertIsInstance(full_df.index, pd.MultiIndex,
            'Full to_dataframe() should return a MultiIndex DataFrame.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex2.to_dataframe('BLAH999999')

    def test_search_fuzzy(self):
        """ Test Case for search_fuzzy() fuzzy search with ranked results. """
#1.) returns a dict
        results = aaindex2.search_fuzzy('PAM')
        self.assertIsInstance(results, dict,
            f'search_fuzzy() should return a dict, got {type(results)}.')
#2.) non-empty results for a known term
        self.assertGreater(len(results), 0,
            'search_fuzzy() should return at least one result for "PAM".')
#3.) respects n limit
        results_n = aaindex2.search_fuzzy('PAM', n=2)
        self.assertLessEqual(len(results_n), 2,
            'search_fuzzy() should return at most n results.')
#4.) TypeError for non-str query
        with self.assertRaises(TypeError):
            aaindex2.search_fuzzy(123)
#5.) ValueError for n < 1
        with self.assertRaises(ValueError):
            aaindex2.search_fuzzy('PAM', n=0)
#6.) ValueError for cutoff outside [0, 1]
        with self.assertRaises(ValueError):
            aaindex2.search_fuzzy('PAM', cutoff=2.0)

    def test_richer_error_messages(self):
        """ Test that __getitem__ raises ValueError with an informative message. """
#1.) completely unknown code raises ValueError
        with self.assertRaises(ValueError):
            aaindex2['BLAH999999']
#2.) error message mentions the queried code
        try:
            aaindex2['BLAH999999']
        except ValueError as e:
            self.assertIn('BLAH999999', str(e),
                f'Error message should mention the queried code, got: {str(e)}.')

    def test_parse_aaindex_deprecated(self):
        """ Test that calling the public parse_aaindex() shim emits a DeprecationWarning. """
        import warnings as _warnings
        with _warnings.catch_warnings(record=True) as w:
            _warnings.simplefilter('always')
            aaindex2.parse_aaindex()
        self.assertTrue(
            any(issubclass(warning.category, DeprecationWarning) for warning in w),
            'parse_aaindex() should emit a DeprecationWarning.')

    def test_instance_version(self):
        """ Test that AAIndex2 instance exposes a version attribute matching __version__. """
#1.) version attribute is a string
        self.assertIsInstance(aaindex2.version, str,
            f'aaindex2.version should be a str, got {type(aaindex2.version)}.')
#2.) version matches the package __version__
        self.assertEqual(aaindex2.version, __version__,
            f'aaindex2.version should match __version__, got {aaindex2.version}.')

    def test_property_setter_validation(self):
        """ Test that data_dir and aaindex_filename setters reject path traversal values. """
        from aaindex._aaindex_matrix import _AAIndexMatrix
        fresh = _AAIndexMatrix.__new__(_AAIndexMatrix)
#1.) relative path in data_dir should raise ValueError
        with self.assertRaises(ValueError):
            fresh.data_dir = '../../etc'
#2.) path with separator in aaindex_filename should raise ValueError
        with self.assertRaises(ValueError):
            fresh.aaindex_filename = '../aaindex2'
#3.) simple names are accepted
        fresh.data_dir = 'data'
        self.assertEqual(fresh.data_dir, 'data')
        fresh.aaindex_filename = 'aaindex2'
        self.assertEqual(fresh.aaindex_filename, 'aaindex2')

    def test_correlation_coefficients_are_floats(self):
        """ Test that correlation_coefficients values are stored as floats, not strings. """
        #find a record that has at least one correlation coefficient
        for code in aaindex2.record_codes():
            record = aaindex2[code]
            cc = record['correlation_coefficients']
            if cc:
                for k, val in cc.items():
                    self.assertIsInstance(val, float,
                        f'correlation_coefficients[{k!r}] in {code} should be float, got {type(val)}.')
                break

if __name__ == '__main__':
    #run all unit tests
    unittest.main(verbosity=2)
