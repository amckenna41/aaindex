################################################################################
################             AAindex3 Module Tests             #################
################################################################################

import unittest
from aaindex import aaindex3, __version__

class AAIndex3_Tests(unittest.TestCase):
    """
    Test suite for testing the aaindex3 module in the aaindex Python software package.

    Test Cases
    ==========
    test_num_records:
        testing to check that the correct number of records are present in the AAi3 database.
    test_records:
        testing that correct record data and contact potential matrix are returned for given
        accession numbers using the __getitem__ function.
    test_records_dot_notation:
        testing record data is accessible via dot notation using the Map class.
    test_get:
        testing the pairwise amino acid lookup, including symmetry, NA handling,
        and invalid input handling.
    test_search:
        testing search functionality returning records whose description contains keywords.
    test_record_codes:
        testing that correct record codes are present in the parsed AAi3 database.
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
        """ Test Case to check the correct number of records are present in the AAi3 database.
        To date, 47 records are present in the database. """
#1.)
        self.assertEqual(aaindex3.num_records(), 47,
            f'Expected 47 records in AAi3, got {aaindex3.num_records()}.')

    def test_records(self):
        """ Test Case to check that correct record data and contact potential matrix are
        returned for given accession numbers using the __getitem__ function. """
        index_code1 = 'TANS760101'
        index_code2 = 'BETM990101'
        index_code3 = 'tans760101'  # lowercase input should be normalised
        index_code4 = 'ABCDEFGH'
        index_code5 = ''
        index_code6 = '123456'
#1.)
        record = aaindex3[index_code1]  # TANS760101

        #verify all expected keys are present in the record
        self.assertIn('description', list(record.keys()), 'description key not found in AAi3 record.')
        self.assertIn('references', list(record.keys()), 'references key not found in AAi3 record.')
        self.assertIn('notes', list(record.keys()), 'notes key not found in AAi3 record.')
        self.assertIn('pmid', list(record.keys()), 'pmid key not found in AAi3 record.')
        self.assertIn('correlation_coefficients', list(record.keys()), 'correlation_coefficients key not found in AAi3 record.')
        self.assertIn('matrix', list(record.keys()), 'matrix key not found in AAi3 record.')
        self.assertIn('row_order', list(record.keys()), 'row_order key not found in AAi3 record.')
        self.assertIn('col_order', list(record.keys()), 'col_order key not found in AAi3 record.')

        self.assertEqual(record['description'], 'Statistical contact potential derived from 25 x-ray protein structures',
            f"Unexpected description for TANS760101, got {record['description']}.")
        self.assertEqual(record['pmid'], '1004017',
            f"Unexpected pmid for TANS760101, got {record['pmid']}.")

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
        self.assertIn('is_symmetric', list(record.keys()), 'is_symmetric key not found in AAi3 record.')
#2.)
        record2 = aaindex3[index_code2]  # BETM990101

        self.assertEqual(record2['description'], 'Modified version of the Miyazawa-Jernigan transfer energy',
            f"Unexpected description for BETM990101, got {record2['description']}.")
        self.assertEqual(record2['pmid'], '10048329',
            f"Unexpected pmid for BETM990101, got {record2['pmid']}.")
        self.assertIsInstance(record2['matrix'], dict,
            f"matrix should be a dict for BETM990101, got {type(record2['matrix'])}.")
#3.)
        #lowercase record code should produce the same result as uppercase
        record3 = aaindex3[index_code3]
        self.assertEqual(record3['description'], record['description'],
            'Lowercase record code should return the same record as uppercase.')
#4.)
        #invalid record codes should raise ValueError
        with self.assertRaises(ValueError):
            aaindex3[index_code4]
        with self.assertRaises(ValueError):
            aaindex3[index_code5]
        with self.assertRaises(ValueError):
            aaindex3[index_code6]

    def test_records_dot_notation(self):
        """ Test Case to check record fields are accessible via dot notation
        using the Map class. """
        index_code1 = 'TANS760101'
        index_code2 = 'BETM990101'
#1.)
        record = aaindex3[index_code1]

        self.assertEqual(record.description, 'Statistical contact potential derived from 25 x-ray protein structures',
            f'Unexpected description via dot notation, got {record.description}.')
        self.assertEqual(record.pmid, '1004017',
            f'Unexpected pmid via dot notation, got {record.pmid}.')
        self.assertTrue(record.references.startswith('Tanaka, S. and Scheraga, H.A.'),
            f'Unexpected references via dot notation, got {record.references}.')
        self.assertIsInstance(record.matrix, dict,
            'matrix via dot notation should be a dict.')
        self.assertIsInstance(record.row_order, list,
            'row_order via dot notation should be a list.')
#2.)
        record2 = aaindex3[index_code2]

        self.assertEqual(record2.description, 'Modified version of the Miyazawa-Jernigan transfer energy',
            f'Unexpected description for BETM990101, got {record2.description}.')
        self.assertEqual(record2.pmid, '10048329',
            f'Unexpected pmid via dot notation, got {record2.pmid}.')

    def test_get(self):
        """ Test Case for the get() method: pairwise lookup, symmetry, NA handling,
        and invalid input handling. """
        index_code1 = 'TANS760101'
        index_code2 = 'BETM990101'
        index_code3 = 'ROBB790102'
#1.)
        #test known scalar values from the contact potential matrix
        self.assertEqual(aaindex3.get(index_code1, 'A', 'A'), -2.6,
            f"Expected A,A = -2.6 for TANS760101, got {aaindex3.get(index_code1, 'A', 'A')}.")
        self.assertEqual(aaindex3.get(index_code1, 'A', 'R'), -3.4,
            f"Expected A,R = -3.4 for TANS760101, got {aaindex3.get(index_code1, 'A', 'R')}.")
        self.assertEqual(aaindex3.get(index_code2, 'A', 'A'), -0.2,
            f"Expected A,A = -0.2 for BETM990101, got {aaindex3.get(index_code2, 'A', 'A')}.")
        self.assertEqual(aaindex3.get(index_code2, 'L', 'V'), -0.8,
            f"Expected L,V = -0.8 for BETM990101, got {aaindex3.get(index_code2, 'L', 'V')}.")
#2.)
        #symmetry: (aa1, aa2) must equal (aa2, aa1)
        self.assertEqual(aaindex3.get(index_code1, 'A', 'R'), aaindex3.get(index_code1, 'R', 'A'),
            'get() should be symmetric: (A,R) should equal (R,A) for TANS760101.')
        self.assertEqual(aaindex3.get(index_code2, 'L', 'V'), aaindex3.get(index_code2, 'V', 'L'),
            'get() should be symmetric: (L,V) should equal (V,L) for BETM990101.')
#3.)
        #NA values should be returned as None
        self.assertIsNone(aaindex3.get(index_code3, 'G', 'A'),
            'Expected None for G,A in ROBB790102 (Gly is not available in this record).')
        self.assertIsNone(aaindex3.get(index_code3, 'A', 'G'),
            'Expected None for A,G in ROBB790102 (symmetric NA lookup).')
#4.)
        #lowercase amino acid codes should be normalised to uppercase
        self.assertEqual(aaindex3.get(index_code1, 'a', 'r'), aaindex3.get(index_code1, 'A', 'R'),
            'get() should accept lowercase amino acid codes.')
#5.)
        #invalid record code should raise ValueError
        with self.assertRaises(ValueError):
            aaindex3.get('BLAH999999', 'A', 'R')
#6.)
        #non-string amino acid input should raise TypeError
        with self.assertRaises(TypeError):
            aaindex3.get(index_code1, 1, 'R')

    def test_search(self):
        """ Test Case for search(), returning records whose description or other text
        fields contain the given keyword(s). """
        description1 = 'contact potential'
        description2 = 'blahblahblah'
        description3 = 'not a description'
        description4 = 1234
        description5 = True
#1.)
        #search for 'contact potential' - at least 3 records match via description
        search1 = aaindex3.search(description1)
        self.assertGreaterEqual(len(search1), 3,
            f'Expected at least 3 records for contact potential search, got {len(search1)}.')
        desc_matched = {k: v for k, v in search1.items() if description1.lower() in v['description'].lower()}
        self.assertEqual(len(desc_matched), 3,
            f'Expected 3 records to match via description, got {len(desc_matched)}.')
#2.)
        #nonsense search terms should return empty dict
        search2 = aaindex3.search(description2)
        self.assertEqual(len(search2), 0,
            f'Expected 0 records for nonsense search, got {len(search2)}.')
#3.)
        search3 = aaindex3.search(description3)
        self.assertEqual(len(search3), 0,
            f'Expected 0 records for nonsense search, got {len(search3)}.')
#4.)
        #non-string/list input should raise TypeError
        with self.assertRaises(TypeError):
            aaindex3.search(description4)
        with self.assertRaises(TypeError):
            aaindex3.search(description5)
#5.) search by PMID should return the matching record
        pmid_search = aaindex3.search('1004017')
        self.assertIn('TANS760101', pmid_search,
            'search() should find TANS760101 when searching by its PMID 1004017.')
#6.) search by exact accession code should return the record itself
        accession_search = aaindex3.search('TANS760101')
        self.assertIn('TANS760101', accession_search,
            'search() should find TANS760101 when searching by its accession code.')

    def test_record_codes(self):
        """ Test Case to check that correct record codes are present in the parsed
        AAi3 database, and that erroneous codes are not. """
        index_code1 = 'TANS760101'
        index_code2 = 'BETM990101'
        index_code3 = 'ROBB790102'
        index_code4 = 'ZHAC000106'
#1.)
        #known record codes should be present
        self.assertIn(index_code1, aaindex3.record_codes(),
            f'Expected {index_code1} in AAi3 record codes.')
        self.assertIn(index_code2, aaindex3.record_codes(),
            f'Expected {index_code2} in AAi3 record codes.')
        self.assertIn(index_code3, aaindex3.record_codes(),
            f'Expected {index_code3} in AAi3 record codes.')
        self.assertIn(index_code4, aaindex3.record_codes(),
            f'Expected {index_code4} in AAi3 record codes.')
#2.)
        #erroneous codes should not be present
        self.assertNotIn('BLAH999999', aaindex3.record_codes(),
            'Expected BLAH999999 to not be in AAi3 record codes.')
        self.assertNotIn('TANS760199', aaindex3.record_codes(),
            'Expected TANS760199 to not be in AAi3 record codes.')

    def test_amino_acids(self):
        """ Test Case to check that only valid amino acid single-letter codes are returned. """
        valid_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
            'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        amino_acids = aaindex3.amino_acids()
#1.)
        self.assertEqual(len(amino_acids), 20,
            f'Expected 20 amino acids, got {len(amino_acids)}.')
        for aa in amino_acids:
            self.assertIn(aa, valid_amino_acids,
                f'Amino acid {aa} not in list of valid amino acids.')

    def test_record_names(self):
        """ Test Case to check that all record descriptions are returned correctly. """
        record_names = aaindex3.record_names()
#1.)
        #number of names must equal number of records
        self.assertEqual(len(record_names), aaindex3.num_records(),
            f'Number of record names should equal num_records(), got {len(record_names)}.')
#2.)
        #known descriptions should be present in the list
        self.assertIn('Statistical contact potential derived from 25 x-ray protein structures', record_names,
            'Expected TANS760101 description in record names.')
        self.assertIn('Modified version of the Miyazawa-Jernigan transfer energy', record_names,
            'Expected BETM990101 description in record names.')

    def test_last_updated(self):
        """ Testing the last updated class attribute matches the known database version. """
        self.assertEqual(aaindex3.last_updated, "February 13, 2017",
            f'Last updated value does not match expected, got {aaindex3.last_updated}.')

    def test_dunder_methods(self):
        """ Test Case for __len__, __contains__, __iter__, and __repr__ dunder methods. """
#1.) __len__
        self.assertEqual(len(aaindex3), 47,
            f'Expected __len__ to return 47, got {len(aaindex3)}.')
#2.) __contains__
        self.assertIn('TANS760101', aaindex3,
            'Expected TANS760101 in aaindex3 via __contains__.')
        self.assertNotIn('BLAH999999', aaindex3,
            'Expected BLAH999999 to not be in aaindex3 via __contains__.')
#3.) __iter__
        codes_from_iter = list(aaindex3)
        self.assertEqual(len(codes_from_iter), 47,
            f'Expected 47 codes from __iter__, got {len(codes_from_iter)}.')
        self.assertIn('TANS760101', codes_from_iter,
            'Expected TANS760101 in __iter__ output.')
#4.) __repr__
        self.assertEqual(repr(aaindex3), "AAIndex3(records=47, last_updated='February 13, 2017')",
            f'Unexpected __repr__ output: {repr(aaindex3)}.')

    def test_values(self):
        """ Test Case for values() which returns the full matrix dict for a record. """
#1.)
        matrix = aaindex3.values('TANS760101')
        self.assertIsInstance(matrix, dict,
            f'values() should return a dict, got {type(matrix)}.')
        self.assertIn('A', matrix,
            'Expected amino acid A as key in values() result.')
        self.assertEqual(matrix['A']['A'], -2.6,
            'Expected A,A = -2.6 in values() result.')
#2.)
        #invalid record code should raise ValueError
        with self.assertRaises(ValueError):
            aaindex3.values('BLAH999999')

    def test_symmetry(self):
        """ Explicit symmetry test: for symmetric records, (aa1, aa2) == (aa2, aa1);
        all AAindex3 records are lower-triangular (symmetric). """
        amino_acids = aaindex3.amino_acids()
        #test symmetry and is_symmetric flag on two representative records
        for code in ['TANS760101', 'BETM990101']:
            record = aaindex3[code]
            self.assertTrue(record['is_symmetric'],
                f'Expected is_symmetric=True for AAindex3 record {code}.')
            for aa1 in amino_acids:
                for aa2 in amino_acids:
                    val_forward = aaindex3.get(code, aa1, aa2)
                    val_reverse = aaindex3.get(code, aa2, aa1)
                    self.assertEqual(val_forward, val_reverse,
                        f'Symmetry violation for {code}: ({aa1},{aa2})={val_forward} != ({aa2},{aa1})={val_reverse}')

    def test_edge_cases(self):
        """ Test edge cases: invalid AA in get(), empty search, non-string input. """
#1.) non-canonical amino acid "Z" should return None from get()
        result = aaindex3.get('TANS760101', 'Z', 'A')
        self.assertIsNone(result,
            'get() with non-canonical amino acid Z should return None.')
#2.) empty string search should match all records
        empty_search = aaindex3.search("")
        self.assertEqual(len(empty_search), aaindex3.num_records(),
            'Empty string search should match all records.')
        self.assertIn('TANS760101', empty_search,
            'Expected TANS760101 in empty string search results.')
#3.) non-string record code for __getitem__ should raise TypeError
        with self.assertRaises(TypeError):
            aaindex3[12345]
        with self.assertRaises(TypeError):
            aaindex3[None]
#4.) whitespace-padded record code should resolve
        record = aaindex3['  TANS760101  ']
        self.assertEqual(record.description,
            'Statistical contact potential derived from 25 x-ray protein structures',
            'Whitespace-padded record code should resolve correctly.')

    def test_to_dict(self):
        """ Test Case for to_dict() exporting records as plain Python dicts. """
        code = 'TANS760101'
#1.) single record export
        d = aaindex3.to_dict(code)
        self.assertIsInstance(d, dict,
            f'to_dict() should return a dict, got {type(d)}.')
        self.assertIn(code, d,
            f'Expected {code} as key in to_dict() result.')
        self.assertIn('matrix', d[code],
            'Expected matrix field in exported record.')
#2.) full database export
        full = aaindex3.to_dict()
        self.assertEqual(len(full), aaindex3.num_records(),
            f'Full to_dict() should contain all {aaindex3.num_records()} records.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex3.to_dict('BLAH999999')

    def test_to_json(self):
        """ Test Case for to_json() serialising records to a JSON string. """
        import json
        code = 'TANS760101'
#1.) single record serialisation returns valid JSON
        j = aaindex3.to_json(code)
        self.assertIsInstance(j, str,
            f'to_json() should return a str, got {type(j)}.')
        parsed = json.loads(j)
        self.assertIn(code, parsed,
            f'Expected {code} as key in parsed to_json() output.')
#2.) full database serialisation is valid JSON containing all records
        full_parsed = json.loads(aaindex3.to_json())
        self.assertEqual(len(full_parsed), aaindex3.num_records(),
            f'Full to_json() should contain all {aaindex3.num_records()} records.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex3.to_json('BLAH999999')

    def test_to_dataframe(self):
        """ Test Case for to_dataframe(). Skipped if pandas is not installed. """
        try:
            import pandas as pd
        except ImportError:
            self.skipTest('pandas not installed — skipping test_to_dataframe.')
        code = 'TANS760101'
#1.) single record returns a DataFrame with amino acids as both index and columns
        df = aaindex3.to_dataframe(code)
        self.assertIsInstance(df, pd.DataFrame,
            f'to_dataframe() should return a DataFrame, got {type(df)}.')
        self.assertGreater(len(df), 0,
            'Single-record to_dataframe() should have at least 1 row.')
        for aa in ['A', 'C']:
            self.assertIn(aa, df.index, f'Expected row {aa} in to_dataframe() result.')
            self.assertIn(aa, df.columns, f'Expected column {aa} in to_dataframe() result.')
#2.) full database returns a MultiIndex DataFrame
        full_df = aaindex3.to_dataframe()
        self.assertIsInstance(full_df.index, pd.MultiIndex,
            'Full to_dataframe() should return a MultiIndex DataFrame.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex3.to_dataframe('BLAH999999')

    def test_search_fuzzy(self):
        """ Test Case for search_fuzzy() fuzzy search with ranked results. """
#1.) returns a dict
        results = aaindex3.search_fuzzy('contact')
        self.assertIsInstance(results, dict,
            f'search_fuzzy() should return a dict, got {type(results)}.')
#2.) non-empty results for a known term
        self.assertGreater(len(results), 0,
            'search_fuzzy() should return at least one result for "contact".')
#3.) respects n limit
        results_n = aaindex3.search_fuzzy('contact', n=2)
        self.assertLessEqual(len(results_n), 2,
            'search_fuzzy() should return at most n results.')
#4.) TypeError for non-str query
        with self.assertRaises(TypeError):
            aaindex3.search_fuzzy(123)
#5.) ValueError for n < 1
        with self.assertRaises(ValueError):
            aaindex3.search_fuzzy('contact', n=0)
#6.) ValueError for cutoff outside [0, 1]
        with self.assertRaises(ValueError):
            aaindex3.search_fuzzy('contact', cutoff=-0.1)

    def test_richer_error_messages(self):
        """ Test that __getitem__ includes a did-you-mean hint on invalid accession. """
#1.) completely unknown code raises ValueError
        with self.assertRaises(ValueError):
            aaindex3['BLAH999999']

    def test_parse_aaindex_deprecated(self):
        """ Test that calling the public parse_aaindex() shim emits a DeprecationWarning. """
        import warnings as _warnings
        with _warnings.catch_warnings(record=True) as w:
            _warnings.simplefilter('always')
            aaindex3.parse_aaindex()
        self.assertTrue(
            any(issubclass(warning.category, DeprecationWarning) for warning in w),
            'parse_aaindex() should emit a DeprecationWarning.')

    def test_search_sorted(self):
        """ Test that search() returns results sorted alphabetically by accession code. """
        results = aaindex3.search('contact potential')
        keys = list(results.keys())
        self.assertEqual(keys, sorted(keys),
            'search() results should be sorted alphabetically by accession code.')

    def test_plot_heatmap(self):
        """ Test Case for plot_heatmap(), which plots a 20x20 matrix as a heatmap.
        Skipped when matplotlib is not installed. """
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            self.skipTest('matplotlib not installed — skipping test_plot_heatmap.')
        code = 'TANS760101'
#1.) returns a matplotlib Axes object
        ax = aaindex3.plot_heatmap(code)
        self.assertIsNotNone(ax,
            'plot_heatmap() should return an Axes object.')
        self.assertEqual(ax.__class__.__name__, 'Axes',
            f'plot_heatmap() should return Axes, got {type(ax)}.')
#2.) title contains the record description
        record = aaindex3[code]
        self.assertTrue(ax.get_title().startswith(record['description'][:20]),
            'plot_heatmap() title should start with the record description.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex3.plot_heatmap('BLAH999999')
        plt.close('all')

    def test_instance_version(self):
        """ Test that AAIndex3 instance exposes a version attribute matching __version__. """
#1.) version attribute is a string
        self.assertIsInstance(aaindex3.version, str,
            f'aaindex3.version should be a str, got {type(aaindex3.version)}.')
#2.) version matches the package __version__
        self.assertEqual(aaindex3.version, __version__,
            f'aaindex3.version should match __version__, got {aaindex3.version}.')

if __name__ == '__main__':
    #run all unit tests
    unittest.main(verbosity=2)
