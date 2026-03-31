################################################################################
################             AAindex3 Module Tests             #################
################################################################################

import os
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
    def setUp(self):
        """ Initialise test data directory. """
        self.test_dir = 'test_data'

        #make test data directory, remove old one if exists
        if os.path.isdir(self.test_dir):
            os.rmdir(self.test_dir)
        os.mkdir(self.test_dir)

    def test_num_records(self):
        """ Test Case to check the correct number of records are present in the AAi3 database.
        To date, 47 records are present in the database. """
#1.)
        self.assertEqual(aaindex3.num_records(), 47,
            'Expected 47 records in AAi3, got {}.'.format(aaindex3.num_records()))

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
            'Unexpected description for TANS760101, got {}.'.format(record['description']))
        self.assertEqual(record['pmid'], '1004017',
            'Unexpected pmid for TANS760101, got {}.'.format(record['pmid']))

        #verify matrix is a nested dict keyed by amino acid
        self.assertIsInstance(record['matrix'], dict,
            'matrix should be a dict, got {}.'.format(type(record['matrix'])))
        self.assertIn('A', record['matrix'],
            'Expected amino acid A as a key in the matrix dict.')
        self.assertIsInstance(record['matrix']['A'], dict,
            'Inner matrix values should be dicts keyed by amino acid.')

        #verify row_order and col_order are lists of 20 amino acids
        self.assertIsInstance(record['row_order'], list,
            'row_order should be a list, got {}.'.format(type(record['row_order'])))
        self.assertEqual(len(record['row_order']), 20,
            'row_order should contain 20 amino acids, got {}.'.format(len(record['row_order'])))
        self.assertEqual(len(record['col_order']), 20,
            'col_order should contain 20 amino acids, got {}.'.format(len(record['col_order'])))
#2.)
        record2 = aaindex3[index_code2]  # BETM990101

        self.assertEqual(record2['description'], 'Modified version of the Miyazawa-Jernigan transfer energy',
            'Unexpected description for BETM990101, got {}.'.format(record2['description']))
        self.assertEqual(record2['pmid'], '10048329',
            'Unexpected pmid for BETM990101, got {}.'.format(record2['pmid']))
        self.assertIsInstance(record2['matrix'], dict,
            'matrix should be a dict for BETM990101, got {}.'.format(type(record2['matrix'])))
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
            'Unexpected description via dot notation, got {}.'.format(record.description))
        self.assertEqual(record.pmid, '1004017',
            'Unexpected pmid via dot notation, got {}.'.format(record.pmid))
        self.assertTrue(record.references.startswith('Tanaka, S. and Scheraga, H.A.'),
            'Unexpected references via dot notation, got {}.'.format(record.references))
        self.assertIsInstance(record.matrix, dict,
            'matrix via dot notation should be a dict.')
        self.assertIsInstance(record.row_order, list,
            'row_order via dot notation should be a list.')
#2.)
        record2 = aaindex3[index_code2]

        self.assertEqual(record2.description, 'Modified version of the Miyazawa-Jernigan transfer energy',
            'Unexpected description for BETM990101, got {}.'.format(record2.description))
        self.assertEqual(record2.pmid, '10048329',
            'Unexpected pmid via dot notation, got {}.'.format(record2.pmid))

    def test_get(self):
        """ Test Case for the get() method: pairwise lookup, symmetry, NA handling,
        and invalid input handling. """
        index_code1 = 'TANS760101'
        index_code2 = 'BETM990101'
        index_code3 = 'ROBB790102'
#1.)
        #test known scalar values from the contact potential matrix
        self.assertEqual(aaindex3.get(index_code1, 'A', 'A'), -2.6,
            'Expected A,A = -2.6 for TANS760101, got {}.'.format(aaindex3.get(index_code1, 'A', 'A')))
        self.assertEqual(aaindex3.get(index_code1, 'A', 'R'), -3.4,
            'Expected A,R = -3.4 for TANS760101, got {}.'.format(aaindex3.get(index_code1, 'A', 'R')))
        self.assertEqual(aaindex3.get(index_code2, 'A', 'A'), -0.2,
            'Expected A,A = -0.2 for BETM990101, got {}.'.format(aaindex3.get(index_code2, 'A', 'A')))
        self.assertEqual(aaindex3.get(index_code2, 'L', 'V'), -0.8,
            'Expected L,V = -0.8 for BETM990101, got {}.'.format(aaindex3.get(index_code2, 'L', 'V')))
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
        """ Test Case for search(), returning records whose description
        contains the given keyword(s). """
        description1 = 'contact potential'
        description2 = 'blahblahblah'
        description3 = 'not a description'
        description4 = 1234
        description5 = True
#1.)
        #search for 'contact potential' - should return 3 matching records
        search1 = aaindex3.search(description1)
        self.assertEqual(len(search1), 3,
            'Expected 3 records for contact potential search, got {}.'.format(len(search1)))
        for index, val in search1.items():
            self.assertIn(description1.lower(), val['description'].lower(),
                'Search keyword not found in returned record description for {}.'.format(index))
#2.)
        #nonsense search terms should return empty dict
        search2 = aaindex3.search(description2)
        self.assertEqual(len(search2), 0,
            'Expected 0 records for nonsense search, got {}.'.format(len(search2)))
#3.)
        search3 = aaindex3.search(description3)
        self.assertEqual(len(search3), 0,
            'Expected 0 records for nonsense search, got {}.'.format(len(search3)))
#4.)
        #non-string/list input should raise TypeError
        with self.assertRaises(TypeError):
            aaindex3.search(description4)
        with self.assertRaises(TypeError):
            aaindex3.search(description5)

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
            'Expected {} in AAi3 record codes.'.format(index_code1))
        self.assertIn(index_code2, aaindex3.record_codes(),
            'Expected {} in AAi3 record codes.'.format(index_code2))
        self.assertIn(index_code3, aaindex3.record_codes(),
            'Expected {} in AAi3 record codes.'.format(index_code3))
        self.assertIn(index_code4, aaindex3.record_codes(),
            'Expected {} in AAi3 record codes.'.format(index_code4))
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
            'Expected 20 amino acids, got {}.'.format(len(amino_acids)))
        for aa in amino_acids:
            self.assertIn(aa, valid_amino_acids,
                'Amino acid {} not in list of valid amino acids.'.format(aa))

    def test_record_names(self):
        """ Test Case to check that all record descriptions are returned correctly. """
        record_names = aaindex3.record_names()
#1.)
        #number of names must equal number of records
        self.assertEqual(len(record_names), aaindex3.num_records(),
            'Number of record names should equal num_records(), got {}.'.format(len(record_names)))
#2.)
        #known descriptions should be present in the list
        self.assertIn('Statistical contact potential derived from 25 x-ray protein structures', record_names,
            'Expected TANS760101 description in record names.')
        self.assertIn('Modified version of the Miyazawa-Jernigan transfer energy', record_names,
            'Expected BETM990101 description in record names.')

    def test_last_updated(self):
        """ Testing the last updated class attribute matches the known database version. """
        self.assertEqual(aaindex3.last_updated, "February 13, 2017",
            'Last updated value does not match expected, got {}.'.format(aaindex3.last_updated))

    def test_dunder_methods(self):
        """ Test Case for __len__, __contains__, __iter__, and __repr__ dunder methods. """
#1.) __len__
        self.assertEqual(len(aaindex3), 47,
            'Expected __len__ to return 47, got {}.'.format(len(aaindex3)))
#2.) __contains__
        self.assertIn('TANS760101', aaindex3,
            'Expected TANS760101 in aaindex3 via __contains__.')
        self.assertNotIn('BLAH999999', aaindex3,
            'Expected BLAH999999 to not be in aaindex3 via __contains__.')
#3.) __iter__
        codes_from_iter = list(aaindex3)
        self.assertEqual(len(codes_from_iter), 47,
            'Expected 47 codes from __iter__, got {}.'.format(len(codes_from_iter)))
        self.assertIn('TANS760101', codes_from_iter,
            'Expected TANS760101 in __iter__ output.')
#4.) __repr__
        self.assertEqual(repr(aaindex3), "AAIndex3(records=47, last_updated='February 13, 2017')",
            'Unexpected __repr__ output: {}.'.format(repr(aaindex3)))

    def test_values(self):
        """ Test Case for values() which returns the full matrix dict for a record. """
#1.)
        matrix = aaindex3.values('TANS760101')
        self.assertIsInstance(matrix, dict,
            'values() should return a dict, got {}.'.format(type(matrix)))
        self.assertIn('A', matrix,
            'Expected amino acid A as key in values() result.')
        self.assertEqual(matrix['A']['A'], -2.6,
            'Expected A,A = -2.6 in values() result.')
#2.)
        #invalid record code should raise ValueError
        with self.assertRaises(ValueError):
            aaindex3.values('BLAH999999')

    def test_symmetry(self):
        """ Explicit symmetry test: for every record, (aa1, aa2) == (aa2, aa1). """
        amino_acids = aaindex3.amino_acids()
        #test symmetry on two representative records
        for code in ['TANS760101', 'BETM990101']:
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

    def tearDown(self):
        """ Remove any test data or directories. """
        os.rmdir('test_data')

if __name__ == '__main__':
    #run all unit tests
    unittest.main(verbosity=2)
