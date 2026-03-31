################################################################################
################             AAindex2 Module Tests             #################
################################################################################

import os
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
    def setUp(self):
        """ Initialise test data directory. """
        self.test_dir = 'test_data'

        #make test data directory, remove old one if exists
        if os.path.isdir(self.test_dir):
            os.rmdir(self.test_dir)
        os.mkdir(self.test_dir)

    def test_num_records(self):
        """ Test Case to check the correct number of records are present in the AAi2 database.
        To date, 94 records are present in the database. """
#1.)
        self.assertEqual(aaindex2.num_records(), 94,
            'Expected 94 records in AAi2, got {}.'.format(aaindex2.num_records()))

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

        self.assertEqual(record['description'], 'The PAM-120 matrix (Altschul, 1991)',
            'Unexpected description for ALTS910101, got {}.'.format(record['description']))
        self.assertEqual(record['pmid'], '1713145 2051488',
            'Unexpected pmid for ALTS910101, got {}.'.format(record['pmid']))
        self.assertEqual(record['notes'], '',
            'Unexpected notes for ALTS910101, got {}.'.format(record['notes']))

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
        record2 = aaindex2[index_code2]  # AZAE970101

        self.assertTrue(record2['description'].startswith('The single residue substitution matrix'),
            'Unexpected description for AZAE970101, got {}.'.format(record2['description']))
        self.assertEqual(record2['pmid'], '9488136',
            'Unexpected pmid for AZAE970101, got {}.'.format(record2['pmid']))
        self.assertIsInstance(record2['matrix'], dict,
            'matrix should be a dict, got {}.'.format(type(record2['matrix'])))
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
            'Unexpected description via dot notation, got {}.'.format(record.description))
        self.assertEqual(record.pmid, '1713145 2051488',
            'Unexpected pmid via dot notation, got {}.'.format(record.pmid))
        self.assertEqual(record.notes, '',
            'Unexpected notes via dot notation, got {}.'.format(record.notes))
        self.assertTrue(record.references.startswith('Altschul, S.F.'),
            'Unexpected references via dot notation, got {}.'.format(record.references))
        self.assertIsInstance(record.matrix, dict,
            'matrix via dot notation should be a dict.')
        self.assertIsInstance(record.row_order, list,
            'row_order via dot notation should be a list.')
#2.)
        record2 = aaindex2[index_code2]

        self.assertTrue(record2.description.startswith('The single residue substitution matrix'),
            'Unexpected description for AZAE970101, got {}.'.format(record2.description))
        self.assertEqual(record2.pmid, '9488136',
            'Unexpected pmid via dot notation, got {}.'.format(record2.pmid))

    def test_get(self):
        """ Test Case for the get() method: pairwise lookup, symmetry,
        and invalid input handling. """
        index_code1 = 'ALTS910101'
        index_code2 = 'AZAE970101'
#1.)
        #test known scalar values from the substitution matrix
        self.assertEqual(aaindex2.get(index_code1, 'A', 'A'), 3.0,
            'Expected A,A = 3.0 for ALTS910101, got {}.'.format(aaindex2.get(index_code1, 'A', 'A')))
        self.assertEqual(aaindex2.get(index_code1, 'A', 'R'), -3.0,
            'Expected A,R = -3.0 for ALTS910101, got {}.'.format(aaindex2.get(index_code1, 'A', 'R')))
        self.assertEqual(aaindex2.get(index_code2, 'L', 'L'), 9.0,
            'Expected L,L = 9.0 for AZAE970101, got {}.'.format(aaindex2.get(index_code2, 'L', 'L')))
        self.assertEqual(aaindex2.get(index_code2, 'V', 'I'), 9.0,
            'Expected V,I = 9.0 for AZAE970101, got {}.'.format(aaindex2.get(index_code2, 'V', 'I')))
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
        """ Test Case for search(), returning records whose description
        contains the given keyword(s). """
        description1 = 'substitution'
        description2 = 'blahblahblah'
        description3 = 'not a description'
        description4 = 1234
        description5 = True
#1.)
        #search for 'substitution' - should return 39 matching records
        search1 = aaindex2.search(description1)
        self.assertEqual(len(search1), 39,
            'Expected 39 records for substitution search, got {}.'.format(len(search1)))
        for index, val in search1.items():
            self.assertIn(description1.lower(), val['description'].lower(),
                'Search keyword not found in returned record description for {}.'.format(index))
#2.)
        #nonsense search terms should return empty dict
        search2 = aaindex2.search(description2)
        self.assertEqual(len(search2), 0,
            'Expected 0 records for nonsense search, got {}.'.format(len(search2)))
#3.)
        search3 = aaindex2.search(description3)
        self.assertEqual(len(search3), 0,
            'Expected 0 records for nonsense search, got {}.'.format(len(search3)))
#4.)
        #non-string/list input should raise TypeError
        with self.assertRaises(TypeError):
            aaindex2.search(description4)
        with self.assertRaises(TypeError):
            aaindex2.search(description5)

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
            'Expected {} in AAi2 record codes.'.format(index_code1))
        self.assertIn(index_code2, aaindex2.record_codes(),
            'Expected {} in AAi2 record codes.'.format(index_code2))
        self.assertIn(index_code3, aaindex2.record_codes(),
            'Expected {} in AAi2 record codes.'.format(index_code3))
        self.assertIn(index_code4, aaindex2.record_codes(),
            'Expected {} in AAi2 record codes.'.format(index_code4))
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
            'Expected 20 amino acids, got {}.'.format(len(amino_acids)))
        for aa in amino_acids:
            self.assertIn(aa, valid_amino_acids,
                'Amino acid {} not in list of valid amino acids.'.format(aa))

    def test_record_names(self):
        """ Test Case to check that all record descriptions are returned correctly. """
        record_names = aaindex2.record_names()
#1.)
        #number of names must equal number of records
        self.assertEqual(len(record_names), aaindex2.num_records(),
            'Number of record names should equal num_records(), got {}.'.format(len(record_names)))
#2.)
        #known descriptions should be present
        self.assertIn('The PAM-120 matrix (Altschul, 1991)', record_names,
            'Expected ALTS910101 description in record names.')
        self.assertTrue(any('single residue substitution matrix' in n for n in record_names),
            'Expected AZAE970101 description in record names.')

    def test_last_updated(self):
        """ Testing the last updated class attribute matches the known database version. """
        self.assertEqual(aaindex2.last_updated, "February 13, 2017",
            'Last updated value does not match expected, got {}.'.format(aaindex2.last_updated))

    def test_dunder_methods(self):
        """ Test Case for __len__, __contains__, __iter__, and __repr__ dunder methods. """
#1.) __len__
        self.assertEqual(len(aaindex2), 94,
            'Expected __len__ to return 94, got {}.'.format(len(aaindex2)))
#2.) __contains__
        self.assertIn('ALTS910101', aaindex2,
            'Expected ALTS910101 in aaindex2 via __contains__.')
        self.assertNotIn('BLAH999999', aaindex2,
            'Expected BLAH999999 to not be in aaindex2 via __contains__.')
#3.) __iter__
        codes_from_iter = list(aaindex2)
        self.assertEqual(len(codes_from_iter), 94,
            'Expected 94 codes from __iter__, got {}.'.format(len(codes_from_iter)))
        self.assertIn('ALTS910101', codes_from_iter,
            'Expected ALTS910101 in __iter__ output.')
#4.) __repr__
        self.assertEqual(repr(aaindex2), "AAIndex2(records=94, last_updated='February 13, 2017')",
            'Unexpected __repr__ output: {}.'.format(repr(aaindex2)))

    def test_values(self):
        """ Test Case for values() which returns the full matrix dict for a record. """
#1.)
        matrix = aaindex2.values('ALTS910101')
        self.assertIsInstance(matrix, dict,
            'values() should return a dict, got {}.'.format(type(matrix)))
        self.assertIn('A', matrix,
            'Expected amino acid A as key in values() result.')
        self.assertEqual(matrix['A']['A'], 3.0,
            'Expected A,A = 3.0 in values() result.')
#2.)
        #invalid record code should raise ValueError
        with self.assertRaises(ValueError):
            aaindex2.values('BLAH999999')

    def test_symmetry(self):
        """ Explicit symmetry test: for every record, (aa1, aa2) == (aa2, aa1). """
        amino_acids = aaindex2.amino_acids()
        #test symmetry on two representative records
        for code in ['ALTS910101', 'AZAE970101']:
            for aa1 in amino_acids:
                for aa2 in amino_acids:
                    val_forward = aaindex2.get(code, aa1, aa2)
                    val_reverse = aaindex2.get(code, aa2, aa1)
                    self.assertEqual(val_forward, val_reverse,
                        f'Symmetry violation for {code}: ({aa1},{aa2})={val_forward} != ({aa2},{aa1})={val_reverse}')

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

    def tearDown(self):
        """ Remove any test data or directories. """
        os.rmdir('test_data')

if __name__ == '__main__':
    #run all unit tests
    unittest.main(verbosity=2)
