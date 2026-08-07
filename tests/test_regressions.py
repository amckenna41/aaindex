################################################################################
################              Bug Regression Tests              ################
################################################################################

import os
import shutil
import tempfile
import unittest
import warnings

from aaindex import aaindex1, aaindex2, aaindex3
from aaindex.aaindex1 import AAIndex1

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'aaindex', 'data')


class MalformedMatrixHeaderTests(unittest.TestCase):
    """B1: record DOSZ010101 has a 'cols rows =' typo in the upstream AAindex2 source.
    Positional parsing read '=' as the entire column order and silently discarded the
    whole 22x20 matrix. """

    def test_dosz010101_matrix_is_intact(self):
        record = aaindex2['DOSZ010101']
        self.assertEqual(len(record['col_order']), 20,
            f"DOSZ010101 should have 20 columns, got {len(record['col_order'])} — "
            f"the malformed 'cols rows =' header was not handled.")
        self.assertNotIn('=', record['col_order'],
            "'=' parsed as an amino acid column; the rows/cols header was read positionally.")
        self.assertIsNotNone(aaindex2.get('DOSZ010101', 'A', 'R'),
            'DOSZ010101 pairwise lookup returned None; matrix was discarded during parsing.')

    def test_dosz010101_is_not_flagged_symmetric(self):
        #the truncated 22x1 matrix made the is_full check fail, wrongly marking it symmetric
        self.assertFalse(aaindex2['DOSZ010101']['is_symmetric'],
            'DOSZ010101 is a full asymmetric matrix and must not be flagged symmetric.')
        self.assertNotEqual(
            aaindex2.get('DOSZ010101', 'A', 'R'),
            aaindex2.get('DOSZ010101', 'R', 'A'),
            'DOSZ010101 transposed lookups should differ for an asymmetric matrix.')

    def test_symmetry_counts_match_documented_figures(self):
        #these exact figures appear in README.md; DOSZ010101 being misparsed shifted them,
        #so pin them here rather than letting the docs drift away from the parser again
        for db, name, total, symmetric in ((aaindex2, 'aaindex2', 94, 67), (aaindex3, 'aaindex3', 47, 44)):
            actual_symmetric = sum(1 for code in db.record_codes() if db[code]['is_symmetric'])
            self.assertEqual(len(db), total, f'{name} should hold {total} records.')
            self.assertEqual(actual_symmetric, symmetric,
                f'{name} should have {symmetric} symmetric / {total - symmetric} non-symmetric '
                f'matrices (README.md quotes these figures), got {actual_symmetric} symmetric.')

    def test_every_record_has_amino_acid_columns(self):
        for db, name in ((aaindex2, 'aaindex2'), (aaindex3, 'aaindex3')):
            for code in db.record_codes():
                col_order = db[code]['col_order']
                self.assertGreaterEqual(len(col_order), 20,
                    f'{name} record {code} has only {len(col_order)} columns.')


class MutationIsolationTests(unittest.TestCase):
    """B2: accessors returned live references into the singleton, so mutating a
    returned record silently corrupted the database for the whole process. """

    def test_getitem_is_isolated(self):
        original = aaindex1['KYTJ820101']['values']['A']
        aaindex1['KYTJ820101']['values']['A'] = 999.0
        self.assertEqual(aaindex1['KYTJ820101']['values']['A'], original,
            '__getitem__ leaked a reference to the shared values dict.')

    def test_matrix_getitem_is_isolated(self):
        original = aaindex2.get('ALTS910101', 'A', 'A')
        aaindex2['ALTS910101']['matrix']['A']['A'] = 42.0
        self.assertEqual(aaindex2.get('ALTS910101', 'A', 'A'), original,
            '__getitem__ leaked a reference to the shared matrix.')

    def test_to_dict_is_isolated(self):
        original = aaindex1['KYTJ820101']['values']['C']
        aaindex1.to_dict('KYTJ820101')['KYTJ820101']['values']['C'] = -1.0
        aaindex1.to_dict()['KYTJ820101']['values']['C'] = -2.0
        self.assertEqual(aaindex1['KYTJ820101']['values']['C'], original,
            'to_dict() returned a shallow copy sharing the nested values dict.')

    def test_search_results_are_isolated(self):
        original = aaindex1['ARGP820101']['values']['A']
        aaindex1.search('hydrophobicity')['ARGP820101']['values']['A'] = -5.0
        self.assertEqual(aaindex1['ARGP820101']['values']['A'], original,
            'search() leaked references to the shared records.')

    def test_category_results_are_isolated(self):
        results = aaindex1.get_record_by_category('charge')
        code = next(iter(results))
        original = aaindex1[code]['values']['A']
        results[code]['values']['A'] = -5.0
        self.assertEqual(aaindex1[code]['values']['A'], original,
            'get_record_by_category() leaked references to the shared records.')

    def test_amino_acids_cache_is_isolated(self):
        for db, name in ((aaindex1, 'aaindex1'), (aaindex2, 'aaindex2'), (aaindex3, 'aaindex3')):
            expected = len(db.amino_acids())
            db.amino_acids().append('ZZ')
            self.assertEqual(len(db.amino_acids()), expected,
                f'{name}.amino_acids() returned the cache list itself.')


class ParserRobustnessTests(unittest.TestCase):
    """B3/B4/B5: malformed input and a damaged cache must degrade with a warning
    rather than raising IndexError / KeyError / JSONDecodeError. """

    def _aaindex1_from_raw(self, raw, categories):
        tmp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp_dir, True)
        os.makedirs(os.path.join(tmp_dir, 'data'))
        with open(os.path.join(tmp_dir, 'data', 'aaindex1'), 'w') as f:
            f.write(raw)
        db = AAIndex1()
        db.aaindex_module_path = tmp_dir
        db.categories = categories
        return db

    def test_short_value_row_warns_instead_of_indexerror(self):
        raw = ('H TEST01\nD short second value row\nI    A/L     R/K     N/M\n'
               '    1.0     2.0     3.0\n   11.0\n//\n')
        db = self._aaindex1_from_raw(raw, {'TEST01': 'charge'})
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            parsed = db._parse_aaindex()
        self.assertIn('TEST01', parsed)
        self.assertEqual(parsed['TEST01']['values']['A'], 1.0)
        self.assertTrue(any('mismatched I-block' in str(w.message) for w in caught),
            'A mismatched I-block should emit a warning.')

    def test_blank_lines_are_skipped_silently(self):
        raw = ('H TEST01\nD ok\nI    A/L     R/K\n    1.0     2.0\n    3.0     4.0\n//\n'
               '\n\n'
               'H TEST02\nD ok too\nI    A/L     R/K\n    5.0     6.0\n    7.0     8.0\n//')
        db = self._aaindex1_from_raw(raw, {'TEST01': 'charge', 'TEST02': 'charge'})
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            parsed = db._parse_aaindex()
        self.assertEqual(sorted(parsed), ['TEST01', 'TEST02'])
        self.assertFalse(any('Unknown field code' in str(w.message) for w in caught),
            'Blank lines should not be reported as unknown field codes.')

    def test_record_missing_from_category_file_warns(self):
        raw = ('H TEST01\nD ok\nI    A/L     R/K\n    1.0     2.0\n    3.0     4.0\n//\n'
               'H TEST02\nD uncategorised\nI    A/L     R/K\n    5.0     6.0\n    7.0     8.0\n//\n')
        db = self._aaindex1_from_raw(raw, {'TEST01': 'charge'})
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            parsed = db._parse_aaindex()
        self.assertEqual(parsed['TEST02']['category'], 'unknown',
            'A record absent from the category file should fall back to "unknown".')
        self.assertTrue(any('no entry in the category file' in str(w.message) for w in caught),
            'An uncategorised record should emit a warning.')

    def test_truncated_cache_falls_back_to_raw_file(self):
        tmp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp_dir, True)
        shutil.copytree(DATA_DIR, os.path.join(tmp_dir, 'data'),
            ignore=shutil.ignore_patterns('*.json'))
        cache_path = os.path.join(tmp_dir, 'data', 'aaindex1.json')
        with open(cache_path, 'w') as f:
            f.write('{"A": {"description": tru')   # simulates an interrupted write
        db = AAIndex1()
        db.aaindex_module_path = tmp_dir
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            self.assertEqual(len(db), 566,
                'A truncated cache should fall back to parsing the raw file.')
        self.assertTrue(any('unreadable' in str(w.message) for w in caught),
            'An unreadable cache should emit a warning.')

    def test_cache_write_is_atomic_and_readable(self):
        tmp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp_dir, True)
        shutil.copytree(DATA_DIR, os.path.join(tmp_dir, 'data'),
            ignore=shutil.ignore_patterns('*.json'))
        db = AAIndex1()
        db.aaindex_module_path = tmp_dir
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            db._parse_aaindex()
        cache_path = os.path.join(tmp_dir, 'data', 'aaindex1.json')
        self.assertTrue(os.path.isfile(cache_path), 'Cache file was not written.')
        self.assertEqual(oct(os.stat(cache_path).st_mode)[-3:], '644',
            'Cache must be world-readable, not the 0600 mkstemp default.')
        leftovers = [f for f in os.listdir(os.path.join(tmp_dir, 'data')) if f.endswith('.tmp')]
        self.assertEqual(leftovers, [], f'Temporary cache files left behind: {leftovers}.')


class ContainsConsistencyTests(unittest.TestCase):
    """B6: `code in db` did not normalise case/whitespace, so it disagreed with db[code]. """

    def test_contains_matches_getitem(self):
        for db, code in ((aaindex1, 'KYTJ820101'), (aaindex2, 'ALTS910101'), (aaindex3, 'TANS760101')):
            for variant in (code, code.lower(), f'  {code}  ', f'  {code.lower()}  '):
                self.assertIn(variant, db,
                    f'{variant!r} resolves via __getitem__ but not via __contains__.')
                db[variant]   # must not raise
            self.assertNotIn('BLAH999999', db)

    def test_contains_with_non_string_is_false(self):
        for db in (aaindex1, aaindex2, aaindex3):
            self.assertNotIn(12345, db, '__contains__ should return False for non-string input.')
            self.assertNotIn(None, db, '__contains__ should return False for non-string input.')


if __name__ == '__main__':
    unittest.main(verbosity=2)
