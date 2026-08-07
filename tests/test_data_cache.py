################################################################################
################          Data Cache Consistency Tests          ################
################################################################################

import json
import os
import shutil
import tempfile
import unittest
import warnings

from aaindex.aaindex1 import AAIndex1
from aaindex.aaindex2 import AAIndex2
from aaindex.aaindex3 import AAIndex3

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'aaindex', 'data')


class DataCacheTests(unittest.TestCase):
    """The committed *.json caches are shipped as package data. If they drift from the raw
    AAindex files a fresh install serves different data than a from-source parse, so this
    asserts they stay in sync. Regenerate with: rm aaindex/data/*.json; python -c "import aaindex"
    """

    def _assert_cache_matches_raw(self, cls, filename):
        cached_path = os.path.join(DATA_DIR, f'{filename}.json')
        self.assertTrue(os.path.isfile(cached_path),
            f'{filename}.json is missing — it must be committed so installs do not re-parse.')
        with open(cached_path) as f:
            cached = json.load(f)

        #parse the raw file in an isolated copy of the data dir so the real cache is untouched
        tmp_dir = tempfile.mkdtemp()
        try:
            shutil.copytree(DATA_DIR, os.path.join(tmp_dir, 'data'),
                ignore=shutil.ignore_patterns('*.json'))
            db = cls()
            db.aaindex_module_path = tmp_dir
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fresh = db._parse_aaindex()
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

        self.assertEqual(sorted(cached), sorted(fresh),
            f'{filename}.json record codes differ from a fresh parse of the raw file.')
        for code in fresh:
            self.assertEqual(
                json.loads(json.dumps(cached[code], sort_keys=True)),
                json.loads(json.dumps(fresh[code], sort_keys=True)),
                f'{filename}.json record {code} differs from a fresh parse of the raw file.')

    def test_aaindex1_cache_matches_raw(self):
        self._assert_cache_matches_raw(AAIndex1, 'aaindex1')

    def test_aaindex2_cache_matches_raw(self):
        self._assert_cache_matches_raw(AAIndex2, 'aaindex2')

    def test_aaindex3_cache_matches_raw(self):
        self._assert_cache_matches_raw(AAIndex3, 'aaindex3')


if __name__ == '__main__':
    unittest.main(verbosity=2)
