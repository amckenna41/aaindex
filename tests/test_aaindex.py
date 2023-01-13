################################################################################
#################             AAIndex Module Tests             #################
################################################################################

import os
import numpy as np
import unittest
import requests
from importlib.metadata import metadata
from aaindex import aaindex1, __version__

class AAIndexTests(unittest.TestCase):

    def setUp(self):
        """ Inititalise AAIndex module variables and test data directory. """
        #make test data directory
        self.test_dir = 'test_data'

        if (os.path.isfile(self.test_dir)):
            os.rmdir(self.test_dir)
        os.mkdir(self.test_dir)

    def test_aaindex_metadata(self):
        """ Testing correct aaindex version and metadata. """
        self.assertEqual(__version__, "0.1.0",
            "aaindex version is not correct, got: {}".format(metadata("aaindex")['version']))
        self.assertEqual(metadata("aaindex")['name'], "aaindex", 
            "aaindex software name is not correct, got: {}".format(metadata("aaindex")['name']))
        # self.assertEqual(metadata("aaindex")['author'], "AJ McKenna, https://github.com/amckenna41", 
        #     "aaindex author is not correct, got: {}".format(metadata("aaindex")['author']))
        self.assertEqual(metadata("aaindex")['author-email'], "amckenna41@qub.ac.uk", 
            "aaindex author email is not correct, got: {}".format(metadata("aaindex")['author-email']))
        self.assertEqual(metadata('aaindex')['home-page'], "https://github.com/amckenna41/aaindex", 
            "aaindex home page url is not correct, got: {}".format(metadata('aaindex')['home-page']))
        # self.assertEqual(metadata('aaindex')['maintainer'], "AJ McKenna", 
        #     "aaindex maintainer is not correct, got: {}".format(metadata('aaindex')['maintainer']))
        self.assertEqual(metadata('aaindex')['license'], "MIT", 
            "aaindex license type is not correct, got: {}".format(metadata('aaindex')['license']))
        # self.assertEqual(metadata('aaindex')['summary'], 
        #     "aaindex is a lightweight Python software package for accessing the data in the various AAIndex databases, \
        #        which represent the physiochemical and biochemical properties of amino acids as numerical indices.", 
        #             "aaindex package summary is not correct, got: {}".format(metadata('aaindex')['summary']))
        # self.assertEqual(metadata('aaindex')['keywords'], 
        #     ["amino acid index", "bioinformatics", "protein engineering", "python", "pypi",
        #         "physiochemical properties", "biochemical properties", "proteins", "psp", "pysar"], 
        #             "aaindex keywords are not correct, got: {}".format(metadata('aaindex')['keywords']))

    @unittest.skip("Don't want to overload the FTP server each time tests are run.")
    def test_download(self):
        """ Test Case to check that the download functionality works for the required AAI1 file.
        The file will firstly be removed from the test directory and then redownloaded, and
        its presence in the directory will pass the test. """
#1.)
        #if AAI1 present in test dir then remove
        if (os.path.isfile(os.path.join(self.TEST_DIR, 'test_aaindex1'))):
            os.remove(os.path.join(self.TEST_DIR, 'test_aaindex1'))

        #assert that OSError exception raised when erroneous directory input to download func
        with self.assertRaises(OSError):
            aaindex1.download_aaindex()
#2.)
        #download AAI1 to local test directory
        aaindex1.download_aaindex()

        #verify download functionality has worked as desired
        self.assertTrue(os.path.isfile(os.path.join(self.TEST_DIR, 'aaindex1')),
            'AAI did not download correctly to the {} directory.'.format(self.TEST_DIR))

    @unittest.skip("Don't want to overload the FTP server each time tests are run.")
    def test_url(self):
        """ Test Case to check that the URL endpoints used for downloading the AAI databases
        return a 200 status code. """

        AA_INDEX1_URL_FTP = "https://www.genome.jp/ftp/db/community/aaindex/aaindex1"
        AA_INDEX1_URL_HTTPS = "https://www.genome.jp/ftp/db/community/aaindex/aaindex1"
        AA_INDEX2_URL = "https://www.genome.jp/ftp/db/community/aaindex/aaindex2"
        AA_INDEX3_URL = "https://www.genome.jp/ftp/db/community/aaindex/aaindex3"
        wrong_AA_INDEX_URL = "https://www.genome.jp/ftp/BLAH/BLAH/BLAH/BLAH"
#1.)
        #test URL endpoints for AAINDEX are active and give a 200 status code
        r = requests.get(AA_INDEX1_URL_FTP, allow_redirects = True)
        self.assertEqual(r.status_code, 200, 'URL not returning Status Code 200.')

        r = requests.get(AA_INDEX1_URL_HTTPS, allow_redirects = True)
        self.assertEqual(r.status_code, 200, 'URL not returning Status Code 200.')

        r = requests.get(AA_INDEX2_URL, allow_redirects = True)
        self.assertEqual(r.status_code, 200, 'URL not returning Status Code 200.')

        r = requests.get(AA_INDEX3_URL, allow_redirects = True)
        self.assertEqual(r.status_code, 200, 'URL not returning Status Code 200.')
#2.)
        r = requests.get(wrong_AA_INDEX_URL, allow_redirects = True)
        self.assertEqual(r.status_code, 404, 'Errenous URL, got 404 page not found error.')

        #maybe try requests.mock

    def test_get_amino_acids(self):
        """ Test Case to check that only valid Amino Acids are used within the AAIndex class. """

        valid_amino_acids = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
            'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        amino_acids = aaindex1.get_amino_acids()
#1.)
        for aa in amino_acids:
            self.assertIn(aa, valid_amino_acids, 'Amino Acid {} not in valid \
                amino acids: {}.'.format(aa, valid_amino_acids))

    def test_num_features(self):
        """ Test Case to check that the correct number of indices/records are present
        in the AAIndex object. To date, 566 indices are present in the database,
        features may be added to the AAI in time so the test takes this into account. """
#1.)
        self.assertTrue(aaindex1.get_num_records() >= 566, 'Incorrect number \
            of features found in AAI, wanted 566, got {}.'.format(aaindex1.get_num_records()))

    def test_feature_size(self):
        """ Test Case to check that the AAI has the correct dimensionality. Check that
        the number of keys in parsed JSON from the AAI is correct. The test takes
        into account if more features are added to the database in the future. """
#1.)
        self.assertTrue(len(aaindex1.record_codes())>=566, 'Incorrect number \
            of features found in AAI, wanted 566, got {}.'.format(len(aaindex1.record_codes())))

    def test_get_record_from_code(self):
        """ Test Case to check that the correct amino acid values are returned for
        the feature/index codes using the __getitem__ function. """
#1.)
        #initialise test feature codes and their correct amino acid values
        feature1 = 'AURR980103'
        feature1_vals = {'A': 1.05, 'L': 0.96, 'R': 0.81, 'K': 0.97, 'N': 0.91,
            'M': 0.99, 'D': 1.39, 'F': 0.95, 'C': 0.6, 'P': 1.05, 'Q': 0.87,
            'S': 0.96, 'E': 1.11, 'T': 1.03, 'G': 1.26, 'W': 1.06, 'H': 1.43,
            'Y': 0.94, 'I': 0.95, 'V': 0.62, '-': 0}
        feature2 = 'FAUJ880113'
        feature2_vals = {'A': 4.76, 'L': 4.79, 'R': 4.3, 'K': 4.27, 'N': 3.64,
            'M': 4.25, 'D': 5.69, 'F': 4.31, 'C': 3.67, 'P': 0.0, 'Q': 4.54,
            'S': 3.83, 'E': 5.48, 'T': 3.87, 'G': 3.77, 'W': 4.75, 'H': 2.84,
            'Y': 4.3, 'I': 4.81, 'V': 4.86, '-': 0}
        feature3 = 'NAGK730101'
        feature3_vals = {'A': 1.29, 'L': 1.23, 'R': 0.83, 'K': 1.23, 'N': 0.77,
            'M': 1.23, 'D': 1.0, 'F': 1.23, 'C': 0.94, 'P': 0.7, 'Q': 1.1,
            'S': 0.78, 'E': 1.54, 'T': 0.87, 'G': 0.72, 'W': 1.06, 'H': 1.29,
            'Y': 0.63, 'I': 0.94, 'V': 0.97, '-': 0}
        feature4 = 'ROBB760105'
        feature4_vals = {'A': -2.3, 'L': 2.3, 'R': 0.4, 'K': -3.3, 'N': -4.1,
            'M': 2.3, 'D': -4.4, 'F': 2.6, 'C': 4.4, 'P': -1.8, 'Q': 1.2,
            'S': -1.7, 'E': -5.0, 'T': 1.3, 'G': -4.2, 'W': -1.0, 'H': -2.5,
            'Y': 4.0, 'I': 6.7, 'V': 6.8, '-': 0}
#1.)
        #get amino acid values for inputted feature/index codes
        record = aaindex1[feature1]

        self.assertIn('description', list(record.keys()))
        self.assertIn('refs', list(record.keys()))
        self.assertIn('notes', list(record.keys()))
        self.assertIn('values', list(record.keys()))

        self.assertEqual(feature1_vals, record['values'], 'Amino acid values gotten \
            from AAI do not match the desired values.')
        self.assertEqual(record['description'], "Normalized positional residue frequency at helix termini N' (Aurora-Rose, 1998)")
        self.assertEqual(record['notes'],'')
        self.assertTrue(record['refs'].startswith("Aurora, R. and Rose, G."))
#2.)
        record = aaindex1[feature2]

        self.assertIn('description', list(record.keys()))
        self.assertIn('refs', list(record.keys()))
        self.assertIn('notes', list(record.keys()))
        self.assertIn('values', list(record.keys()))

        self.assertEqual(feature2_vals, record['values'], 'Amino acid values gotten \
            from AAI do not match the desired values.')
        self.assertEqual(record['description'], 'pK-a(RCOOH) (Fauchere et al., 1988)')
        self.assertEqual(record['notes'],'')
        self.assertTrue(record['refs'].startswith("Fauchere, J.L., Charton, M"))
#3.)
        record = aaindex1[feature3]

        self.assertIn('description', list(record.keys()))
        self.assertIn('refs', list(record.keys()))
        self.assertIn('notes', list(record.keys()))
        self.assertIn('values', list(record.keys()))

        self.assertEqual(feature3_vals, record['values'], 'Amino acid values gotten \
            from AAI do not match the desired values.')
        self.assertEqual(record['description'], 'Normalized frequency of alpha-helix (Nagano, 1973)')
        self.assertEqual(record['notes'],'')
        self.assertTrue(record['refs'].startswith("Nagano, K."))
#4.)
        record = aaindex1[feature4]
        self.assertIn('description', list(record.keys()))
        self.assertIn('refs', list(record.keys()))
        self.assertIn('notes', list(record.keys()))
        self.assertIn('values', list(record.keys()))

        self.assertEqual(feature4_vals, record['values'], 'Amino acid values gotten \
            from AAI do not match the desired values.')
        self.assertEqual(record['description'], 'Information measure for extended (Robson-Suzuki, 1976)')
        self.assertEqual(record['notes'],'')
        self.assertTrue(record['refs'].startswith("Robson, B. and Suzuki, E. 'Conformational properties"))
#5.)
        #testing value error raised when errenous indices put into function
        feature5 = "ABCDEFGH"
        feature6 = "123456"
        feature7 = "blahblahblah"

        with self.assertRaises(ValueError):
            feature_vals = aaindex1[feature5]

        with self.assertRaises(ValueError):
            feature_vals = aaindex1[feature6]

        with self.assertRaises(ValueError):
            feature_vals = aaindex1[feature7]

    def test_get_record_from_name(self):
        """ 

        """
        feature1 = 'YUTK870101'
        # Unfolding Gibbs energy in water, pH7.0 (Yutani et al., 1987)
        feature2 = 'ZIMJ680103'
        # Polarity (Zimmerman et al., 1968)
        feature3 = 'ROSM880103'
        # Loss of Side chain hydropathy by helix formation (Roseman, 1988)

        error_feature1 = 'blahblahblah'
        error_feature2 = 'not a record name'
        error_feature3 = 'also not a record name'

    def test_record_codes(self):
        """ Test Case to check that correct feature/record codes are in the parsed JSON
        of all AAI records. Also testing that random erroneous feature codes are
        not present. """

        #testing actual index codes
        feature1 = 'VHEG790101'
        feature2 = 'PONP800107'
        feature3 = 'OOBM770102'
        feature4 = 'NADH010101'
#1.)
        #testing index codes are in the AAI1
        self.assertIn(feature1, aaindex1.record_codes(), 'Index {} \
            not found in list of available indices.'.format(feature1))
        self.assertIn(feature2, aaindex1.record_codes(), 'Index {} \
            not found in list of available indices'.format(feature2))
        self.assertIn(feature3, aaindex1.record_codes(), 'Index {} \
            not found in list of available indices'.format(feature3))
        self.assertIn(feature4, aaindex1.record_codes(), 'Index {} \
            not found in list of available indices'.format(feature4))
#2.)
        #testing bogus index codes
        feature5 = 'ABC1234'
        feature6 = 'ABC12345'
        feature7 = 'ABC123456'
        feature8 = 'ABC1234567'

        #testing errenous index codes are not in the AAI1
        self.assertNotIn(feature5, aaindex1.record_codes(), 'Index {} \
            erroneously found in list of available indices.'.format(feature5))
        self.assertNotIn(feature6, aaindex1.record_codes(), 'Index {} \
            erroneously found in list of available indices.'.format(feature6))
        self.assertNotIn(feature7, aaindex1.record_codes(), 'Index {} \
            erroneously found in list of available indices.'.format(feature7))
        self.assertNotIn(feature8, aaindex1.record_codes(), 'Index {} \
            erroneously found in list of available indices.'.format(feature8))

    def test_category(self):
        """ Test Case to check that the correct category value is returned for each index. """

        feature1 = 'TANS770109'
        feature2 = 'RADA880107'
        feature3 = 'ZIMJ680103'
        feature4 = 'WOLS870102'
        err_feature1 = "blahblahblah"
        err_feature2 = "12345"
        err_feature3 = 100.902
#1.)
        cat1 = aaindex1.get_category_from_record(feature1)
        cat2 = aaindex1.get_category_from_record(feature2)
        cat3 = aaindex1.get_category_from_record(feature3)
        cat4 = aaindex1.get_category_from_record(feature4)

        self.assertEqual(cat1, 'sec_struct', 'Incorrect category returned, got {} \
            instead of {} '.format(cat1, 'sec_struct'))
        self.assertEqual(cat2, 'hydrophobic','Incorrect category returned, got {} \
            instead of {} '.format(cat2, 'hydrophobic'))
        self.assertEqual(cat3, 'polar','Incorrect category returned, got {} \
            instead of {} '.format(cat3, 'polar'))
        self.assertEqual(cat4, 'meta','Incorrect category returned, got {} \
            instead of {} '.format(cat4, 'meta'))
#2.)
        with (self.assertRaises(ValueError)):
            cat5 = aaindex1.get_category_from_record(err_feature1)

        with (self.assertRaises(ValueError)):
            cat6 = aaindex1.get_category_from_record(err_feature2)

        with (self.assertRaises(TypeError)):
            cat7 = aaindex1.get_category_from_record(err_feature3)

    def test_aaindex_names(self):
        """ Test Case to check that the correct names/description is returned for each index. """

        names = aaindex1.get_record_names()

        name1 = names[50]
        name2 = names[140]
        name3 = names[368]
        name4 = names[560]

        self.assertEqual(len(names), aaindex1.get_num_records(), 'Total number \
            of names in the AAI database should equal the number of records')
        self.assertEqual(name1, 'alpha-NH chemical shifts (Bundi-Wuthrich, 1979)',
            'Incorrect index name returned, got {} instead of {}'.format(name1,
                'alpha-NH chemical shifts (Bundi-Wuthrich, 1979)'))
        self.assertTrue(name2.startswith('Number of full nonbonding'),
            'Incorrect index name returned, got {} instead of {}'.format(name2,
                'Number of full nonbonding'))
        self.assertIn('alpha-helix in alpha/beta class', name3, 'String {} not found in index name: \
            {}'.format('alpha-helix in alpha/beta class',name3))
        self.assertTrue(name4.endswith('Buriability (Zhou-Zhou, 2004)'), 'Index does \
            not end with {}, got {}'.format('Buriability (Zhou-Zhou, 2004)', name4))

    def test_aaindex_refs(self):
        """ Test Case to check that the correct references are returned for each index. """

        feature1 = 'VELV850101'
        feature2 = 'QIAN880139'
        feature3 = 'NAKH900112'
        feature4 = 'LEVM760103'
        feature5 = "blah_blah"
        feature6 = 123
#1.)
        ref1 = aaindex1[feature1]['refs']
        ref2 = aaindex1[feature2]['refs']
        ref3 = aaindex1[feature3]['refs']
        ref4 = aaindex1[feature4]['refs']

        self.assertTrue(ref1.startswith('Veljkovic'), 'Incorrect index reference, \
            got {} instead of {} '.format(ref1,'' ))
        self.assertIn("globular proteins using neural network models", ref2)
        self.assertIn("hydrophobicity of amino acid composition of mitochondrial proteins", ref3)
        self.assertTrue(ref4.endswith('Nucleic Acids Res. 28, 374 (2000).'))
#2.)
        with (self.assertRaises(ValueError)):
            ref5 = aaindex1[feature5]['refs']

        with (self.assertRaises(TypeError)):
            ref6 = aaindex1[feature6]['refs']

    def test_aaindex_encoding(self):
        """ Testing amino acids one-hot encoding. """
#1.)
        encoding = aaindex1.get_amino_acids_encoding()
        self.assertEqual(encoding.shape, (20,20))
        self.assertIsInstance(encoding, np.ndarray)

    def tearDown(self):
        """ Remove any test data or directories. """
        #remove test data directory
        os.rmdir('test_data')

if __name__ == '__main__':
    #run all unit tests
    unittest.main(verbosity=2)