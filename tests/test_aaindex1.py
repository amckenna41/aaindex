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
        self.test_dir = 'test_data'

        #make test data directory, remove old one if exists
        if (os.path.isfile(self.test_dir)):
            os.rmdir(self.test_dir)
        os.mkdir(self.test_dir)

    def test_aaindex_metadata(self):
        """ Testing correct aaindex version and metadata. """
        self.assertEqual(__version__, "1.0.4",
            "aaindex version is not correct, got: {}".format(metadata("aaindex")['version']))
        self.assertEqual(metadata("aaindex")['name'], "aaindex", 
            "aaindex software name is not correct, got: {}".format(metadata("aaindex")['name']))
        self.assertEqual(metadata("aaindex")['author'], "AJ McKenna, https://github.com/amckenna41", 
            "aaindex author is not correct, got: {}".format(metadata("aaindex")['author']))
        self.assertEqual(metadata("aaindex")['author-email'], "amckenna41@qub.ac.uk", 
            "aaindex author email is not correct, got: {}".format(metadata("aaindex")['author-email']))
        self.assertEqual(metadata('aaindex')['home-page'], "https://github.com/amckenna41/aaindex", 
            "aaindex home page url is not correct, got: {}".format(metadata('aaindex')['home-page']))
        self.assertEqual(metadata('aaindex')['maintainer'], "AJ McKenna", 
            "aaindex maintainer is not correct, got: {}".format(metadata('aaindex')['maintainer']))
        self.assertEqual(metadata('aaindex')['license'], "MIT", 
            "aaindex license type is not correct, got: {}".format(metadata('aaindex')['license']))
        self.assertEqual(metadata('aaindex')['summary'], 
            "aaindex is a lightweight Python software package for accessing the data in the various AAIndex databases,"
            " which represent the physiochemical and biochemical properties of amino acids as numerical indices.", 
                    "aaindex package summary is not correct, got: {}".format(metadata('aaindex')['summary']))
        self.assertEqual(metadata('aaindex')['keywords'], 
            "amino acid index,bioinformatics,protein engineering,python,pypi,physiochemical properties,"
            "biochemical properties,proteins,psp,pysar", 
                "aaindex keywords are not correct, got: {}".format(metadata('aaindex')['keywords']))

    @unittest.skip("Don't want to overload the FTP server each time tests are run.")
    def test_download(self):
        """ Test Case to check that the download functionality works for the required AAIndex files.
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
            'AAI1 did not download correctly to the {} directory.'.format(self.TEST_DIR))

    @unittest.skip("Don't want to overload the FTP server each time tests are run.")
    def test_url(self):
        """ Test Case to check that the URL endpoints used for downloading the AAI1 databases
        return a 200 status code. """

        AA_INDEX1_URL_FTP = "https://www.genome.jp/ftp/db/community/aaindex/aaindex1"
        AA_INDEX1_URL_HTTPS = "https://www.genome.jp/ftp/db/community/aaindex/aaindex1"
        AA_INDEX2_URL = "https://www.genome.jp/ftp/db/community/aaindex/aaindex2"
        AA_INDEX3_URL = "https://www.genome.jp/ftp/db/community/aaindex/aaindex3"
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

    def test_amino_acids(self):
        """ Test Case to check that only valid Amino Acids are used within the AAIndex class. """
        valid_amino_acids = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
            'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        amino_acids = aaindex1.amino_acids()
#1.)
        for aa in amino_acids:
            self.assertIn(aa, valid_amino_acids, 'Amino Acid {} not in valid \
                amino acids: {}.'.format(aa, valid_amino_acids))

    def test_num_records(self):
        """ Test Case to check that the correct number of indices/records are present
        in the AAIndex object. To date, 566 indices are present in the database. The 
        test takes into account if more features are added to the database in the future. """
#1.)
        self.assertTrue(aaindex1.num_records() >= 566, 'Incorrect number \
            of features found in AAI1, wanted 566, got {}.'.format(aaindex1.num_records()))

    def test_feature_size(self):
        """ Test Case to check that the AAI1 has the correct dimensionality. Check that
        the number of keys in parsed JSON from the AAI1 is correct. The test takes
        into account if more features are added to the database in the future. """
#1.)
        self.assertTrue(len(aaindex1.record_codes()) >= 566, 'Incorrect number \
            of features found in AAI1, wanted 566, got {}.'.format(len(aaindex1.record_codes())))

    def test_records(self):
        """ Test Case to check that the correct record data and amino acid values 
        are returned for the accessiom number/index codes using the 
        __getitem__ function. """
#1.)
        #initialise test feature codes and their correct amino acid values
        index_code1 = 'AURR980103'
        expected_index_code1_vals = {'A': 1.05, 'L': 0.96, 'R': 0.81, 'K': 0.97, 'N': 0.91,
            'M': 0.99, 'D': 1.39, 'F': 0.95, 'C': 0.6, 'P': 1.05, 'Q': 0.87,
            'S': 0.96, 'E': 1.11, 'T': 1.03, 'G': 1.26, 'W': 1.06, 'H': 1.43,
            'Y': 0.94, 'I': 0.95, 'V': 0.62, '-': 0}

        index_code2 = 'FINA770101'
        expected_index_code2_vals = {'-': 0, 'A': 1.08, 'C': 0.95, 'D': 0.85, 'E': 1.15, 
            'F': 1.1, 'G': 0.55, 'H': 1.0, 'I': 1.05, 'K': 1.15, 'L': 1.25, 'M': 1.15, 
            'N': 0.85, 'P': 0.71, 'Q': 0.95, 'R': 1.05, 'S': 0.75, 'T': 0.75, 
            'V': 0.95, 'W': 1.1, 'Y': 1.1}

        index_code3 = 'NAGK730101'
        expected_index_code3_vals = {'A': 1.29, 'L': 1.23, 'R': 0.83, 'K': 1.23, 'N': 0.77,
            'M': 1.23, 'D': 1.0, 'F': 1.23, 'C': 0.94, 'P': 0.7, 'Q': 1.1,
            'S': 0.78, 'E': 1.54, 'T': 0.87, 'G': 0.72, 'W': 1.06, 'H': 1.29,
            'Y': 0.63, 'I': 0.94, 'V': 0.97, '-': 0}

        index_code4 = 'TANS770103'
        expected_index_code4_vals = {'-': 0, 'A': 0.79, 'C': 1.268, 'D': 0.53, 'E': 0.643, 
            'F': 1.052, 'G': 0.725, 'H': 0.864, 'I': 1.361, 'K': 0.735, 'L': 1.111, 
            'M': 1.092, 'N': 0.832, 'P': 1.249, 'Q': 1.038, 'R': 1.087, 'S': 1.093, 
            'T': 1.214, 'V': 1.428, 'W': 1.114, 'Y': 1.34}

        index_code5 = "ABCDEFGH"
        index_code6 = "123456"
        index_code7 = "blahblahblah"
#1.)
        #get amino acid values for inputted feature/index codes
        record = aaindex1[index_code1]

        self.assertIn('description', list(record.keys()), 
            "description key not found in aaindex record object.")
        self.assertIn('references', list(record.keys()), 
            "reference key not found in aaindex record object.")
        self.assertIn('notes', list(record.keys()), 
            "notes key not found in aaindex record object.")
        self.assertIn('values', list(record.keys()), 
            "values key not found in aaindex record object.")
        self.assertIn('category', list(record.keys()), 
            "category key not found in aaindex record object.")
        self.assertIn('pmid', list(record.keys()), 
            "pmid key not found in aaindex record object.")
        self.assertIn('correlation_coefficients', list(record.keys()), 
            "correlation_coefficients key not found in aaindex record object.")

        self.assertEqual(record['description'], "Normalized positional residue frequency at helix termini N' (Aurora-Rose, 1998)",
            "Description from AA1 record do not match the expected value.")
        self.assertEqual(record['notes'], '',
            "Notes from AA1 record do not match the expected value.")
        self.assertTrue(record['references'].startswith("Aurora, R. and Rose, G."),
            "Reference from AA1 record do not match the expected value.")
        self.assertEqual(record['category'], 'sec_struct', 
            "Category from AA1 record do not match the expected value.")
        self.assertEqual(record['pmid'], '9514257', 
            "PMID from AA1 record do not match the expected value.")
        self.assertEqual(record['correlation_coefficients'], {}, 
            "Correlation Coefficients from AA1 record do not match the expected value.")
        self.assertEqual(record['values'], expected_index_code1_vals, 
            "Values for individual amino acids do not match expected.")
#2.)
        record = aaindex1[index_code2]

        self.assertIn('description', list(record.keys()), 
            "description key not found in aaindex record object.")
        self.assertIn('references', list(record.keys()), 
            "reference key not found in aaindex record object.")
        self.assertIn('notes', list(record.keys()), 
            "notes key not found in aaindex record object.")
        self.assertIn('values', list(record.keys()), 
            "values key not found in aaindex record object.")
        self.assertIn('category', list(record.keys()), 
            "category key not found in aaindex record object.")
        self.assertIn('pmid', list(record.keys()), 
            "pmid key not found in aaindex record object.")
        self.assertIn('correlation_coefficients', list(record.keys()), 
            "correlation_coefficients key not found in aaindex record object.")

        self.assertEqual(record['description'], 'Helix-coil equilibrium constant (Finkelstein-Ptitsyn, 1977)',
            "Description from AA1 record do not match the expected value.")
        self.assertEqual(record['notes'], '', 
            "Notes from AA1 record do not match the expected value.")
        self.assertTrue(record['references'].startswith("Finkelstein, A.V. and Ptitsyn"), 
            "Reference from AA1 record do not match the expected value.")
        self.assertEqual(record['category'], 'sec_struct', 
            "Category from AA1 record do not match the expected value.")
        self.assertEqual(record['pmid'], '843599', 
            "PMID from AA1 record do not match the expected value.")
        self.assertEqual(record['correlation_coefficients'], {'AURR980109': '0.802', 
            'AURR980113': '0.849', 'AURR980114': '0.875', 'KANM800103': '0.823', 
            'MAXF760101': '0.810', 'PTIO830101': '0.826', 'QIAN880106': '0.810', 
            'QIAN880107': '0.814', 'SUEM840101': '0.883'}, 
            "Correlation Coefficients from AA1 record do not match the expected value.")
        self.assertEqual(record['values'], expected_index_code2_vals, 
            "Values for individual amino acids do not match expected.")
#3.)
        record = aaindex1[index_code3]

        self.assertIn('description', list(record.keys()), 
            "description key not found in aaindex record object.")
        self.assertIn('references', list(record.keys()), 
            "reference key not found in aaindex record object.")
        self.assertIn('notes', list(record.keys()), 
            "notes key not found in aaindex record object.")
        self.assertIn('values', list(record.keys()), 
            "values key not found in aaindex record object.")
        self.assertIn('category', list(record.keys()), 
            "category key not found in aaindex record object.")
        self.assertIn('pmid', list(record.keys()), 
            "pmid key not found in aaindex record object.")
        self.assertIn('correlation_coefficients', list(record.keys()), 
            "correlation_coefficients key not found in aaindex record object.")

        self.assertEqual(record['description'], 'Normalized frequency of alpha-helix (Nagano, 1973)',
            "Description from AA1 record do not match the expected value.")
        self.assertEqual(record['notes'], '',
            "Notes from AA1 record do not match the expected value.")
        self.assertTrue(record['references'].startswith("Nagano, K."),
            "Reference from AA1 record do not match the expected value.")
        self.assertEqual(record['category'], 'sec_struct', 
            "Category from AA1 record do not match the expected value.")
        self.assertEqual(record['pmid'], '4728695', 
            "PMID from AA1 record do not match the expected value.")
        self.assertEqual(record['correlation_coefficients'], 
            {'BURA740101': '0.883', 'CHOP780201': '0.886', 'CRAJ730101': '0.925', 
            'GEIM800101': '0.912', 'GEIM800104': '0.828', 'ISOY800101': '0.862', 
            'KANM800101': '0.883', 'LEVM780101': '0.894', 'LEVM780104': '0.918', 
            'MAXF760101': '0.877', 'NAGK730103': '-0.870', 'PALJ810101': '0.953', 
            'PALJ810102': '0.876', 'PRAM900102': '0.894', 'RACS820108': '0.820', 
            'ROBB760101': '0.910', 'TANS770101': '0.925'}, 
            "Correlation Coefficients from AA1 record do not match the expected value.")
        self.assertEqual(record['values'], expected_index_code3_vals, 
            "Values for individual amino acids do not match expected.")
#4.)
        record = aaindex1[index_code4]

        self.assertIn('description', list(record.keys()), 
            "description key not found in aaindex record object.")
        self.assertIn('references', list(record.keys()), 
            "reference key not found in aaindex record object.")
        self.assertIn('notes', list(record.keys()), 
            "notes key not found in aaindex record object.")
        self.assertIn('values', list(record.keys()), 
            "values key not found in aaindex record object.")
        self.assertIn('category', list(record.keys()), 
            "category key not found in aaindex record object.")
        self.assertIn('pmid', list(record.keys()), 
            "pmid key not found in aaindex record object.")
        self.assertIn('correlation_coefficients', list(record.keys()), 
            "correlation_coefficients key not found in aaindex record object.")

        self.assertEqual(record['description'], 'Normalized frequency of extended structure (Tanaka-Scheraga, 1977)',
            "Description from AA1 record do not match the expected value.")
        self.assertEqual(record['notes'], '', 
            "Notes from AA1 record do not match the expected value.")
        self.assertTrue(record['references'].startswith("Tanaka, S. and Scheraga, H.A."),
            "Reference from AA1 record do not match the expected value.")
        self.assertEqual(record['category'], 'sec_struct', 
            "Category from AA1 record do not match the expected value.")
        self.assertEqual(record['pmid'], '557155', 
            "PMID from AA1 record do not match the expected value.")
        self.assertEqual(record['correlation_coefficients'], 
            {'GEIM800105': '0.850', 'ISOY800102': '0.929', 'MAXF760102': '0.891', 
            'PALJ810103': '0.824', 'RACS820111': '0.841', 'ROBB760105': '0.871', 
            'WOEC730101': '-0.806'}, 
            "Correlation Coefficients from AA1 record do not match the expected value.")
        self.assertEqual(record['values'], expected_index_code4_vals, 
            "Values for individual amino acids do not match expected.")
#5.)
        #testing value error raised when erroneous indices sought from object
        with self.assertRaises(ValueError):
            feature_vals = aaindex1[index_code5]

        with self.assertRaises(ValueError):
            feature_vals = aaindex1[index_code6]

        with self.assertRaises(ValueError):
            feature_vals = aaindex1[index_code7]

    def test_records_dot_notation(self):
            """ Test Case to check that the correct record data and amino acid values 
            are returned for the accessiom number/index codes using the 
            __getitem__ function, accessing the elements using the dot notation created
            from the Map class. """
    #1.)
            #initialise test feature codes and their correct amino acid values
            index_code1 = 'KARP850103'
            expected_index_code1_vals = {'-': 0, 'A': 0.892, 'C': 0.925, 'D': 0.932, 
                'E': 0.933, 'F': 0.914, 'G': 0.923, 'H': 0.894, 'I': 0.872, 'K': 1.057, 
                'L': 0.921, 'M': 0.804, 'N': 0.93, 'P': 0.932, 'Q': 0.885, 'R': 0.901, 
                'S': 0.923, 'T': 0.934, 'V': 0.913, 'W': 0.803, 'Y': 0.837}

            index_code2 = 'PALJ810109'
            expected_index_code2_vals = {'-': 0, 'A': 1.15, 'C': 1.03, 'D': 1.0, 'E': 1.37, 
                'F': 0.92, 'G': 0.64, 'H': 0.95, 'I': 0.99, 'K': 1.2, 'L': 1.22, 'M': 1.45, 
                'N': 0.87, 'P': 0.72, 'Q': 1.43, 'R': 1.06, 'S': 0.84, 'T': 0.97, 'V': 0.82, 
                'W': 1.11, 'Y': 0.72}

            index_code3 = 'ROBB760112'
            expected_index_code3_vals = {'-': 0, 'A': -2.5, 'C': -4.7, 'D': 0.0, 'E': -4.4, 
                'F': -4.1, 'G': 4.9, 'H': 1.6, 'I': -3.3, 'K': -0.8, 'L': -2.0, 'M': -4.1, 
                'N': 4.6, 'P': 5.8, 'Q': -0.5, 'R': -1.2, 'S': 2.5, 'T': 1.7, 'V': -3.5, 
                'W': 1.2, 'Y': -0.6}

            index_code4 = 'TANS770103'
            expected_index_code4_vals = {'-': 0, 'A': 0.79, 'C': 1.268, 'D': 0.53, 'E': 0.643, 
                'F': 1.052, 'G': 0.725, 'H': 0.864, 'I': 1.361, 'K': 0.735, 'L': 1.111, 
                'M': 1.092, 'N': 0.832, 'P': 1.249, 'Q': 1.038, 'R': 1.087, 'S': 1.093, 
                'T': 1.214, 'V': 1.428, 'W': 1.114, 'Y': 1.34}

            index_code5 = "ABCDEFGH"
            index_code6 = "123456"
            index_code7 = "blahblahblah"
    #1.)
            #get amino acid values for inputted feature/index codes
            record = aaindex1[index_code1]

            self.assertIn('description', list(record.keys()), 
                "description key not found in aaindex record object.")
            self.assertIn('references', list(record.keys()), 
                "reference key not found in aaindex record object.")
            self.assertIn('notes', list(record.keys()), 
                "notes key not found in aaindex record object.")
            self.assertIn('values', list(record.keys()), 
                "values key not found in aaindex record object.")
            self.assertIn('category', list(record.keys()), 
                "category key not found in aaindex record object.")
            self.assertIn('pmid', list(record.keys()), 
                "pmid key not found in aaindex record object.")
            self.assertIn('correlation_coefficients', list(record.keys()), 
                "correlation_coefficients key not found in aaindex record object.")

            self.assertEqual(record.description, "Flexibility parameter for two rigid neighbors (Karplus-Schulz, 1985)",
                "Description from AA1 record do not match the expected value.")
            self.assertEqual(record.notes, '',
                "Notes from AA1 record do not match the expected value.")
            self.assertTrue(record.references.startswith("Karplus, P.A. and Schulz, G.E"),
                "Reference from AA1 record do not match the expected value.")
            self.assertEqual(record.category, 'flexibility', 
                "Category from AA1 record do not match the expected value.")
            self.assertEqual(record.pmid, '', 
                "PMID from AA1 record do not match the expected value.")
            self.assertEqual(record.correlation_coefficients, {}, 
                "Correlation Coefficients from AA1 record do not match the expected value.")
            self.assertEqual(record.values, expected_index_code1_vals, 
                "Values for individual amino acids do not match expected.")
    #2.)
            record = aaindex1[index_code2]

            self.assertIn('description', list(record.keys()), 
                "description key not found in aaindex record object.")
            self.assertIn('references', list(record.keys()), 
                "reference key not found in aaindex record object.")
            self.assertIn('notes', list(record.keys()), 
                "notes key not found in aaindex record object.")
            self.assertIn('values', list(record.keys()), 
                "values key not found in aaindex record object.")
            self.assertIn('category', list(record.keys()), 
                "category key not found in aaindex record object.")
            self.assertIn('pmid', list(record.keys()), 
                "pmid key not found in aaindex record object.")
            self.assertIn('correlation_coefficients', list(record.keys()), 
                "correlation_coefficients key not found in aaindex record object.")

            self.assertEqual(record.description, 'Normalized frequency of alpha-helix in alpha/beta class (Palau et al., 1981)',
                "Description from AA1 record do not match the expected value.")
            self.assertEqual(record.notes, '', 
                "Notes from AA1 record do not match the expected value.")
            self.assertTrue(record.references.startswith("Palau, J., Argos, P. and Puigdomenech"), 
                "Reference from AA1 record do not match the expected value.")
            self.assertEqual(record.category, 'sec_struct', 
                "Category from AA1 record do not match the expected value.")
            self.assertEqual(record.pmid, '7118409', 
                "PMID from AA1 record do not match the expected value.")
            self.assertEqual(record.correlation_coefficients, 
                {'AURR980112': '0.817', 'CHOP780201': '0.814', 
                'CRAJ730101': '0.811', 'GEIM800101': '0.816', 'GEIM800104': '0.937', 'ISOY800101': '0.874', 
                'KANM800101': '0.849', 'LEVM780101': '0.898', 'LEVM780104': '0.819', 'MAXF760101': '0.876', 
                'PALJ810102': '0.864', 'PRAM900102': '0.898', 'ROBB760101': '0.805'}, 
                "Correlation Coefficients from AA1 record do not match the expected value.")
            self.assertEqual(record.values, expected_index_code2_vals, 
                "Values for individual amino acids do not match expected.")
    #3.)
            record = aaindex1[index_code3]

            self.assertIn('description', list(record.keys()), 
                "description key not found in aaindex record object.")
            self.assertIn('references', list(record.keys()), 
                "reference key not found in aaindex record object.")
            self.assertIn('notes', list(record.keys()), 
                "notes key not found in aaindex record object.")
            self.assertIn('values', list(record.keys()), 
                "values key not found in aaindex record object.")
            self.assertIn('category', list(record.keys()), 
                "category key not found in aaindex record object.")
            self.assertIn('pmid', list(record.keys()), 
                "pmid key not found in aaindex record object.")
            self.assertIn('correlation_coefficients', list(record.keys()), 
                "correlation_coefficients key not found in aaindex record object.")

            self.assertEqual(record.description, 'Information measure for coil (Robson-Suzuki, 1976)',
                "Description from AA1 record do not match the expected value.")
            self.assertEqual(record.notes, '',
                "Notes from AA1 record do not match the expected value.")
            self.assertTrue(record.references.startswith("Robson, B. and Suzuki, E"),
                "Reference from AA1 record do not match the expected value.")
            self.assertEqual(record.category, 'sec_struct', 
                "Category from AA1 record do not match the expected value.")
            self.assertEqual(record.pmid, '1003471', 
                "PMID from AA1 record do not match the expected value.")
            self.assertEqual(record.correlation_coefficients, 
                {'CHOP780211': '0.841', 'ISOY800103': '0.807', 'PALJ810115': '0.885', 
                'QIAN880132': '0.800', 'QIAN880133': '0.814'}, 
                "Correlation Coefficients from AA1 record do not match the expected value.")
            self.assertEqual(record.values, expected_index_code3_vals, 
                "Values for individual amino acids do not match expected.")
    #4.)
            record = aaindex1[index_code4]

            self.assertIn('description', list(record.keys()), 
                "description key not found in aaindex record object.")
            self.assertIn('references', list(record.keys()), 
                "reference key not found in aaindex record object.")
            self.assertIn('notes', list(record.keys()), 
                "notes key not found in aaindex record object.")
            self.assertIn('values', list(record.keys()), 
                "values key not found in aaindex record object.")
            self.assertIn('category', list(record.keys()), 
                "category key not found in aaindex record object.")
            self.assertIn('pmid', list(record.keys()), 
                "pmid key not found in aaindex record object.")
            self.assertIn('correlation_coefficients', list(record.keys()), 
                "correlation_coefficients key not found in aaindex record object.")

            self.assertEqual(record.description, 'Normalized frequency of extended structure (Tanaka-Scheraga, 1977)',
                "Description from AA1 record do not match the expected value.")
            self.assertEqual(record.notes, '', 
                "Notes from AA1 record do not match the expected value.")
            self.assertTrue(record.references.startswith("Tanaka, S. and Scheraga, H.A."),
                "Reference from AA1 record do not match the expected value.")
            self.assertEqual(record.category, 'sec_struct', 
                "Category from AA1 record do not match the expected value.")
            self.assertEqual(record.pmid, '557155', 
                "PMID from AA1 record do not match the expected value.")
            self.assertEqual(record.correlation_coefficients, 
                {'GEIM800105': '0.850', 'ISOY800102': '0.929', 'MAXF760102': '0.891', 'PALJ810103': '0.824', 
                'RACS820111': '0.841', 'ROBB760105': '0.871', 'WOEC730101': '-0.806'}, 
                "Correlation Coefficients from AA1 record do not match the expected value.")
            self.assertEqual(record.values, expected_index_code4_vals, 
                "Values for individual amino acids do not match expected.")
    #5.)
            #testing value error raised when erroneous indices sought from object
            with self.assertRaises(ValueError):
                feature_vals = aaindex1[index_code5]

            with self.assertRaises(ValueError):
                feature_vals = aaindex1[index_code6]

            with self.assertRaises(ValueError):
                feature_vals = aaindex1[index_code7]

    def test_search(self):
        """ Testing search functionality where 1 or more aaindex1 records can be returned 
        if their description contains user-inputted keywords. """

        description1 = "Gibbs Energy"
        description2 = 'Polarity'
        description3 = 'Hydrophobic'
        description4 = 'pi-helices'
        description5 = 'blahblahblah'
        description6 = 'not a description'
        description7 = 1234
        description8 = False
#1.)
        description1_search = aaindex1.search(description1)
        expected_description1_indices = ['YUTK870101', 'YUTK870102', 'YUTK870103', 'YUTK870104']

        self.assertEqual(len(description1_search), 4, 
            "Expected there to be 4 records returned from search, got {}.".format(len(description1_search)))
        for index, val in description1_search.items():
            self.assertIn(index, expected_description1_indices, 
                "Index record {} expected to be found in search results.".format(index))
            self.assertTrue(description1.lower() in description1_search[index]['description'].lower(), 
                "Expected description to be in record of returned search results: {}".format(description1_search[index]))
#2.)
        description2_search = aaindex1.search(description2)
        expected_description2_indices = ['GRAR740102', 'RADA880108', 'ZIMJ680103']
      
        self.assertEqual(len(description2_search), 3, 
            "Expected there to be 3 records returned from search, got {}.".format(len(description2_search)))
        for index, val in description2_search.items():
            self.assertIn(index, expected_description2_indices,
                "Index record {} expected to be found in search results.".format(index))
            self.assertTrue(description2.lower() in description2_search[index]['description'].lower(),
                "Expected description to be in record of returned search results: {}.".format(description2_search[index]))
#3.)
        description3_search = aaindex1.search(description3)
        expected_description3_indices = ['ARGP820101', 'BLAS910101', 'CASG920101', 'CIDH920101', 'CIDH920102', 'CIDH920103', 
        'CIDH920104', 'CIDH920105', 'COWR900101', 'EISD840101', 'EISD860102', 'EISD860103', 'ENGD860101', 'FASG890101', 
        'FAUJ830101', 'GOLD730101', 'JOND750101', 'JURD980101', 'KIDA850101', 'LEVM760101', 'MANP780101', 'PONP800101', 
        'PONP800102', 'PONP800103', 'PONP800104', 'PONP800105', 'PONP800106', 'PONP930101', 'PRAM900101', 'SWER830101', 
        'WILM950101', 'WILM950102', 'WILM950103', 'WILM950104', 'WOLR790101', 'ZIMJ680101']
 
        self.assertEqual(len(description3_search), 36, 
            "Expected there to be 36 records returned from search, got {}.".format(len(description3_search)))
        for index, val in description3_search.items():
            self.assertIn(index, expected_description3_indices, 
                "Index record {} expected to be found in search results.".format(index))
            self.assertTrue(description3.lower() in description3_search[index]['description'].lower(),
                "Expected description to be in record of returned search results: {}.".format(description3_search[index]))
#4.)
        description4_search = aaindex1.search(description4)
        expected_description4_indices = ['FODM020101']

        self.assertEqual(len(description4_search), 1, 
            "Expected there to be 1 records returned from search, got {}.".format(len(description4_search)))
        for index, val in description4_search.items():
            self.assertIn(index, expected_description4_indices, 
                "Index record {} expected to be found in search results.".format(index))
            self.assertTrue(description4.lower() in description4_search[index]['description'].lower(),
                "Expected description to be in record of returned search results.")
#5.)
        description5_search = aaindex1.search(description5)
        self.assertEqual(len(description5_search), 0, 
            "Expected there to be 0 records returned from search, got {}.".format(len(description5_search)))
#6.)
        description6_search = aaindex1.search(description6)
        self.assertEqual(len(description6_search), 0, 
            "Expected there to be 0 records returned from search, got {}.".format(len(description6_search)))
#7.)    
        with (self.assertRaises(TypeError)):
            description7_search = aaindex1.search(description7)
#8.)
        with (self.assertRaises(TypeError)):
            description8_search = aaindex1.search(description8)

    def test_record_codes(self):
        """ Test Case to check that correct feature/record codes are in the parsed JSON
        of all AAI1 records. Also testing that closely named erroneous records codes are
        not present. """

        #testing actual index codes
        index_code1 = 'VHEG790101'
        index_code2 = 'PONP800107'
        index_code3 = 'OOBM770102'
        index_code4 = 'NADH010101'
#1.)
        #testing index codes that are in the AAI1
        self.assertIn(index_code1, aaindex1.record_codes(), 'Index {} \
            not found in list of available indices.'.format(index_code1))
        self.assertIn(index_code2, aaindex1.record_codes(), 'Index {} \
            not found in list of available indices'.format(index_code2))
        self.assertIn(index_code3, aaindex1.record_codes(), 'Index {} \
            not found in list of available indices'.format(index_code3))
        self.assertIn(index_code4, aaindex1.record_codes(), 'Index {} \
            not found in list of available indices'.format(index_code4))
#2.)
        #testing closely named index codes
        index_code5 = 'VHEG790500'
        index_code6 = 'P0NP800107'
        index_code7 = 'OObm770102'
        index_code8 = 'NADW210101'

        #testing errenous index codes are not in the AAI1
        self.assertNotIn(index_code5, aaindex1.record_codes(), 'Index {} \
            erroneously found in list of available indices.'.format(index_code5))
        self.assertNotIn(index_code6, aaindex1.record_codes(), 'Index {} \
            erroneously found in list of available indices.'.format(index_code6))
        self.assertNotIn(index_code7, aaindex1.record_codes(), 'Index {} \
            erroneously found in list of available indices.'.format(index_code7))
        self.assertNotIn(index_code8, aaindex1.record_codes(), 'Index {} \
            erroneously found in list of available indices.'.format(index_code8))

    def test_aaindex_names(self):
        """ Test Case to check that the correct names/description is returned for each index. """
        
        record_names = aaindex1.record_names()
        record_name1 = record_names[50]
        record_name2 = record_names[140]
        record_name3 = record_names[368]
        record_name4 = record_names[560]
#1.)
        self.assertEqual(len(record_names), aaindex1.num_records(), 
            'Total number of names in the AAI1 database should equal the number of records.')
#2.)
        self.assertEqual(record_name1, 'alpha-NH chemical shifts (Bundi-Wuthrich, 1979)',
            'Incorrect index name returned, got {} instead of {}'.format(record_name1,
                'alpha-NH chemical shifts (Bundi-Wuthrich, 1979).'))
#3.)
        self.assertTrue(record_name2.startswith('Number of full nonbonding'),
            'Incorrect index name returned, got {} instead of {}.'.format(record_name2,
                'Number of full nonbonding'))
#4.)
        self.assertIn('alpha-helix in alpha/beta class', record_name3, 
            'String {} not found in index name: {}'.format(
                'alpha-helix in alpha/beta class', record_name3))
#5.)
        self.assertTrue(record_name4.endswith('Buriability (Zhou-Zhou, 2004)'), 
            'Index does not end with {}, got {}'.format(
                'Buriability (Zhou-Zhou, 2004)', record_name4))

    def test_aaindex_refs(self):
        """ Test Case to check that the correct references are returned for each index. """

        feature1 = 'VELV850101'
        feature2 = 'QIAN880139'
        feature3 = 'NAKH900112'
        feature4 = 'LEVM760103'
        feature5 = "blah_blah"
        feature6 = 123
#1.)
        ref1 = aaindex1[feature1]['references']
        ref2 = aaindex1[feature2]['references']
        ref3 = aaindex1[feature3]['references']
        ref4 = aaindex1[feature4]['references']

        self.assertTrue(ref1.startswith('Veljkovic'), 'Incorrect index reference, \
            got {} instead of {} '.format(ref1,'' ))
        self.assertIn("globular proteins using neural network models", ref2)
        self.assertIn("hydrophobicity of amino acid composition of mitochondrial proteins", ref3)
        self.assertTrue(ref4.endswith('J. Mol. Biol. 104, 59-107 (1976) (Gly missing)'), "{}".format(ref4))
#2.)
        with (self.assertRaises(ValueError)):
            ref5 = aaindex1[feature5]['references']

        with (self.assertRaises(TypeError)):
            ref6 = aaindex1[feature6]['references']

    def test_last_updates(self):
        """ Testing the last updated class attribute which is the latest version of the database. """
        self.assertEqual(aaindex1.last_updated, "February 13, 2017", 
            "Last updated value does not match expected, got {}.".format(aaindex1.last_updated))

    def tearDown(self):
        """ Remove any test data or directories. """
        os.rmdir('test_data')

if __name__ == '__main__':
    #run all unit tests
    unittest.main(verbosity=2)