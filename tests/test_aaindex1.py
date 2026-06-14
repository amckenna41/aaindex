################################################################################
################             AAindex1 Module Tests             #################
################################################################################

import unittest
from unittest.mock import patch
from importlib.metadata import metadata
from aaindex import aaindex1, __version__

class AAIndex1_Tests(unittest.TestCase):
    """
    Test suite for testing the aaindex1 module in the aaindex Python software package.

    Test Cases
    ==========
    test_aaindex_metadata:
        testing correct aaindex version and metadata.
    test_num_records:
        testing to check that the correct number of indices/records are present
        in the AAindex object. 
    test_records:
        testing to check that the correct record data and amino acid values 
        are returned for the accessiom number/index codes using the __getitem__ 
        function.
    test_records_dot_notation:
        testing to check that the correct record data and amino acid values 
        are returned for the accessiom number/index codes using the __getitem__ 
        function, accessing the elements using the dot notation created from 
        the Map class.
    test_search:
        testing search functionality where 1 or more aaindex1 records can be returned 
        if their description contains user-inputted keywords.
    test_record_codes:
        testing to check that correct feature/record codes are in the parsed JSON
        of all AAi1 records.
    test_record_names:
        testing to check that the correct names/description is returned for each index.
    test_last_updated:
        testing hard-coded date value that the AAindex database were last updated.
    test_values:
        testing the values() method returns the correct amino acid values dict for a given record.
    test_get_record_by_category:
        testing get_record_by_category() returns all records belonging to a given category.
    test_dunder_methods:
        testing __len__, __contains__, __iter__, and __repr__ dunder methods.
    test_repr_format:
        testing __repr__ returns the expected format string.
    """
    def test_aaindex_metadata(self):
        """ Testing correct aaindex version and metadata. """
        self.assertEqual(__version__, "1.3.0",
            f"aaindex version is not correct, got: {metadata('aaindex')['version']}.")
        self.assertEqual(metadata("aaindex")['name'], "aaindex", 
            f"aaindex software name is not correct, got: {metadata('aaindex')['name']}.")
        self.assertEqual(metadata("aaindex")['author-email'], "AJ McKenna <amckenna41@qub.ac.uk>", 
            f"aaindex author-email is not correct, got: {metadata('aaindex')['author-email']}.")
        self.assertEqual(metadata('aaindex')['maintainer'], "AJ McKenna", 
            f"aaindex maintainer is not correct, got: {metadata('aaindex')['maintainer']}.")
        self.assertEqual(metadata('aaindex')['License-Expression'], "MIT", 
            f"aaindex license type is not correct, got: {metadata('aaindex')['License-Expression']}.")
        self.assertEqual(metadata('aaindex')['keywords'], 
            "amino acid index,aaindex,bioinformatics,protein engineering,python,pypi,"
            "physicochemical properties,biochemical properties,proteins,protein structure prediction,pysar", 
                f"aaindex keywords are not correct, got:\n{metadata('aaindex')['keywords']}.")

    def test_amino_acids(self):
        """ Test Case to check that only valid Amino Acids are used within the AAindex class. """
        valid_amino_acids = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
            'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        amino_acids = aaindex1.amino_acids()
#1.)
        for aa in amino_acids:
            self.assertIn(aa, valid_amino_acids, f'Amino Acid {aa} not in list of valid amino acids:\n{valid_amino_acids}.')

    def test_num_records(self):
        """ Test Case to check that the correct number of indices/records are present
        in the AAindex object. To date, 566 indices are present in the database. """
#1.)
        self.assertEqual(aaindex1.num_records(), 566, 
            f'Incorrect number of records found in AAi1, expected 566, got {aaindex1.num_records()}.')

    def test_records(self):
        """ Test Case to check that the correct record data and amino acid values 
        are returned for the accessiom number/index codes using the __getitem__ function. """
#1.)
        #initialise test records codes and their correct and expected amino acid values
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

        index_code3 = 'nagk730101' #function can accept lowercase record codes
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
        index_code8 = ""
#1.)
        #get amino acid values for inputted feature/index codes
        record = aaindex1[index_code1] #AURR980103

        self.assertIn('description', list(record.keys()), "description key not found in aaindex record object.")
        self.assertIn('references', list(record.keys()), "references key not found in aaindex record object.")
        self.assertIn('notes', list(record.keys()), "notes key not found in aaindex record object.")
        self.assertIn('values', list(record.keys()), "values key not found in aaindex record object.")
        self.assertIn('category', list(record.keys()), "category key not found in aaindex record object.")
        self.assertIn('pmid', list(record.keys()), "pmid key not found in aaindex record object.")
        self.assertIn('correlation_coefficients', list(record.keys()), "correlation_coefficients key not found in aaindex record object.")

        self.assertEqual(record['description'], "Normalized positional residue frequency at helix termini N' (Aurora-Rose, 1998)",
            f"Description from AAi1 record does not match the expected value, got {record['description']}.")
        self.assertEqual(record['notes'], '', 
            f"Notes from AAi1 record does not match the expected value, got {record['notes']}.")
        self.assertEqual(record['references'], "Aurora, R. and Rose, G. 'Helix capping' Protein Science 7, 21-38 (1998)", 
            f"Reference from AAi1 record do not match the expected value, got {record['references']}.")
        self.assertEqual(record['category'], 'sec_struct', 
            f"Category from AAi1 record does not match the expected value, got {record['category']}.")
        self.assertEqual(record['pmid'], '9514257', 
            f"PMID from AAi1 record does not match the expected value, got {record['pmid']}.")
        self.assertEqual(record['correlation_coefficients'], {}, 
            f"Correlation Coefficients from AAi1 record do not match the expected value, got {record['correlation_coefficients']}.")
        self.assertEqual(record['values'], expected_index_code1_vals, 
            f"Values for individual amino acids do not match expected, got\n{record['values']}.")
        for aa, val in list(record['values'].items()):
            if (aa == '-'):
                continue
            self.assertIsInstance(val, (float, type(None)),
                f"Values in record values not of correct datatype, expected float or None, got {type(val)}.")

        self.assertIn('description', list(record.keys()), "description key not found in aaindex record object.")
        self.assertIn('references', list(record.keys()), "references key not found in aaindex record object.")
        self.assertIn('notes', list(record.keys()), "notes key not found in aaindex record object.")
        self.assertIn('values', list(record.keys()), "values key not found in aaindex record object.")
        self.assertIn('category', list(record.keys()), "category key not found in aaindex record object.")
        self.assertIn('pmid', list(record.keys()), "pmid key not found in aaindex record object.")
        self.assertIn('correlation_coefficients', list(record.keys()), "correlation_coefficients key not found in aaindex record object.")
#2.)
        record = aaindex1[index_code2] #FINA770101

        self.assertEqual(record['description'], 'Helix-coil equilibrium constant (Finkelstein-Ptitsyn, 1977)',
            f"Description from AAi1 record does not match the expected value, got {record['description']}.")
        self.assertEqual(record['notes'], '',
            f"Notes from AAi1 record does not match the expected value, got {record['notes']}.")
        self.assertEqual(record['references'], "Finkelstein, A.V. and Ptitsyn, O.B. 'Theory of protein molecule self-organization. "
                          "II. A comparison of calculated thermodynamic parameters of local secondary structures with experiments' "
                          "Biopolymers 16, 497-524 (1977) (Pro 0.096)", 
            f"Reference from AAi1 record do not match the expected value, got {record['references']}.")
        self.assertEqual(record['category'], 'sec_struct', 
            f"Category from AAi1 record does not match the expected value, got {record['category']}.")
        self.assertEqual(record['pmid'], '843599', 
            f"PMID from AAi1 record does not match the expected value, got {record['pmid']}.")
        self.assertEqual(record['correlation_coefficients'], {'AURR980109': 0.802,
            'AURR980113': 0.849, 'AURR980114': 0.875, 'KANM800103': 0.823,
            'MAXF760101': 0.81, 'PTIO830101': 0.826, 'QIAN880106': 0.81,
            'QIAN880107': 0.814, 'SUEM840101': 0.883},
            f"Correlation Coefficients from AAi1 record do not match the expected value, got {record['correlation_coefficients']}.")
        self.assertEqual(record['values'], expected_index_code2_vals, 
            f"Values for individual amino acids do not match expected, got\n{record['values']}.")
        for aa, val in list(record['values'].items()):
            if (aa == '-'):
                continue
            self.assertIsInstance(val, (float, type(None)),
                f"Values in record values not of correct datatype, expected float or None, got {type(val)}.")
#3.)
        record = aaindex1[index_code3] #nagk730101

        self.assertIn('description', list(record.keys()), "description key not found in aaindex record object.")
        self.assertIn('references', list(record.keys()), "references key not found in aaindex record object.")
        self.assertIn('notes', list(record.keys()), "notes key not found in aaindex record object.")
        self.assertIn('values', list(record.keys()), "values key not found in aaindex record object.")
        self.assertIn('category', list(record.keys()), "category key not found in aaindex record object.")
        self.assertIn('pmid', list(record.keys()), "pmid key not found in aaindex record object.")
        self.assertIn('correlation_coefficients', list(record.keys()), "correlation_coefficients key not found in aaindex record object.")
        
        self.assertEqual(record['description'], 'Normalized frequency of alpha-helix (Nagano, 1973)',
            f"Description from AAi1 record does not match the expected value, got {record['description']}.")
        self.assertEqual(record['notes'], '', 
            f"Notes from AAi1 record does not match the expected value, got {record['notes']}.")
        self.assertEqual(record['references'], "Nagano, K. 'Local analysis of the mechanism of protein folding. "
                        "I. Prediction of helices, loops, and beta-structures from primary structure' J. Mol. Biol. 75, 401-420 (1973)",
            f"Reference from AAi1 record do not match the expected value, got {record['references']}.")
        self.assertEqual(record['category'], 'sec_struct',
            f"Category from AAi1 record does not match the expected value, got {record['category']}.")
        self.assertEqual(record['pmid'], '4728695', 
            f"PMID from AAi1 record does not match the expected value, got {record['pmid']}.")
        self.assertEqual(record['correlation_coefficients'],
            {'BURA740101': 0.883, 'CHOP780201': 0.886, 'CRAJ730101': 0.925,
            'GEIM800101': 0.912, 'GEIM800104': 0.828, 'ISOY800101': 0.862,
            'KANM800101': 0.883, 'LEVM780101': 0.894, 'LEVM780104': 0.918,
            'MAXF760101': 0.877, 'NAGK730103': -0.87, 'PALJ810101': 0.953,
            'PALJ810102': 0.876, 'PRAM900102': 0.894, 'RACS820108': 0.82,
            'ROBB760101': 0.91, 'TANS770101': 0.925},
            f"Correlation Coefficients from AAi1 record do not match the expected value, got {record['correlation_coefficients']}.")
        self.assertEqual(record['values'], expected_index_code3_vals, 
            f"Values for individual amino acids do not match expected, got\n{record['values']}.")
        for aa, val in list(record['values'].items()):
            if (aa == '-'):
                continue
            self.assertIsInstance(val, (float, type(None)),
                f"Values in record values not of correct datatype, expected float or None, got {type(val)}.")
#4.)
        record = aaindex1[index_code4] #TANS770103

        self.assertIn('description', list(record.keys()), "description key not found in aaindex record object.")
        self.assertIn('references', list(record.keys()), "reference key not found in aaindex record object.")
        self.assertIn('notes', list(record.keys()), "notes key not found in aaindex record object.")
        self.assertIn('values', list(record.keys()), "values key not found in aaindex record object.")
        self.assertIn('category', list(record.keys()), "category key not found in aaindex record object.")
        self.assertIn('pmid', list(record.keys()), "pmid key not found in aaindex record object.")
        self.assertIn('correlation_coefficients', list(record.keys()), "correlation_coefficients key not found in aaindex record object.")

        self.assertEqual(record['description'], 'Normalized frequency of extended structure (Tanaka-Scheraga, 1977)',
            f"Description from AAi1 record does not match the expected value, got {record['description']}.")
        self.assertEqual(record['notes'], '', 
            f"Notes from AAi1 record does not match the expected value, got {record['notes']}.")
        self.assertEqual(record['references'], "Tanaka, S. and Scheraga, H.A. 'Statistical mechanical treatment of protein"
                        " conformation. 5. A multiphasic model for specific-sequence copolymers of amino acids' Macromolecules"
                        " 10, 9-20 (1977) Recalculated by Kidera as normalized frequencies", 
            f"Reference from AAi1 record do not match the expected value, got {record['references']}.")
        self.assertEqual(record['category'], 'sec_struct', 
            f"Category from AAi1 record does not match the expected value, got {record['category']}.")
        self.assertEqual(record['pmid'], '557155', 
            f"PMID from AAi1 record does not match the expected value, got {record['pmid']}.")
        self.assertEqual(record['correlation_coefficients'],
            {'GEIM800105': 0.85, 'ISOY800102': 0.929, 'MAXF760102': 0.891,
            'PALJ810103': 0.824, 'RACS820111': 0.841, 'ROBB760105': 0.871,
            'WOEC730101': -0.806},
            f"Correlation Coefficients from AAi1 record do not match the expected value, got {record['correlation_coefficients']}.")
        self.assertEqual(record['values'], expected_index_code4_vals,
            f"Values for individual amino acids do not match expected, got\n{record['values']}.")
        for aa, val in list(record['values'].items()):
            if (aa == '-'):
                continue
            self.assertIsInstance(val, (float, type(None)),
                f"Values in record values not of correct datatype, expected float or None, got {type(val)}.")
#5.)
        #testing value error raised when erroneous indices sought from object
        with self.assertRaises(ValueError):
            aaindex1[index_code5]
        with self.assertRaises(ValueError):
            aaindex1[index_code6]
        with self.assertRaises(ValueError):
            aaindex1[index_code7]
        with self.assertRaises(ValueError):
            aaindex1[index_code8]

    def test_records_dot_notation(self):
            """ Test Case to check that the correct record data and amino acid values 
            are returned for the accessiom number/index codes using the __getitem__ 
            function, accessing the elements using the dot notation created from the 
            Map class. """
            #initialise test record codes and their correct and expected amino acid values
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
            index_code8 = ""
#1.)
            #get amino acid values for inputted record/index codes
            record = aaindex1[index_code1]

            self.assertEqual(record.description, "Flexibility parameter for two rigid neighbors (Karplus-Schulz, 1985)",
                f"Description from AAi1 record does not match the expected value, got {record['description']}.")
            self.assertEqual(record.notes, '',
                f"Notes from AAi1 record does not match the expected value, got {record['notes']}.")
            self.assertEqual(record.references, "Karplus, P.A. and Schulz, G.E. 'Prediction of chain flexibility in proteins' Naturwiss. 72, 212-213 (1985)",
                f"Reference from AAi1 record do not match the expected value, got {record['references']}.")
            self.assertEqual(record.category, 'flexibility', 
                f"Category from AAi1 record does not match the expected value, got {record['category']}.")
            self.assertEqual(record.pmid, '', 
                f"PMID from AAi1 record does not match the expected value, got {record['pmid']}.")
            self.assertEqual(record.correlation_coefficients, {}, 
                f"Correlation Coefficients from AAi1 record do not match the expected value, got {record['correlation_coefficients']}.")
            self.assertEqual(record.values, expected_index_code1_vals, 
                f"Values for individual amino acids do not match expected, got\n{record['values']}.")
#2.)
            record = aaindex1[index_code2]

            self.assertEqual(record.description, 'Normalized frequency of alpha-helix in alpha/beta class (Palau et al., 1981)',
                f"Description from AAi1 record does not match the expected value, got {record['description']}.")
            self.assertEqual(record.notes, '', 
                f"Notes from AAi1 record does not match the expected value, got {record['notes']}.")
            self.assertTrue(record.references.startswith("Palau, J., Argos, P. and Puigdomenech"), 
                f"Reference from AAi1 record do not match the expected value, got {record['references']}.")
            self.assertEqual(record.category, 'sec_struct', 
                f"Category from AAi1 record does not match the expected value, got {record['category']}.")
            self.assertEqual(record.pmid, '7118409', 
                f"PMID from AAi1 record does not match the expected value, got {record['pmid']}.")
            self.assertEqual(record.correlation_coefficients,
                {'AURR980112': 0.817, 'CHOP780201': 0.814,
                'CRAJ730101': 0.811, 'GEIM800101': 0.816, 'GEIM800104': 0.937, 'ISOY800101': 0.874,
                'KANM800101': 0.849, 'LEVM780101': 0.898, 'LEVM780104': 0.819, 'MAXF760101': 0.876,
                'PALJ810102': 0.864, 'PRAM900102': 0.898, 'ROBB760101': 0.805},
                f"Correlation Coefficients from AAi1 record do not match the expected value, got {record['correlation_coefficients']}.")
            self.assertEqual(record.values, expected_index_code2_vals, 
                f"Values for individual amino acids do not match expected, got\n{record['values']}.")
#3.)
            record = aaindex1[index_code3]

            self.assertEqual(record.description, 'Information measure for coil (Robson-Suzuki, 1976)',
                f"Description from AAi1 record does not match the expected value, got {record['description']}.")
            self.assertEqual(record.notes, '',
                f"Notes from AAi1 record does not match the expected value, got {record['notes']}.")
            self.assertTrue(record.references.startswith("Robson, B. and Suzuki, E"),
                f"Reference from AAi1 record do not match the expected value, got {record['references']}.")
            self.assertEqual(record.category, 'sec_struct', 
                f"Category from AAi1 record does not match the expected value, got {record['category']}.")
            self.assertEqual(record.pmid, '1003471', 
                f"PMID from AAi1 record does not match the expected value, got {record['pmid']}.")
            self.assertEqual(record.correlation_coefficients,
                {'CHOP780211': 0.841, 'ISOY800103': 0.807, 'PALJ810115': 0.885,
                'QIAN880132': 0.8, 'QIAN880133': 0.814},
                f"Correlation Coefficients from AAi1 record do not match the expected value, got {record['correlation_coefficients']}.")
            self.assertEqual(record.values, expected_index_code3_vals, 
                f"Values for individual amino acids do not match expected, got\n{record['values']}.")
#4.)
            record = aaindex1[index_code4]

            self.assertEqual(record.description, 'Normalized frequency of extended structure (Tanaka-Scheraga, 1977)',
                f"Description from AAi1 record does not match the expected value, got {record['description']}.")
            self.assertEqual(record.notes, '', 
                f"Notes from AAi1 record does not match the expected value, got {record['notes']}.")
            self.assertTrue(record.references.startswith("Tanaka, S. and Scheraga, H.A."),
                f"Reference from AAi1 record do not match the expected value, got {record['references']}.")
            self.assertEqual(record.category, 'sec_struct', 
                f"Category from AAi1 record does not match the expected value, got {record['category']}.")
            self.assertEqual(record.pmid, '557155', 
                f"PMID from AAi1 record does not match the expected value, got {record['pmid']}.")
            self.assertEqual(record.correlation_coefficients,
                {'GEIM800105': 0.85, 'ISOY800102': 0.929, 'MAXF760102': 0.891, 'PALJ810103': 0.824,
                'RACS820111': 0.841, 'ROBB760105': 0.871, 'WOEC730101': -0.806},
                f"Correlation Coefficients from AAi1 record do not match the expected value, got {record['correlation_coefficients']}.")
            self.assertEqual(record.values, expected_index_code4_vals, 
                f"Values for individual amino acids do not match expected, got\n{record['values']}.")
#5.)
            #testing value error raised when erroneous indices sought from object
            with self.assertRaises(ValueError):
                aaindex1[index_code5]
            with self.assertRaises(ValueError):
                aaindex1[index_code6]
            with self.assertRaises(ValueError):
                aaindex1[index_code7]
            with self.assertRaises(ValueError):
                aaindex1[index_code8]

    def test_search(self):
        """ Testing search functionality where 1 or more aaindex1 records will be returned 
            if their description contains user-inputted keywords. """ 
        description1 = "Gibbs Energy"
        description2 = 'Polarity'
        description3 = 'Hydrophobic'
        description4 = ['pi-helices', 'Retention coefficient', 'p-Values']
        description5 = 'blahblahblah'
        description6 = 'not a description'
        description7 = 1234
        description8 = 1.24
        description9 = False
#1.)
        description1_search = aaindex1.search(description1)
        expected_description1_indices = ['YUTK870101', 'YUTK870102', 'YUTK870103', 'YUTK870104']

        self.assertGreaterEqual(len(description1_search), 4, 
            f"Expected at least 4 records returned from search, got {len(description1_search)}.")
        for index in expected_description1_indices:
            self.assertIn(index, description1_search,
                f"Expected {index} to be in description1 search results.")
        for index, val in description1_search.items():
            self.assertTrue(description1.lower() in description1_search[index]['description'].lower()
                or description1.lower() in index.lower()
                or description1.lower() in description1_search[index]['pmid'].lower()
                or description1.lower() in description1_search[index]['references'].lower()
                or description1.lower() in description1_search[index]['notes'].lower(),
                f"Expected search keyword in some field of returned record {index}.")
#2.)
        description2_search = aaindex1.search(description2)
        expected_description2_indices = ['GRAR740102', 'RADA880108', 'ZIMJ680103']
      
        self.assertGreaterEqual(len(description2_search), 3, 
            f"Expected at least 3 records returned from search, got {len(description2_search)}.")
        for index in expected_description2_indices:
            self.assertIn(index, description2_search,
                f"Expected {index} to be in description2 search results.")
#3.)
        description3_search = aaindex1.search(description3)
        expected_description3_indices = ['ARGP820101', 'BLAS910101', 'CASG920101', 'CIDH920101', 'CIDH920102', 'CIDH920103', 
            'CIDH920104', 'CIDH920105', 'COWR900101', 'EISD840101', 'EISD860102', 'EISD860103', 'ENGD860101', 'FASG890101', 
            'FAUJ830101', 'GOLD730101', 'JOND750101', 'JURD980101', 'KIDA850101', 'LEVM760101', 'MANP780101', 'PONP800101', 
            'PONP800102', 'PONP800103', 'PONP800104', 'PONP800105', 'PONP800106', 'PONP930101', 'PRAM900101', 'SWER830101', 
            'WILM950101', 'WILM950102', 'WILM950103', 'WILM950104', 'WOLR790101', 'ZIMJ680101']
 
        self.assertGreaterEqual(len(description3_search), 36, 
            f"Expected at least 36 records returned from search, got {len(description3_search)}.")
        for index in expected_description3_indices:
            self.assertIn(index, description3_search,
                f"Expected {index} to be in description3 search results.")
#4.)
        description4_search = aaindex1.search(description4)
        expected_description4_indices = ['FODM020101', 'BROC820101', 'BROC820102', 'GUOD860101', 'MEEJ800101', 'MEEJ800102', \
                                         'MEEJ810101', 'MEEJ810102', 'PARS000101', 'PARS000102']

        self.assertGreaterEqual(len(description4_search), 10, 
            f"Expected at least 10 records returned from search, got {len(description4_search)}.")
        for index in expected_description4_indices:
            self.assertIn(index, description4_search,
                f"Expected {index} to be in description4 search results.")
#5.)
        description5_search = aaindex1.search(description5)
        self.assertEqual(len(description5_search), 0, 
            f"Expected there to be 0 records returned from search, got {len(description5_search)}.")
#6.)
        description6_search = aaindex1.search(description6)
        self.assertEqual(len(description6_search), 0, 
            f"Expected there to be 0 records returned from search, got {len(description6_search)}.")
#7.)    
        with (self.assertRaises(TypeError)):
            aaindex1.search(description7)
            aaindex1.search(description8)
            aaindex1.search(description9)
#8.) search by PMID — should find the record whose pmid field matches
        pmid_search = aaindex1.search('9514257')
        self.assertIn('AURR980103', pmid_search,
            'search() should find AURR980103 when searching by its PMID 9514257.')
#9.) search by exact accession code — should find the record itself
        accession_search = aaindex1.search('CHOP780206')
        self.assertIn('CHOP780206', accession_search,
            'search() should find CHOP780206 when searching by its accession code.')

    def test_na_values(self):
        """ Test that amino acids with no measured value are stored as None,
        not as the string 'NA' or the number 0. """
#1.) AVBF000101 has P=NA in the raw database
        record = aaindex1['AVBF000101']
        self.assertIsNone(record['values']['P'],
            f"Expected None for P in AVBF000101 (NA in source), got {record['values']['P']}.")
#2.) non-NA values in the same record should still be floats
        self.assertIsInstance(record['values']['A'], float,
            f"Expected float for A in AVBF000101, got {type(record['values']['A'])}.")
#3.) the '-' gap placeholder should be 0.0 (float), not 0 (int)
        self.assertEqual(record['values']['-'], 0.0,
            "Gap placeholder '-' should be 0.0 (float).")
        self.assertIsInstance(record['values']['-'], float,
            "Gap placeholder '-' should be stored as float, not int.")
#4.) arithmetic on a None value should raise TypeError, not silently return 0
        na_val = record['values']['P']
        with self.assertRaises(TypeError):
            _ = na_val + 1.0

    def test_record_codes(self):
        """ Test Case to check that correct feature/record codes are in the parsed JSON
            of all AAi1 records. Also testing that closely named erroneous records codes 
            are not present. """
        #testing actual index codes
        index_code1 = 'VHEG790101'
        index_code2 = 'PONP800107'
        index_code3 = 'OOBM770102'
        index_code4 = 'NADH010101'
#1.)
        #testing index codes that are in the AAi1
        self.assertIn(index_code1, aaindex1.record_codes(), f'Expected record {index_code1} to be in list of available indices.')
        self.assertIn(index_code2, aaindex1.record_codes(), f'Expected record {index_code2} to be in list of available indices.')
        self.assertIn(index_code3, aaindex1.record_codes(), f'Expected record {index_code3} to be in list of available indices.')
        self.assertIn(index_code4, aaindex1.record_codes(), f'Expected record {index_code4} to be in list of available indices.')
#2.)
        #testing closely named index codes
        index_code5 = 'VHEG790500'
        index_code6 = 'P0NP800107'
        index_code7 = 'OObm770102'
        index_code8 = 'NADW210101'

        #testing errenous index codes are not in the AAi1
        self.assertNotIn(index_code5, aaindex1.record_codes(), f'Expected record {index_code5} to not be in list of available indices.')
        self.assertNotIn(index_code6, aaindex1.record_codes(), f'Expected record {index_code6} to not be in list of available indices.')
        self.assertNotIn(index_code7, aaindex1.record_codes(), f'Expected record {index_code7} to not be in list of available indices.')
        self.assertNotIn(index_code8, aaindex1.record_codes(), f'Expected record {index_code8} to not be in list of available indices.')

    def test_record_names(self):
        """ Test Case to check that the correct names/description are returned for each record. """
        record_names = aaindex1.record_names()
        record_name1 = record_names[50]
        record_name2 = record_names[140]
        record_name3 = record_names[368]
        record_name4 = record_names[560]
#1.)
        self.assertEqual(len(record_names), aaindex1.num_records(), 
            f'Total number of names in the AAi1 database should equal the number of records, got {len(record_names)}.')
#2.)
        self.assertEqual(record_name1, 'alpha-NH chemical shifts (Bundi-Wuthrich, 1979)', 
            f'Incorrect record name returned, expected alpha-NH chemical shifts (Bundi-Wuthrich, 1979). got {record_name1}.')
#3.)
        self.assertTrue(record_name2.startswith('Number of full nonbonding'), 
            f'Incorrect record name returned, expected Number of full nonbonding got {record_name2}.')
#4.)
        self.assertIn('alpha-helix in alpha/beta class', record_name3, 
            f'Expected string alpha-helix in alpha/beta class to be found in record name: {record_name3}.')
#5.)
        self.assertEqual(record_name4, 'Buriability (Zhou-Zhou, 2004)', 
            f'Expected record name to be Buriability (Zhou-Zhou, 2004), got {record_name4}.')

    def test_last_updated(self):
        """ Testing the last updated class attribute which is the latest version of the database, 
            according to https://www.genome.jp/aaindex/. """
        self.assertEqual(aaindex1.last_updated, "February 13, 2017", 
            f"Last updated value does not match expected, got {aaindex1.last_updated}.")

    def test_values(self):
        """ Test Case to check that values() returns the correct amino acid values dict
        for a given record code. """
        index_code1 = 'AURR980103'
        expected_vals = {'A': 1.05, 'L': 0.96, 'R': 0.81, 'K': 0.97, 'N': 0.91,
            'M': 0.99, 'D': 1.39, 'F': 0.95, 'C': 0.6, 'P': 1.05, 'Q': 0.87,
            'S': 0.96, 'E': 1.11, 'T': 1.03, 'G': 1.26, 'W': 1.06, 'H': 1.43,
            'Y': 0.94, 'I': 0.95, 'V': 0.62, '-': 0}
#1.)
        vals = aaindex1.values(index_code1)
        self.assertIsInstance(vals, dict,
            f'values() should return a dict, got {type(vals)}.')
        self.assertEqual(vals, expected_vals,
            f'Values returned do not match expected, got {vals}.')
#2.)
        #verify values() matches direct record attribute access
        self.assertEqual(vals, aaindex1[index_code1]['values'],
            'values() result should match direct record[values] access.')
#3.)
        #invalid record code should raise ValueError
        with self.assertRaises(ValueError):
            aaindex1.values('ABCDEFGH')

    def test_get_record_by_category(self):
        """ Test Case for get_record_by_category(), which returns all records
        belonging to a given category string. """
#1.)
        #sec_struct category should be non-empty and contain known record
        sec_struct_records = aaindex1.get_record_by_category('sec_struct')
        self.assertIsInstance(sec_struct_records, dict,
            'get_record_by_category() should return a dict.')
        self.assertGreater(len(sec_struct_records), 0,
            'Expected non-empty results for sec_struct category.')
        self.assertIn('AURR980103', sec_struct_records,
            'Expected AURR980103 to be in sec_struct category.')
        self.assertIn('FINA770101', sec_struct_records,
            'Expected FINA770101 to be in sec_struct category.')
#2.)
        #flexibility category should be non-empty and contain known record
        flexibility_records = aaindex1.get_record_by_category('flexibility')
        self.assertGreater(len(flexibility_records), 0,
            'Expected non-empty results for flexibility category.')
        self.assertIn('KARP850103', flexibility_records,
            'Expected KARP850103 to be in flexibility category.')
#3.)
        #category lookup is case-insensitive
        self.assertEqual(
            aaindex1.get_record_by_category('SEC_STRUCT'),
            aaindex1.get_record_by_category('sec_struct'),
            'get_record_by_category() should be case-insensitive.')
#4.)
        #non-existent category should raise ValueError
        with self.assertRaises(ValueError):
            aaindex1.get_record_by_category('nonexistent_category')
#5.)
        #non-string input should raise TypeError
        with self.assertRaises(TypeError):
            aaindex1.get_record_by_category(123)

    def test_dunder_methods(self):
        """ Test Case for __len__, __contains__, __iter__, and __repr__ dunder methods. """
#1.) __len__
        self.assertEqual(len(aaindex1), 566,
            f'Expected __len__ to return 566, got {len(aaindex1)}.')
#2.) __contains__
        self.assertIn('AURR980103', aaindex1,
            'Expected AURR980103 in aaindex1 via __contains__.')
        self.assertNotIn('BLAH999999', aaindex1,
            'Expected BLAH999999 to not be in aaindex1 via __contains__.')
#3.) __iter__
        codes_from_iter = list(aaindex1)
        self.assertEqual(len(codes_from_iter), 566,
            f'Expected 566 codes from __iter__, got {len(codes_from_iter)}.')
        self.assertIn('AURR980103', codes_from_iter,
            'Expected AURR980103 in __iter__ output.')
#4.) __repr__
        self.assertEqual(repr(aaindex1), "AAIndex1(records=566, last_updated='February 13, 2017')",
            f'Unexpected __repr__ output: {repr(aaindex1)}.')

    def test_edge_cases(self):
        """ Test edge cases: invalid amino acid key, empty search, whitespace input. """
#1.) accessing a non-existent field should raise AttributeError via Map's __getattr__
        record = aaindex1['AURR980103']
        with self.assertRaises(AttributeError):
            _ = record.nonexistent_field
#2.) values for a valid record should not contain key 'Z' (non-canonical AA)
        vals = aaindex1.values('AURR980103')
        self.assertNotIn('Z', vals,
            'Non-canonical amino acid Z should not be in values dict.')
#3.) empty string search should match all records (every description contains "")
        empty_search = aaindex1.search("")
        self.assertEqual(len(empty_search), aaindex1.num_records(),
            'Empty string search should match all records.')
#4.) non-string input for __getitem__ should raise TypeError
        with self.assertRaises(TypeError):
            aaindex1[12345]
        with self.assertRaises(TypeError):
            aaindex1[None]
#5.) record code with leading/trailing whitespace should still resolve
        record = aaindex1['  AURR980103  ']
        self.assertEqual(record.description,
            "Normalized positional residue frequency at helix termini N' (Aurora-Rose, 1998)",
            'Whitespace-padded record code should resolve correctly.')

    def test_repr_format(self):
        """ Test that __repr__ returns the expected format string. """
        r = repr(aaindex1)
        self.assertTrue(r.startswith('AAIndex1('),
            f'Expected __repr__ to start with AAIndex1(, got {r}.')
        self.assertIn('records=', r,
            'Expected records= in __repr__ output.')
        self.assertIn('last_updated=', r,
            'Expected last_updated= in __repr__ output.')
        self.assertTrue(r.endswith(')'),
            'Expected __repr__ to end with ).')

    def test_parse_aaindex_io_error(self):
        """ Test that _parse_aaindex() raises IOError when the data file cannot be opened. """
        with patch('builtins.open', side_effect=IOError('mocked io error')):
            with self.assertRaises(IOError):
                aaindex1._parse_aaindex()

    def test_parse_aaindex_deprecated(self):
        """ Test that calling the public parse_aaindex() shim emits a DeprecationWarning. """
        import warnings as _warnings
        with _warnings.catch_warnings(record=True) as w:
            _warnings.simplefilter('always')
            aaindex1.parse_aaindex()
        self.assertTrue(
            any(issubclass(warning.category, DeprecationWarning) for warning in w),
            'parse_aaindex() should emit a DeprecationWarning.')

    def test_parse_categories_deprecated(self):
        """ Test that calling the public parse_categories() shim emits a DeprecationWarning. """
        import warnings as _warnings
        with _warnings.catch_warnings(record=True) as w:
            _warnings.simplefilter('always')
            aaindex1.parse_categories()
        self.assertTrue(
            any(issubclass(warning.category, DeprecationWarning) for warning in w),
            'parse_categories() should emit a DeprecationWarning.')

    def test_get_all_categories(self):
        """ Test that get_all_categories() returns the expected structure. """
        cats = aaindex1.get_all_categories()
#1.) should return a non-empty dict
        self.assertIsInstance(cats, dict,
            f'get_all_categories() should return a dict, got {type(cats)}.')
        self.assertGreater(len(cats), 0,
            'get_all_categories() should return a non-empty dict.')
#2.) known record codes should map to known category strings
        known_categories = {'sec_struct', 'hydrophobic', 'flexibility', 'charge',
                            'alpha_helix', 'beta_sheet', 'residue_prop', 'other'}
        self.assertIn('AURR980103', cats,
            'Expected AURR980103 in categories dict.')
        self.assertIn(cats['AURR980103'], known_categories,
            f"Category for AURR980103 should be a known category, got {cats['AURR980103']}.")
#3.) number of category entries should equal number of records
        self.assertEqual(len(cats), aaindex1.num_records(),
            f'Number of category entries should equal num_records(), got {len(cats)}.')

    def test_encode(self):
        """ Test Case for encode(), which maps a protein sequence to a numeric vector
        using a given AAindex1 record. """
        index_code = 'KYTJ820101'  # Kyte-Doolittle hydrophobicity
#1.) basic encoding — known values for A, C, D, E
        result = aaindex1.encode('ACDE', index_code)
        self.assertIsInstance(result, list,
            f'encode() should return a list, got {type(result)}.')
        self.assertEqual(len(result), 4,
            f'Encoded vector should have length 4, got {len(result)}.')
        self.assertAlmostEqual(result[0], 1.8,
            msg=f'Expected A=1.8 for {index_code}, got {result[0]}.')
        self.assertAlmostEqual(result[1], 2.5,
            msg=f'Expected C=2.5 for {index_code}, got {result[1]}.')
        self.assertAlmostEqual(result[2], -3.5,
            msg=f'Expected D=-3.5 for {index_code}, got {result[2]}.')
        self.assertAlmostEqual(result[3], -3.5,
            msg=f'Expected E=-3.5 for {index_code}, got {result[3]}.')
#2.) single amino acid
        single = aaindex1.encode('G', index_code)
        self.assertEqual(len(single), 1,
            'encode() of single AA should return a length-1 list.')
#3.) empty sequence returns empty list
        empty = aaindex1.encode('', index_code)
        self.assertEqual(empty, [],
            'encode() of empty string should return [].')
#4.) lowercase sequence is normalised
        lower = aaindex1.encode('acde', index_code)
        self.assertEqual(lower, result,
            'encode() should normalise lowercase to uppercase.')
#5.) non-canonical amino acid maps to None
        unknown = aaindex1.encode('Z', index_code)
        self.assertIsNone(unknown[0],
            'Non-canonical amino acid Z should map to None.')
#6.) gap character uses gap_value (default 0.0)
        gapped = aaindex1.encode('-', index_code)
        self.assertEqual(gapped[0], 0.0,
            'Gap character should map to the default gap_value 0.0.')
        gapped_custom = aaindex1.encode('-', index_code, gap_value=-1.0)
        self.assertEqual(gapped_custom[0], -1.0,
            'Gap character should map to the custom gap_value -1.0.')
#7.) NA values in source data encode as None
        na_record = 'AVBF000101'  # has NA for proline
        na_result = aaindex1.encode('P', na_record)
        self.assertIsNone(na_result[0],
            'encode() should return None for an NA-valued amino acid.')
#8.) non-string sequence raises TypeError
        with self.assertRaises(TypeError):
            aaindex1.encode(123, index_code)
#9.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex1.encode('ACG', 'BLAH999999')

    def test_get_correlated_indices(self):
        """ Test Case for get_correlated_indices(), which traverses correlation chains. """
        seed = 'KYTJ820101'  # has several correlated neighbours
#1.) depth=1 returns at least one correlated record
        result = aaindex1.get_correlated_indices(seed, depth=1)
        self.assertIsInstance(result, dict,
            f'get_correlated_indices() should return a dict, got {type(result)}.')
        self.assertGreater(len(result), 0,
            'Expected at least one correlated record for KYTJ820101 at depth=1.')
#2.) seed itself is not in the output
        self.assertNotIn(seed, result,
            'Seed record should not appear in get_correlated_indices() output.')
#3.) all returned records exist in the database
        for code in result:
            self.assertIn(code, aaindex1,
                f'Correlated record {code} not found in AAindex1.')
#4.) correlation values are floats
        for code, coeff in result.items():
            self.assertIsInstance(coeff, float,
                f'Correlation coefficient for {code} should be float, got {type(coeff)}.')
#5.) depth=2 returns >= depth-1 results (i.e. at least as many as depth=1)
        result_d2 = aaindex1.get_correlated_indices(seed, depth=2)
        self.assertGreaterEqual(len(result_d2), len(result),
            'depth=2 should return at least as many results as depth=1.')
#6.) min_correlation=1.0 should return very few (or zero) results
        strict = aaindex1.get_correlated_indices(seed, min_correlation=1.0)
        self.assertIsInstance(strict, dict,
            'get_correlated_indices() with min_correlation=1.0 should return a dict.')
#7.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex1.get_correlated_indices('BLAH999999')
#8.) invalid min_correlation type raises TypeError
        with self.assertRaises(TypeError):
            aaindex1.get_correlated_indices(seed, min_correlation='high')
#9.) out-of-range min_correlation raises ValueError
        with self.assertRaises(ValueError):
            aaindex1.get_correlated_indices(seed, min_correlation=1.5)
#10.) depth < 1 raises ValueError
        with self.assertRaises(ValueError):
            aaindex1.get_correlated_indices(seed, depth=0)

    def test_compare_indices(self):
        """ Test Case for compare_indices(), which computes Pearson r between two records. """
        code1 = 'KYTJ820101'  # Kyte-Doolittle hydrophobicity
        code2 = 'EISD840101'  # Eisenberg normalised consensus hydrophobicity
#1.) result is a float in [-1, 1]
        r = aaindex1.compare_indices(code1, code2)
        self.assertIsNotNone(r,
            f'compare_indices({code1}, {code2}) should not return None.')
        self.assertIsInstance(r, float,
            f'compare_indices() should return a float, got {type(r)}.')
        self.assertGreaterEqual(r, -1.0,
            f'Pearson r must be >= -1, got {r}.')
        self.assertLessEqual(r, 1.0,
            f'Pearson r must be <= 1, got {r}.')
#2.) a record correlated with itself is 1.0
        self_r = aaindex1.compare_indices(code1, code1)
        self.assertAlmostEqual(self_r, 1.0, places=9,
            msg='compare_indices() of a record with itself should return 1.0.')
#3.) KYTJ820101 and EISD840101 are known high-correlation neighbours
        self.assertGreater(r, 0.8,
            f'Expected r > 0.8 between {code1} and {code2}, got {r}.')
#4.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex1.compare_indices('BLAH999999', code2)
        with self.assertRaises(ValueError):
            aaindex1.compare_indices(code1, 'BLAH999999')

    def test_to_dict(self):
        """ Test Case for to_dict(), exporting records as plain Python dicts. """
        code = 'AURR980103'
#1.) single record export
        d = aaindex1.to_dict(code)
        self.assertIsInstance(d, dict,
            f'to_dict() should return a dict, got {type(d)}.')
        self.assertIn(code, d,
            f'Expected {code} as key in to_dict() result.')
        self.assertIn('description', d[code],
            'Expected description field in exported record.')
        self.assertIn('values', d[code],
            'Expected values field in exported record.')
#2.) full database export
        full = aaindex1.to_dict()
        self.assertEqual(len(full), aaindex1.num_records(),
            f'Full to_dict() should contain all {aaindex1.num_records()} records.')
        self.assertIn(code, full,
            f'Expected {code} in full to_dict() export.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex1.to_dict('BLAH999999')

    def test_to_json(self):
        """ Test Case for to_json(), serialising records to a JSON string. """
        import json
        code = 'AURR980103'
#1.) single record serialisation returns valid JSON
        j = aaindex1.to_json(code)
        self.assertIsInstance(j, str,
            f'to_json() should return a str, got {type(j)}.')
        parsed = json.loads(j)
        self.assertIn(code, parsed,
            f'Expected {code} as key in parsed to_json() output.')
#2.) full database serialisation is valid JSON
        full_json = aaindex1.to_json()
        full_parsed = json.loads(full_json)
        self.assertEqual(len(full_parsed), aaindex1.num_records(),
            f'Full to_json() should contain all {aaindex1.num_records()} records.')
#3.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex1.to_json('BLAH999999')

    def test_to_dataframe(self):
        """ Test Case for to_dataframe(), exporting values as a pandas DataFrame.
        Skipped when pandas is not installed. """
        try:
            import pandas as pd
        except ImportError:
            self.skipTest('pandas not installed — skipping test_to_dataframe.')
        code = 'AURR980103'
#1.) single record export
        df = aaindex1.to_dataframe(code)
        self.assertIsInstance(df, pd.DataFrame,
            f'to_dataframe() should return a DataFrame, got {type(df)}.')
        self.assertEqual(len(df), 1,
            'Single-record to_dataframe() should have 1 row.')
        self.assertIn(code, df.index,
            f'Expected {code} in DataFrame index.')
#2.) full database export has one row per record
        full_df = aaindex1.to_dataframe()
        self.assertEqual(len(full_df), aaindex1.num_records(),
            f'Full to_dataframe() should have {aaindex1.num_records()} rows.')
#3.) amino acid columns are present
        for aa in ['A', 'C', 'D', 'E']:
            self.assertIn(aa, full_df.columns,
                f'Expected column {aa} in to_dataframe() output.')
#4.) invalid record code raises ValueError
        with self.assertRaises(ValueError):
            aaindex1.to_dataframe('BLAH999999')

    def test_valid_categories(self):
        """ Test that valid_categories() returns a sorted list of all distinct categories. """
        cats = aaindex1.valid_categories()
#1.) return type
        self.assertIsInstance(cats, list,
            f'valid_categories() should return a list, got {type(cats)}.')
#2.) non-empty
        self.assertGreater(len(cats), 0,
            'valid_categories() should not be empty.')
#3.) sorted
        self.assertEqual(cats, sorted(cats),
            'valid_categories() should return a sorted list.')
#4.) contains known categories
        for expected in ('sec_struct', 'hydrophobic'):
            self.assertIn(expected, cats,
                f'Expected category "{expected}" in valid_categories().')
#5.) no duplicates
        self.assertEqual(len(cats), len(set(cats)),
            'valid_categories() should not contain duplicate entries.')

    def test_search_fuzzy(self):
        """ Test Case for search_fuzzy() fuzzy search with ranked results. """
#1.) returns a dict
        results = aaindex1.search_fuzzy('hydrophobicity')
        self.assertIsInstance(results, dict,
            f'search_fuzzy() should return a dict, got {type(results)}.')
#2.) non-empty results for a known term
        self.assertGreater(len(results), 0,
            'search_fuzzy() should return at least one result for "hydrophobicity".')
#3.) respects n limit
        results_n = aaindex1.search_fuzzy('hydrophobicity', n=3)
        self.assertLessEqual(len(results_n), 3,
            'search_fuzzy() should return at most n results.')
#4.) cutoff filters results
        results_high = aaindex1.search_fuzzy('hydrophobicity', cutoff=0.99)
        results_low = aaindex1.search_fuzzy('hydrophobicity', cutoff=0.0)
        self.assertLessEqual(len(results_high), len(results_low),
            'Higher cutoff should return fewer or equal results.')
#5.) TypeError for non-str query
        with self.assertRaises(TypeError):
            aaindex1.search_fuzzy(123)
#6.) ValueError for n < 1
        with self.assertRaises(ValueError):
            aaindex1.search_fuzzy('hydrophobicity', n=0)
#7.) ValueError for cutoff outside [0, 1]
        with self.assertRaises(ValueError):
            aaindex1.search_fuzzy('hydrophobicity', cutoff=1.5)

    def test_map_repr_str(self):
        """ Test that Map.__repr__ and __str__ produce a structured AAIndexRecord() string. """
        record = aaindex1['AURR980103']
        repr_str = repr(record)
        str_str = str(record)
#1.) repr starts with AAIndexRecord(
        self.assertTrue(repr_str.startswith('AAIndexRecord('),
            f'repr should start with "AAIndexRecord(", got: {repr_str[:30]}.')
#2.) contains description and pmid keys
        self.assertIn('description:', repr_str,
            'repr should contain "description:".')
        self.assertIn('pmid:', repr_str,
            'repr should contain "pmid:".')
#3.) values dict is summarised not fully expanded
        self.assertIn('values:', repr_str,
            'repr should mention "values:".')
        self.assertIn('entries)', repr_str,
            'repr should show entry count for the values dict.')
#4.) str and repr produce the same output
        self.assertEqual(repr_str, str_str,
            '__repr__ and __str__ should produce identical output.')

    def test_richer_error_messages(self):
        """ Test that __getitem__ includes a did-you-mean hint on invalid accession. """
#1.) completely unknown code raises ValueError
        with self.assertRaises(ValueError):
            aaindex1['BLAH999999']
#2.) near-miss code produces a hint mentioning possible matches
        try:
            aaindex1['AURR980100']  # close to AURR980103
        except ValueError as e:
            err_msg = str(e)
            self.assertIn('Did you mean', err_msg,
                f'Expected did-you-mean hint in error message, got: {err_msg}.')

    def test_batch_encode(self):
        """ Test Case for batch_encode(), which encodes multiple sequences at once. """
        index_code = 'KYTJ820101'
#1.) basic batch — same result as calling encode() per sequence
        seqs = ['ACDE', 'FG']
        result = aaindex1.batch_encode(seqs, index_code)
        self.assertIsInstance(result, list,
            f'batch_encode() should return a list, got {type(result)}.')
        self.assertEqual(len(result), 2,
            'batch_encode() should return one vector per sequence.')
        self.assertEqual(result[0], aaindex1.encode('ACDE', index_code),
            'batch_encode() vector should match encode() for the same sequence.')
        self.assertEqual(result[1], aaindex1.encode('FG', index_code),
            'batch_encode() vector should match encode() for the same sequence.')
#2.) empty sequence list returns empty list
        self.assertEqual(aaindex1.batch_encode([], index_code), [],
            'batch_encode() of empty list should return [].')
#3.) single-element list
        single = aaindex1.batch_encode(['ACDE'], index_code)
        self.assertEqual(len(single), 1,
            'batch_encode() of one sequence should return a length-1 list.')
#4.) record is looked up once: all vectors use the same values
        result2 = aaindex1.batch_encode(['A', 'A', 'A'], index_code)
        self.assertTrue(all(v == result2[0] for v in result2),
            'All single-AA vectors for the same amino acid should be equal.')
#5.) gap character uses gap_value
        gapped = aaindex1.batch_encode(['-'], index_code, gap_value=-99.0)
        self.assertEqual(gapped[0][0], -99.0,
            'Gap character should use the custom gap_value.')
#6.) non-canonical amino acid maps to None
        unknown = aaindex1.batch_encode(['Z'], index_code)
        self.assertIsNone(unknown[0][0],
            'Non-canonical amino acid Z should map to None in batch_encode().')
#7.) lowercase sequences are normalised
        lower = aaindex1.batch_encode(['acde'], index_code)
        self.assertEqual(lower, aaindex1.batch_encode(['ACDE'], index_code),
            'batch_encode() should normalise lowercase to uppercase.')
#8.) TypeError if sequences is not a list
        with self.assertRaises(TypeError):
            aaindex1.batch_encode('ACDE', index_code)
#9.) TypeError if an element is not a str
        with self.assertRaises(TypeError):
            aaindex1.batch_encode(['ACDE', 123], index_code)
#10.) ValueError for unknown record code
        with self.assertRaises(ValueError):
            aaindex1.batch_encode(['ACDE'], 'BLAH999999')

    def test_instance_version(self):
        """ Test that AAIndex1, AAIndex2, and AAIndex3 instances expose a version attribute. """
        from aaindex import aaindex2, aaindex3
#1.) version attribute is a string
        self.assertIsInstance(aaindex1.version, str,
            f'aaindex1.version should be a str, got {type(aaindex1.version)}.')
        self.assertIsInstance(aaindex2.version, str,
            f'aaindex2.version should be a str, got {type(aaindex2.version)}.')
        self.assertIsInstance(aaindex3.version, str,
            f'aaindex3.version should be a str, got {type(aaindex3.version)}.')
#2.) version matches the package __version__
        self.assertEqual(aaindex1.version, __version__,
            f'aaindex1.version should match __version__, got {aaindex1.version}.')
        self.assertEqual(aaindex2.version, __version__,
            f'aaindex2.version should match __version__, got {aaindex2.version}.')
        self.assertEqual(aaindex3.version, __version__,
            f'aaindex3.version should match __version__, got {aaindex3.version}.')

    def test_parse_categories_path_validation(self):
        """ Test that _parse_categories() rejects relative paths with directory components. """
#1.) relative traversal path should raise ValueError
        with self.assertRaises(ValueError):
            aaindex1._parse_categories('../../etc/passwd')
#2.) path with subdirectory component should raise ValueError
        with self.assertRaises(ValueError):
            aaindex1._parse_categories('subdir/file.txt')
#3.) bare filename is accepted (resolved against data dir)
        result = aaindex1._parse_categories('aaindex_to_category.txt')
        self.assertIsInstance(result, dict,
            '_parse_categories() with bare filename should return a dict.')

    def test_property_setter_validation(self):
        """ Test that data_dir and aaindex_filename setters reject path traversal values. """
        from aaindex.aaindex1 import AAIndex1 as _AAIndex1
        fresh = _AAIndex1.__new__(_AAIndex1)
#1.) relative path in data_dir should raise ValueError
        with self.assertRaises(ValueError):
            fresh.data_dir = '../../etc'
#2.) path with separator in aaindex_filename should raise ValueError
        with self.assertRaises(ValueError):
            fresh.aaindex_filename = '../aaindex1'
#3.) simple names are accepted
        fresh.data_dir = 'data'
        self.assertEqual(fresh.data_dir, 'data')
        fresh.aaindex_filename = 'aaindex1'
        self.assertEqual(fresh.aaindex_filename, 'aaindex1')

    def test_correlation_coefficients_are_floats(self):
        """ Test that correlation_coefficients values are stored as floats, not strings. """
        record = aaindex1['FINA770101']
        cc = record['correlation_coefficients']
        self.assertGreater(len(cc), 0,
            'FINA770101 should have at least one correlation coefficient.')
        for code, val in cc.items():
            self.assertIsInstance(val, float,
                f'correlation_coefficients[{code!r}] should be float, got {type(val)}.')

    def test_gap_placeholder_is_float(self):
        """ Test that the gap placeholder '-' is stored as float 0.0, not int 0. """
        record = aaindex1['AURR980103']
        gap_val = record['values']['-']
        self.assertIsInstance(gap_val, float,
            f"Gap placeholder should be float, got {type(gap_val)}.")
        self.assertEqual(gap_val, 0.0,
            f"Gap placeholder should equal 0.0, got {gap_val}.")


if __name__ == '__main__':
    #run all unit tests
    unittest.main(verbosity=2)