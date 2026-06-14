# aaindex - Python package for working with the AAindex database (https://www.genome.jp/aaindex/) <a name="TOP"></a>
<p align="center">
  <img src="https://raw.githubusercontent.com/amckenna41/aaindex/main/images/aaindex_logo.png" height="400" />
</p>

[![AAindex](https://img.shields.io/pypi/v/aaindex)](https://pypi.org/project/aaindex/)
[![pytest](https://github.com/amckenna41/aaindex/workflows/Building%20and%20Testing/badge.svg)](https://github.com/amckenna41/aaindex/actions?query=workflowBuilding%20and%20Testing)
[![PythonV](https://img.shields.io/pypi/pyversions/aaindex?logo=2)](https://pypi.org/project/aaindex/)
[![Platforms](https://img.shields.io/badge/platforms-linux%2C%20macOS%2C%20Windows-green)](https://pypi.org/project/aaindex/)
[![Documentation Status](https://readthedocs.org/projects/aaindex/badge/?version=latest)](https://aaindex.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)
[![Issues](https://img.shields.io/github/issues/amckenna41/aaindex)](https://github.com/amckenna41/aaindex/issues)
<!-- [![Size](https://img.shields.io/github/repo-size/amckenna41/aaindex)](https://github.com/amckenna41/aaindex) -->
<!-- [![codecov](https://codecov.io/gh/amckenna41/aaindex/branch/main/graph/badge.svg?token=SM2ZKPN8PZ)](https://codecov.io/gh/amckenna41/aaindex) -->

Table of Contents
-----------------

  * [Introduction](#introduction)
  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Documentation](#documentation)
  * [Tests](#tests)
  * [Directories](#directories)
  * [Contact](#contact)
  * [License](#license)
  * [References](#References)

Introduction
------------
The AAindex is a database of numerical indices representing various physicochemical, structural and biochemical properties of amino acids and pairs of amino acids 🧬. The AAindex consists of three sections: AAindex1 for the amino acid index of 20 numerical values, AAindex2 for the amino acid mutation matrix and AAindex3 for the statistical protein contact potentials. All data are derived from published literature [[1]](#references). 

This `aaindex` Python software package is a very lightweight way of accessing the data represented in the various AAindex databases, requiring no additional external library installations. Any record within the 3 databases and their associated data/numerical indices can be accessed in one simple command. The package supports all three AAindex databases: AAindex1 (amino acid property indices), AAindex2 (substitution matrices), and AAindex3 (contact potential matrices).

* 💻 A quick Colab notebook demo of `aaindex` is available [here][demo]. 
* 📝 A **Medium** article that dives deeper into the AAindex and the `aaindex` software itself is available [here][medium].

<!-- <strong>A demo of the software is available [here](https://colab.research.google.com/drive/1dccV_n1BRMiU8W13F9PPXbSaFzvOdQLC?usp=sharing). </strong> -->
<!-- <strong>A medium article outlining the background and usage of the aaindex database and software is available [here]().</strong> -->

Background
----------

**AAindex1:**

The AAindex1 section currently contains 566 amino acid indices representing the various physicochemical, structural and biochemical properties of amino acids. Each entry consists of an accession number, a short description on the index, the reference information, notes, PMID (pubmed ID) and the numerical values for the property of 20 amino acids. In addition, it contains neighbour information; namely, the cross-links to other entries with an absolute value for the correlation coefficient of 0.8 or larger. With the links the user can identify a set of entries describing similar properties An example of the format of an AAindex1 record can be seen within the [aaindex](https://github.com/amckenna41/aaindex/tree/main/aaindex) folder [[1]](#references).

  ```python
  ************************************************************************
  *                                                                      *
  * H Accession number                                                   *
  * D Data description                                                   *
  * R Pub med article ID (PMID)                                          *
  * A Author(s)                                                          *
  * T Title of the article                                               *
  * J Journal reference                                                  *
  * * Comment or missing                                                 *
  * C Accession numbers of similar entries with the correlation          *
  *   coefficients of 0.8 (-0.8) or more (less).                         *
  *   Notice: The correlation coefficient is calculated with zeros       *
  *   filled for missing values.                                         *
  * I Amino acid index data in the following order                       *
  *   Ala    Arg    Asn    Asp    Cys    Gln    Glu    Gly    His    Ile *
  *   Leu    Lys    Met    Phe    Pro    Ser    Thr    Trp    Tyr    Val *
  * //                                                                   *
  ************************************************************************
  ```

**AAindex2:**

The AAindex2 section currently contains 66 amino acid mutation matrices: 47 symmetric matrices and 19 non-symmetric matrices. The format of the entry is almost the same as that of AAindex1 except that it contains 210 numerical values (20 diagonal and 20 × 19/2 off-diagonal elements) for a symmetric matrix and 400 or more numerical values for a non-symmetric matrix (some matrices include a gap or distinguish two states of cysteine). An example of the format of an AAindex2 record can be seen within the [aaindex](https://github.com/amckenna41/aaindex/tree/main/aaindex) folder.

**AAindex3:**

The AAindex3 section contains 47 statistical protein contact potentials and follows the same record format to that of the AAindex2. An example of the format of an AAindex3 record can be seen within the [aaindex](https://github.com/amckenna41/aaindex/tree/main/aaindex) folder.

<!-- ![alt text](https://raw.githubusercontent.com/amckenna41/aaindex/main/images/aaindex_example.png) -->

  ```python  
  ************************************************************************
  *                                                                      *
  * Each entry has the following format.                                 *
  *                                                                      *
  * H Accession number                                                   *
  * D Data description                                                   *
  * R PMID                                                               *
  * A Author(s)                                                          *
  * T Title of the article                                               *
  * J Journal reference                                                  *
  * * Comment or missing                                                 *
  * M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV           *
  *   AA                                                                 *
  *   AR RR                                                              *
  *   AN RN NN                                                           *
  *   AD RD ND DD                                                        *
  *   AC RC NC DC CC                                                     *
  *   AQ RQ NQ DQ CQ QQ                                                  *
  *   AE RE NE DE CE QE EE                                               *
  *   AG RG NG DG CG QG EG GG                                            *
  *   AH RH NH DH CH QH EH GH HH                                         *
  *   AI RI NI DI CI QI EI GI HI II                                      *
  *   AL RL NL DL CL QL EL GL HL IL LL                                   *
  *   AK RK NK DK CK QK EK GK HK IK LK KK                                *
  *   AM RM NM DM CM QM EM GM HM IM LM KM MM                             *
  *   AF RF NF DF CF QF EF GF HF IF LF KF MF FF                          *
  *   AP RP NP DP CP QP EP GP HP IP LP KP MP FP PP                       *
  *   AS RS NS DS CS QS ES GS HS IS LS KS MS FS PS SS                    *
  *   AT RT NT DT CT QT ET GT HT IT LT KT MT FT PT ST TT                 *
  *   AW RW NW DW CW QW EW GW HW IW LW KW MW FW PW SW TW WW              *
  *   AY RY NY DY CY QY EY GY HY IY LY KY MY FY PY SY TY WY YY           *
  *   AV RV NV DV CV QV EV GV HV IV LV KV MV FV PV SV TV WV YV VV        *
  * //                                                                   *
  ************************************************************************
```

Installation
-----------------
Install the latest version of `aaindex` using pip:

```bash
pip3 install aaindex --upgrade
```

Install by cloning the repository:

```bash
git clone https://github.com/amckenna41/aaindex.git
cd aaindex
pip install .
```

Usage
-----
The `aaindex` package is made up of three modules for each AAindex database, with each having a Python class of the same name, when importing the package you should import the required database module:

```python
from aaindex import aaindex1
# from aaindex import aaindex2
# from aaindex import aaindex3
```

## AAIndex1 Usage

### Get record from AAindex1
The AAindex1 class offers diverse functionalities for obtaining any element from any record in the database. The records are imported from a parsed json in the data folder of the package. You can search for a particular record by its record code/accession number or its name/description. You can also get the record category, references, notes, correlation coefficients, PMID and importantly its associated amino acid values:
```python
from aaindex import aaindex1

full_record = aaindex1['CHOP780206']   #get full AAI record
''' full_record ->
{'category': 'sec_struct', 
'correlation_coefficients': {}, 
'description': 'Normalized frequency of N-terminal non helical region (Chou-Fasman, 1978b)', 
'notes': '', 
'pmid': '364941', 
'references': "Chou, P.Y. and Fasman, G.D. 'Prediction of the secondary structure of proteins from their amino acid sequence' Adv. Enzymol. 47, 45-148 (1978)", 
'values': {'-': 0, 'A': 0.7, 'C': 0.65, 'D': 0.98, 'E': 1.04, 'F': 0.93, 'G': 1.41, 'H': 1.22, 'I': 0.78, 'K': 1.01, 'L': 0.85, 'M': 0.83, 'N': 1.42, 'P': 1.1, 'Q': 0.75, 'R': 0.34, 'S': 1.55, 'T': 1.09, 'V': 0.75, 'W': 0.62, 'Y': 0.99}}
'''

#get individual elements of AAindex record
record_values = aaindex1['CHOP780206']['values'] 
record_values = aaindex1['CHOP780206'].values
#'values': {'-': 0, 'A': 0.7, 'C': 0.65, 'D': 0.98, 'E': 1.04, 'F': 0.93, 'G': 1.41, 'H': 1.22, 'I': 0.78, 'K': 1.01, 'L': 0.85, 'M': 0.83, 'N': 1.42, 'P': 1.1, 'Q': 0.75, 'R': 0.34, 'S': 1.55, 'T': 1.09, 'V': 0.75, 'W': 0.62, 'Y': 0.99}

record_description = aaindex1['CHOP780206']['description']
record_description = aaindex1['CHOP780206'].description
#'description': 'Normalized frequency of N-terminal non helical region (Chou-Fasman, 1978b)'

record_references = aaindex1['CHOP780206']['references']
record_references = aaindex1['CHOP780206'].references
#'references': "Chou, P.Y. and Fasman, G.D. 'Prediction of the secondary structure of proteins from their amino acid sequence' Adv. Enzymol. 47, 45-148 (1978)"

record_notes = aaindex1['CHOP780206']['notes']
record_notes = aaindex1['CHOP780206'].notes
#""

record_correlation_coefficients = aaindex1['CHOP780206']['correlation_coefficients']
record_correlation_coefficients = aaindex1['CHOP780206'].correlation_coefficients
#{}

record_pmid = aaindex1['CHOP780206']['pmid']  
record_pmid = aaindex1['CHOP780206'].pmid
#364941

record_category = aaindex1['CHOP780206']['category']
record_category = aaindex1['CHOP780206'].category
#sec_struct
```

### Get total number of AAindex1 records
```python
aaindex1.num_records()
```

### Get list of all AAindex1 record codes
```python
aaindex1.record_codes()
```

### Get list of all AAindex1 record names
```python
aaindex1.record_names()
```

### Get amino acid values for a record
```python
# Shortcut to retrieve only the values dict without fetching the full record
aaindex1.values('CHOP780206')
# {'-': 0, 'A': 0.7, 'C': 0.65, 'D': 0.98, ...}
```

### Search records by keyword
```python
# Search with a single keyword (case-insensitive)
aaindex1.search('hydrophobicity')   # dict of matching records

# Search with multiple keywords — returns records matching any of the terms
aaindex1.search(['hydrophobicity', 'charge'])   # dict of matching records
```

### Get records by category
> **Note:** `get_record_by_category()` is available on `AAIndex1` only. `AAIndex2` and `AAIndex3` records do not include a category field.
```python
# Retrieve all records belonging to a given category (case-insensitive)
aaindex1.get_record_by_category('sec_struct')   # dict of matching records
aaindex1.get_record_by_category('hydrophobicity')
```

### Get list of amino acid single-letter codes
```python
aaindex1.amino_acids()
# ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
```

### Encode a protein sequence as a numeric feature vector
```python
# Map each amino acid in a sequence to its index value
aaindex1.encode('ACDEFGHIKLMNPQRSTVWY', 'KYTJ820101')
# [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, -1.3, -3.9, 3.8, ...]

# NA values encode as None; gaps use gap_value (default 0.0)
aaindex1.encode('ACGP-', 'AVBF000101', gap_value=0.0)
```

### Find correlated indices
```python
# Get all indices directly correlated with a seed (|r| >= 0.8 by default)
aaindex1.get_correlated_indices('KYTJ820101')
# {'EISD840101': 0.949, 'JOND750101': 0.838, ...}

# Expand to 2-hop neighbours
aaindex1.get_correlated_indices('KYTJ820101', depth=2)

# Raise the threshold
aaindex1.get_correlated_indices('KYTJ820101', min_correlation=0.9)
```

### Compare two indices (Pearson r)
```python
# Compute Pearson correlation between any two AAindex1 records
aaindex1.compare_indices('KYTJ820101', 'EISD840101')   # e.g. 0.949
aaindex1.compare_indices('KYTJ820101', 'KYTJ820101')   # 1.0
```

### Export to dict, JSON, or pandas DataFrame
```python
# Export one record or the full database as a Python dict
aaindex1.to_dict('KYTJ820101')   # single record
aaindex1.to_dict()               # all 566 records

# Export as a JSON string
aaindex1.to_json('KYTJ820101')
aaindex1.to_json()

# Export amino acid values as a pandas DataFrame (pandas must be installed)
aaindex1.to_dataframe('KYTJ820101')  # DataFrame with 1 row, amino acid columns
aaindex1.to_dataframe()              # DataFrame with 566 rows

# AAindex2 / AAindex3 export the pairwise matrix
aaindex2.to_dict('ALTS910101')
aaindex2.to_dataframe('ALTS910101')  # 20×20 DataFrame
aaindex2.to_dataframe()              # MultiIndex DataFrame (record, row_aa)
```

### Built-in protocol support
```python
# Check membership
'CHOP780206' in aaindex1   # True

# Get total number of records
len(aaindex1)              # 566

# Iterate over all accession numbers
for record_code in aaindex1:
    print(record_code)
```

## AAIndex2 Usage
```python
from aaindex import aaindex2

# Get number of records, record codes, and record names
aaindex2.num_records()            # 94
aaindex2.record_codes()           # sorted list of all accession numbers
aaindex2.record_names()           # list of all record descriptions

# Get a full record by accession number
record = aaindex2['ALTS910101']
record.description                # 'The PAM-120 matrix (Altschul, 1991)'
record.matrix                     # nested dict of 20x20 substitution scores

# Look up a pairwise substitution score (symmetric)
aaindex2.get('ALTS910101', 'A', 'R')   # -3.0
aaindex2.get('ALTS910101', 'R', 'A')   # -3.0

# Get just the matrix dict for a record
aaindex2.values('ALTS910101')

# Get list of amino acid single-letter codes
aaindex2.amino_acids()   # ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Search records by keyword — accepts a single string or a list of keywords
aaindex2.search('substitution')              # dict of matching records
aaindex2.search(['substitution', 'PAM'])     # records matching any keyword

# Built-in protocol support
'ALTS910101' in aaindex2   # True
len(aaindex2)              # 94
for record_code in aaindex2:
    print(record_code)
```

## AAIndex3 Usage
```python
from aaindex import aaindex3

# Get number of records, record codes, and record names
aaindex3.num_records()            # 47
aaindex3.record_codes()           # sorted list of all accession numbers
aaindex3.record_names()           # list of all record descriptions

# Get a full record by accession number
record = aaindex3['TANS760101']
record.description                # 'Statistical contact potential derived from 25 x-ray protein structures'
record.matrix                     # nested dict of 20x20 contact potential scores

# Look up a pairwise contact potential (symmetric)
aaindex3.get('TANS760101', 'A', 'A')   # -2.6
aaindex3.get('TANS760101', 'A', 'R')   # -3.4

# Get just the matrix dict for a record
aaindex3.values('TANS760101')

# Get list of amino acid single-letter codes
aaindex3.amino_acids()   # ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Search records by keyword — accepts a single string or a list of keywords
aaindex3.search('contact potential')                    # dict of matching records
aaindex3.search(['contact potential', 'statistical'])   # records matching any keyword

# Built-in protocol support
'TANS760101' in aaindex3   # True
len(aaindex3)              # 47
for record_code in aaindex3:
    print(record_code)
```

Documentation 📖
----------------
Full API documentation is available on [Read the Docs][readthedocs].


Tests 🧪
--------
To run all tests, from the main `aaindex` folder run:
```
python3 -m unittest discover tests
```

Directories 📁
--------------
* `/tests` - unit and integration tests for `aaindex` package.
* `/aaindex` - source code and all required external data files for package.
* `/images` - images used throughout README.
* `/docs` - `aaindex` documentation.
 

Contact ✉️
---------
If you have any questions or comments, please contact amckenna41@qub.ac.uk or raise an issue on the [Issues][Issues] tab.

License
-------
Distributed under the MIT License. See `LICENSE` for more details.  

References
----------
\[1\]: Shuichi Kawashima, Minoru Kanehisa, AAindex: Amino Acid index database, Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Page 374, https://doi.org/10.1093/nar/28.1.374 <br>
\[2\]: https://www.genome.jp/aaindex/ <br>
\[3\]: Nakai, K., Kidera, A., and Kanehisa, M.; Cluster analysis of amino acid indices for prediction of protein structure and function. Protein Eng. 2, 93-100 (1988). [PMID:3244698] <br>
\[4\]: Tomii, K. and Kanehisa, M.; Analysis of amino acid indices and mutation matrices for sequence comparison and structure prediction of proteins. Protein Eng. 9, 27-36 (1996). [PMID:9053899] <br>
\[5\]: Kawashima, S., Ogata, H., and Kanehisa, M.; AAindex: amino acid index database. Nucleic Acids Res. 27, 368-369 (1999). [PMID:9847231] <br>
\[6\]: Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. Nucleic Acids Res. 28, 374 (2000). [PMID:10592278] <br>
\[7\]: Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., and Kanehisa, M.; AAindex: amino acid index database, progress report 2008. Nucleic Acids Res. 36, D202-D205 (2008). [PMID:17998252] 

[<img src="https://img.shields.io/github/stars/amckenna41/aaindex?color=green&label=star%20it%20on%20GitHub" width="132" height="20" alt="Star it on GitHub">](https://github.com/amckenna41/aaindex)


<a href="https://www.buymeacoffee.com/amckenna41" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>


[Back to top](#TOP)

[python]: https://www.python.org/downloads/release/python-360/
[aaindex]: https://github.com/amckenna41/aaindex
[requests]: https://requests.readthedocs.io/en/latest/
[numpy]: https://numpy.org/
[PyPi]: https://pypi.org/project/aaindex/
[demo]: https://colab.research.google.com/drive/1dccV_n1BRMiU8W13F9PPXbSaFzvOdQLC?usp=sharing
[medium]: https://medium.com/@ajmckenna69/aaindex-a63de37ec118
[Issues]: https://github.com/amckenna41/aaindex/issues
[readthedocs]: https://aaindex.readthedocs.io/en/latest/
[changelog]: https://github.com/amckenna41/aaindex/blob/main/CHANGELOG.md