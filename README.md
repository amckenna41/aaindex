## Python package for working with the AAindex database (https://www.genome.jp/aaindex/) <a name="TOP"></a>
<p align="center">
  <img src="https://raw.githubusercontent.com/amckenna41/aaindex/main/images/aaindex_logo.png" />
</p>

[![AAindex](https://img.shields.io/pypi/v/aaindex)](https://pypi.org/project/aaindex/)
[![pytest](https://github.com/amckenna41/aaindex/workflows/Building%20and%20Testing/badge.svg)](https://github.com/amckenna41/aaindex/actions?query=workflowBuilding%20and%20Testing)
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/amckenna41/aaindex/tree/main.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/amckenna41/aaindex/tree/main)
[![PythonV](https://img.shields.io/pypi/pyversions/aaindex?logo=2)](https://pypi.org/project/aaindex/)
[![Platforms](https://img.shields.io/badge/platforms-linux%2C%20macOS%2C%20Windows-green)](https://pypi.org/project/aaindex/)
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
  * [Tests](#tests)
  * [Contact](#contact)
  * [License](#license)
  * [References](#References)

Introduction
------------
The AAindex is a database of numerical indices representing various physicochemical, structural and biochemical properties of amino acids and pairs of amino acids 🧬. The AAindex consists of three sections: AAindex1 for the amino acid index of 20 numerical values, AAindex2 for the amino acid mutation matrix and AAindex3 for the statistical protein contact potentials. All data are derived from published literature [[1]](#references). 

This `aaindex` Python software package is a very lightweight way of accessing the data represented in the various AAindex databases, requiring no additional external library installations. Any record within the 3 databases and their associated data/numerical indices can be accessed in one simple command. Currently the software supports the AAindex1 database with plans to include the AAindex 2 & 3 in the future.

<strong>A demo of the software is available [here](https://colab.research.google.com/drive/1dccV_n1BRMiU8W13F9PPXbSaFzvOdQLC?usp=sharing). </strong>

<!-- ### Format of AAindex1 record
![alt text](https://raw.githubusercontent.com/amckenna41/aaindex/main/images/aaindex_example.png)

```
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
``` -->

Installation
-----------------
Install the latest version of `aaindex` using pip:

```bash
pip3 install aaindex --upgrade
```

Install by cloning the repository:

```bash
git clone https://github.com/amckenna41/aaindex.git
python3 setup.py install
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
{'category': 'sec_struct', 'correlation_coefficients': {}, 
'description': 'Normalized frequency of N-terminal non helical region (Chou-Fasman, 1978b)', 'notes': '', 'pmid': '364941', 
'references': "Chou, P.Y. and Fasman, G.D. 'Prediction of the secondary structure of proteins from their amino acid sequence' Adv. Enzymol. 47, 45-148 (1978)", 'values': {'-': 0, 'A': 0.7, 'C': 0.65, 'D': 0.98, 'E': 1.04, 'F': 0.93, 'G': 1.41, 'H': 1.22, 'I': 0.78, 'K': 1.01, 'L': 0.85, 'M': 0.83, 'N': 1.42, 'P': 1.1, 'Q': 0.75, 'R': 0.34, 'S': 1.55, 'T': 1.09, 'V': 0.75, 'W': 0.62, 'Y': 0.99}}
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

record_correlation_coefficient = aaindex1['CHOP780206']['correlation_coefficient']
record_correlation_coefficient = aaindex1['CHOP780206'].correlation_coefficient
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

## AAIndex2 Usage
```python
# from aaindex import aaindex1
from aaindex import aaindex2
# from aaindex import aaindex3
```
## AAIndex3 Usage
```python
# from aaindex import aaindex1
# from aaindex import aaindex2
from aaindex import aaindex3
```

Directories 📁
--------------
* `/tests` - unit and integration tests for `aaindex` package.
* `/aaindex` - source code and all required external data files for package.
* `/images` - images used throughout README.
* `/docs` - `aaindex` documentation.
 
Tests 🧪
--------
To run all tests, from the main `aaindex` folder run:
```
python3 -m unittest discover tests
```

Contact ✉️
---------
If you have any questions or comments, please contact amckenna41@qub.ac.uk or raise an issue on the [Issues][Issues] tab.

License
-------
Distributed under the MIT License. See `LICENSE` for more details.  

References
----------
\[1\]: Shuichi Kawashima, Minoru Kanehisa, AAindex: Amino Acid index database, Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Page 374, https://doi.org/10.1093/nar/28.1.374 <br>
\[2\]: https://www.genome.jp/aaindex/ 

[Back to top](#TOP)

[python]: https://www.python.org/downloads/release/python-360/
[aaindex]: https://github.com/amckenna41/aaindex
[requests]: https://requests.readthedocs.io/en/latest/
[numpy]: https://numpy.org/
[PyPi]: https://pypi.org/project/aaindex/
[demo]: https://colab.research.google.com/drive/1dccV_n1BRMiU8W13F9PPXbSaFzvOdQLC?usp=sharing
[Issues]: https://github.com/amckenna41/aaindex/issues