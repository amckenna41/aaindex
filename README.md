## Python package for working with AAIndex database (https://www.genome.jp/aaindex/) <a name="TOP"></a>
<p align="center">
  <img src="https://raw.githubusercontent.com/amckenna41/aaindex/main/images/aaindex_logo.png" />
</p>

[![AAIndex](https://img.shields.io/pypi/v/aaindex)](https://pypi.org/project/aaindex/)
[![pytest](https://github.com/amckenna41/aaindex/workflows/Building%20and%20Testing%20%F0%9F%90%8D/badge.svg)](https://github.com/amckenna41/aaindex/actions?query=workflowBuilding%20and%20Testing%20%F0%9F%90%8D)
[![Build](https://img.shields.io/github/workflow/status/amckenna41/aaindex/Deploy%20to%20PyPI%20%F0%9F%93%A6)](https://github.com/amckenna41/aaindex/actions)
[![Platforms](https://img.shields.io/badge/platforms-linux%2C%20macOS%2C%20Windows-green)](https://pypi.org/project/aaindex/)
[![PythonV](https://img.shields.io/pypi/pyversions/aaindex?logo=2)](https://pypi.org/project/aaindex/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)
[![Issues](https://img.shields.io/github/issues/amckenna41/aaindex)](https://github.com/amckenna41/aaindex/issues)
[![Size](https://img.shields.io/github/repo-size/amckenna41/aaindex)](https://github.com/amckenna41/aaindex)
[![Commits](https://img.shields.io/github/commit-activity/w/amckenna41/aaindex)](https://github.com/amckenna41/aaindex)

Table of Contents
-----------------

  * [Introduction](#introduction)
  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Tests](#tests)
  * [Contact](#contact)
  * [References](#References)

Introduction
------------
As stated on the AAIndex website - The AAindex is a database of numerical indices representing various physicochemical and biochemical properties of amino acids and pairs of amino acids. AAindex consists of three sections now: AAindex1 for the amino acid index of 20 numerical values, AAindex2 for the amino acid mutation matrix and AAindex3 for the statistical protein contact potentials. All data are derived from published literature [[1]](#references). 

This aaindex Python software package is a very lightweight way of accessing the data represented in the various AAIndex databases. Minimal requirements and external libraries are required to use the package and any record and its associated data can be accessed in one line. Currently the software supports the AAIndex1 database with plans to include AAIndex 2 & 3 in the future. The format of an AAIndex1 record can be seen below.
### Format of AAIndex record
![alt text](https://raw.githubusercontent.com/amckenna41/aaindex/main/images/aaindex_example.png)

Requirements
------------
* [Python][python] >= 3.6
* [numpy][numpy] >= 1.16.0
* [sklearn][sklearn] >= 0.24
* [requests][requests] >= 2.24.0

Installation
-----------------
Install the latest version of aaindex using pip:

```bash
pip3 install aaindex
```

Install by cloning repository:

```bash
git clone https://github.com/amckenna41/aaindex.git
python3 setup.py install
```
Usage
-----

The AAIndex module is made up of the AAIndex class which itself has all the functions/attributes of the package, so when importing the module you have to import the class as well.

```python
from aaindex.aaindex import aaindex

```

### Get record from AAIndex module
The AAIndex class offers diverse functionalities for obtaining any element from any record in the database. Each record is stored in json format in a class attribute called <em>aaindex_json</em>. You can search for a particular record by its index/record code, description or reference. You can also get the index category, and importantly its associated amino acid values.

```python
from aaindex.aaindex import aaindex

full_record = aaindex['CHOP780206']   #get full AAI record
''' Above statement will return -> 
{'description': 'Normalized frequency of N-terminal non helical region (Chou-Fasman, 1978b)', 'notes': '', 'refs': "Chou, P.Y. and Fasman, G.D. 'Prediction of the secondary structure of proteins from their amino acid sequence' Adv. Enzymol. 47, 45-148 (1978); Kawashima, S. and Kanehisa, M.                     'AAindex: amino acid index database.'  Nucleic Acids Res. 28, 374 (2000).", 'values': {'-': 0, 'A': 0.7, 'C': 0.65, 'D': 0.98, 'E': 1.04, 'F': 0.93, 'G': 1.41, 'H': 1.22, 'I': 0.78, 'K': 1.01, 'L': 0.85, 'M': 0.83, 'N': 1.42, 'P': 1.1, 'Q': 0.75, 'R': 0.34, 'S': 1.55, 'T': 1.09, 'V': 0.75, 'W': 0.62, 'Y': 0.99}}
'''
#get individual elements of AAIndex record
record_values = aaindex['CHOP780206']['values']
record_description = aaindex['CHOP780206']['description']
record_references = aaindex['CHOP780206']['refs']

```
### Get category from AAIndex record 
```python

"""
Categories: 
Each AAIndex record is classified into 1 of 8 categories: Charge, Composition, Flexibility, Geometry, Hydrophobic, Meta, 
Observable, Polar and Secondary Structure. The record categories are parsed from the aaindex_categories.txt file and can be accessed for each record via the get_category_from_record() function.
"""
category = aaindex.get_category_from_record('CHOP780206')

```

### Get total number of AAIndex records
```python
#get total number of records in AAI database
print(aaindex.get_num_records())

```

### Get list of all AAIndex record names
```python
#get list of all AAIndex record names
print(aaindex.get_record_names())

```

Directories
-----------
* `/tests` - unit and integration tests for aaindex package.
* `/aaindex` - source code and all required external data files for package.
* `/images` - images used throughout README.

Tests
-----
To run all tests, from the main aaindex folder run:
```
python3 -m unittest discover
```

To run main test module, from the main aaindex folder run:
```
python -m unittest tests.test_aaindex -v
```

Contact
-------
If you have any questions or comments, please contact amckenna41@qub.ac.uk or raise an issue on the [Issues][Issues] tab.

References
----------
\[1\]: https://www.genome.jp/aaindex/

[Back to top](#TOP)

[python]: https://www.python.org/downloads/release/python-360/
[numpy]: https://numpy.org/
[pandas]: https://pandas.pydata.org/
[sklearn]: https://scikit-learn.org/stable/
[requests]: https://docs.python-requests.org/en/latest/
[Issues]: https://github.com/amckenna41/pySAR/issues

