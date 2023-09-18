<!-- [![CircleCI](https://circleci.com/gh/amckenna41/aaindex.svg?style=svg&circle-token=d860bb64668be19d44f106841b80eb47a8b7e7e8)](https://app.circleci.com/pipelines/github/amckenna41/aaindex) -->
[![AAindex](https://img.shields.io/pypi/v/aaindex)](https://pypi.org/project/aaindex/)
[![PythonV](https://img.shields.io/pypi/pyversions/aaindex?logo=2)](https://pypi.org/project/aaindex/)
[![pytest](https://github.com/amckenna41/aaindex/workflows/Building%20and%20Testing/badge.svg)](https://github.com/amckenna41/aaindex/actions?query=workflowBuilding%20and%20Testing)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)

*****
AAindex1 Usage
--------------
```python
from aaindex import aaindex1
# from aaindex import aaindex2
# from aaindex import aaindex3
```

### Get record from AAindex1
The AAindex1 class offers diverse functionalities for obtaining any element from any record in the database. The records are imported from a parsed json <em>aaindex_json</em> in the data folder of the package. You can search for a particular record by its index/record code, description or reference. You can also get the index category, and importantly its associated amino acid values:
```python
from aaindex import aaindex1

full_record = aaindex1['CHOP780206']   #get full AAI record
''' full_record ->
{'category': 'sec_struct', 'correlation_coefficients': {}, 'description': 'Normalized frequency of N-terminal non helical region (Chou-Fasman, 1978b)', 'notes': '', 'pmid': '364941', 'references': "Chou, P.Y. and Fasman, G.D. 'Prediction of the secondary structure of proteins from their amino acid sequence' Adv. Enzymol. 47, 45-148 (1978)", 'values': {'-': 0, 'A': 0.7, 'C': 0.65, 'D': 0.98, 'E': 1.04, 'F': 0.93, 'G': 1.41, 'H': 1.22, 'I': 0.78, 'K': 1.01, 'L': 0.85, 'M': 0.83, 'N': 1.42, 'P': 1.1, 'Q': 0.75, 'R': 0.34, 'S': 1.55, 'T': 1.09, 'V': 0.75, 'W': 0.62, 'Y': 0.99}}
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

### Get total number of AAindex records
```python
#get total number of records in AAI database
aaindex1.num_records()
```

### Get list of all AAindex record names
```python
#get list of all AAindex record names
aaindex1.record_names()
```

AAindex2 Usage
--------------
```python
# from aaindex import aaindex1
from aaindex import aaindex2
# from aaindex import aaindex3
```
AAindex3 Usage
--------------
```python
# from aaindex import aaindex1
# from aaindex import aaindex2
from aaindex import aaindex3
```

References
----------
\[1\]: Shuichi Kawashima, Minoru Kanehisa, AAindex: Amino Acid index database, Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Page 374, https://doi.org/10.1093/nar/28.1.374 <br>
\[2\]: https://www.genome.jp/aaindex/ 