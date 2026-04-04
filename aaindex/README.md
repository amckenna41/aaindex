# aaindex


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

References
----------
\[1\]: Shuichi Kawashima, Minoru Kanehisa, AAindex: Amino Acid index database, Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Page 374, https://doi.org/10.1093/nar/28.1.374 <br>
\[2\]: https://www.genome.jp/aaindex/ 

[readthedocs]: https://aaindex.readthedocs.io/en/latest/