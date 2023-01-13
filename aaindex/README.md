
AAIndex1 Usage
--------------

The AAIndex module is made up of the AAIndex class which itself has all the functions/attributes of the package, so when importing the module you have to import the class as well.

```python
from aaindex import aaIndex
```

### Get record from AAIndex module
The AAIndex class offers diverse functionalities for obtaining any element from any record in the database. Each record is stored in json format in a class attribute called <em>aaindex_json</em>. You can search for a particular record by its index/record code, description or reference. You can also get the index category, and importantly its associated amino acid values.

```python
from aaindex import aaIndex

full_record = aaIndex['CHOP780206']   #get full AAI record
''' Above statement will return -> 
{'description': 'Normalized frequency of N-terminal non helical region (Chou-Fasman, 1978b)', 'notes': '', 'refs': "Chou, P.Y. and Fasman, G.D. 'Prediction of the secondary structure of proteins from their amino acid sequence' Adv. Enzymol. 47, 45-148 (1978); Kawashima, S. and Kanehisa, M.                     'AAindex: amino acid index database.'  Nucleic Acids Res. 28, 374 (2000).", 'values': {'-': 0, 'A': 0.7, 'C': 0.65, 'D': 0.98, 'E': 1.04, 'F': 0.93, 'G': 1.41, 'H': 1.22, 'I': 0.78, 'K': 1.01, 'L': 0.85, 'M': 0.83, 'N': 1.42, 'P': 1.1, 'Q': 0.75, 'R': 0.34, 'S': 1.55, 'T': 1.09, 'V': 0.75, 'W': 0.62, 'Y': 0.99}}
'''
#get individual elements of AAIndex record
record_values = aaIndex['CHOP780206']['values']
record_description = aaIndex['CHOP780206']['description']
record_references = aaIndex['CHOP780206']['refs']

```
### Get category from AAIndex record 
```python

"""
Categories: 
Each AAIndex record is classified into 1 of 8 categories: Charge, Composition, Flexibility, Geometry, Hydrophobic, Meta, 
Observable, Polar and Secondary Structure. The record categories are parsed from the aaindex_categories.txt file and can be accessed for each record via the get_category_from_record() function.
"""
category = aaIndex.get_category_from_record('CHOP780206')

```

### Get total number of AAIndex records
```python
#get total number of records in AAI database
print(aaIndex.get_num_records())

```

### Get list of all AAIndex record names
```python
#get list of all AAIndex record names
print(aaIndex.get_record_names())

```

AAIndex2 Usage
--------------

AAIndex3 Usage
--------------