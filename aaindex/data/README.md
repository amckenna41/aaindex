# Data files required for aaindex package

* `aaindex1` - numerical amino acid indices for aaindex1 database.
* `aaindex2` - amino acid mutation matrices for aaindex2 database.
* `aaindex3` - amino acid pair-wise contact potentials for aaindex3 database.
* `aaindex1.json` - aaindex1 database in parsed JSON format.
* `aaindex2.json` - aaindex2 database in parsed JSON format.
* `aaindex3.json` - aaindex3 database in parsed JSON format.
* `aaindex_to_category.txt` - text file matching each numerical index from the aaindex1 to its associated category.

The `.json` files are parsed from the raw databases above and are committed so a fresh
install never has to re-parse them or write into `site-packages`. Regenerate them with:

```bash
rm aaindex/data/*.json && python -c "import aaindex"
```

`tests/test_data_cache.py` fails if they drift from the raw databases.