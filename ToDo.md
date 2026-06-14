- [ ] Visualise function in class where you can do histogram/heat map etc of each amino acid and corresponding value.
- [ ] Use aaindex API URLs search_url = 'https://www.genome.jp/dbget-bin/www_bfind_sub?locale=en&serv=gn&keywords=charge&page=1&max_hit=0&dbkey=aaindex' , record_url = 'https://www.genome.jp/dbget-bin/www_bget?aaindex:KLEP840101'
- [ ] Pull request to biopython for aaindex.
- [X] Add readthedocs badge.
- [ ] Add emojis to readme - test if they show up on Mac and Windows.
- [X] Add documentation section in readme.
- [X] Test last updated - self.last_updated = "February 13, 2017" #as shown on homepage (https://www.genome.jp/aaindex/)
- [X] Remove one-hot encoding from aaindex, and sci-kit learn library. 
- [X] Below aaindex record image in readme - add description of each record key.
- [X] Unit tests to check that aaindex1 can be accessed using dot notation. 
- [X] Double check usage in aaindex/readme.md
- [ ] Remove all camel casing function names/vars, change to underscores and lowercase (https://peps.python.org/pep-0008/#function-and-variable-names).
- [X] Update dependancies so non-specific pandas/numpy version required.
- [X] Think about removing download_using and download function all together.
- [X] Add and fix codecov and coverage installs from github workflow.
- [X] Update workflows.
- [X] Can aaindex1 object accept lowercase version of index - add unit tests of lowercase indices.
- [X] Remove delayed==0.11.0b1 from requirments.txt.
- [X] Update setup.cfg.
- [X] Change all comment underlining from "------" to "=======".
- [X] Add methods comments to aaindex class and test classes.
- [X] Mention its lightweight and requires no external python libraries to use.
- [X] Change all refs of AAIndex to AAindex.
- [X] Go over all modules etc.
- [X] Go over main readme.
- [X] Read over code.
- [X] Read over tests.
- [X] Change physiochemical to physicochemical.
- [X] In /docs folder, add readme with process of how to build and upload the docs to the
- [ ] https://github.com/tadorfer/AAIndex/blob/master/raw_to_df.ipynb
- [X] Parameter typing https://docs.python.org/3/library/typing.html
- [ ] Update workflows
- [ ] Double check spelling using code-spell-checker extension

----- Done
bug / fixMissing value (NA) handling is inconsistent

Several AAindex1 entries use NA for amino acids with no measured value (e.g. Pro in helix propensity scales). These are currently stored as the string "NA" in the JSON. Arithmetic on those values silently produces TypeError at runtime with no guidance to the caller.

bug / fixsearch() searches description only — misses accession codes and PMID

Searching for "364941" (a known PMID) or for an exact accession number like "CHOP780206" returns an empty dict, even though those are valid record identifiers. Search should cover all text fields.

bug / fixget_record_by_category() only exists on AAindex1

AAindex2 and AAindex3 records have no category field, and the method is absent from the shared _AAIndexMatrix base. This is fine, but it is undocumented — the README usage sections for AAindex2/3 do not mention this limitation, leading users to expect parity.

bug / fixAAindex2 non-symmetric matrices silently return wrong values for transposed pairs

The 19 non-symmetric matrices store different scores for (A→R) vs (R→A). The current get() implementation always looks up matrix[aa1][aa2] without checking matrix symmetry, so the transpose lookup is wrong for non-symmetric records.

----- Done

featureSequence encoding / feature-vector generation

The most common downstream use-case for AAindex1 is encoding a protein sequence into a numeric vector for ML pipelines. One convenience method would cover the vast majority of users:

featureCorrelation graph / similar-index lookup

Every AAindex1 record already stores correlation_coefficients linking it to similar indices (|r| ≥ 0.8). Expose a traversal helper so users can find a cluster of related properties without manually following chains.

featureExport to dict, JSON, and optionally pandas DataFrame

Allow exporting a record's values (or the full database) to common interchange formats. The pandas path should be optional (guarded by try/import) to preserve the zero-dependency promise.

featurePairwise comparison / distance matrix between indices

Researchers often want to know how redundant two indices are — i.e. their Pearson r — without reading the pre-stored correlation table (which only covers |r| ≥ 0.8).

----- Done

DXImprove __repr__ and __str__ for records

Currently a record prints as a raw dict. A structured repr makes interactive (REPL) work far more pleasant:

DXFuzzy / partial search with ranking

Current search() is exact substring match. A ranked fuzzy search (difflib or simple token scoring, zero deps) would help users find records when they don't remember the exact wording:

DXCategory enumeration / validation

Valid category strings are not documented or enumerable. A user who passes 'hydrophobic' instead of 'hydrophobicity' gets an empty dict with no feedback.

DXRicher error messages on invalid accession codes

Passing an unrecognised key raises a generic KeyError. A custom exception with a hint (did-you-mean, correct DB suggestion) would save significant debugging time.

DXComplete the camelCase → snake_case rename (open ToDo item)

The ToDo lists this as an open task. Any remaining camelCase public names (if still present after v1.2.0) should be deprecated with a shim and removed in the next major version. The deprecation warnings make migration gradual and auditable.

perfLazy loading — don't parse all 566 records on import

The module-level singleton parses and loads the full JSON on first import. For scripts that only ever access 2–3 records, this is wasteful. A lazy-loading approach (parse on first __getitem__ or __iter__) would cut cold-start time measurably.

API designparse_aaindex() and parse_categories() are public but should be private

Both methods appear in the public API docs and __all__, but they're implementation details — called internally on first load and never needed by callers. Exposing them invites misuse (calling them manually re-parses and overwrites the JSON cache). They should be prefixed _parse_aaindex() and _parse_categories() and removed from __all__. A v1 → v2 deprecation cycle with a warning can handle the rename cleanly.
----