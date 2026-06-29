# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.1]

### Added
- `AAIndex1.multi_encode(sequence, record_codes, gap_value=0.0, unknown_value=None)` — encode a single protein sequence against multiple AAindex1 records in one call. Returns a dict mapping each accession number to its encoded vector. Avoids repeated record lookups and is the recommended approach for multi-index feature extraction in ML pipelines.
- `AAIndex1.encode()` and `AAIndex1.batch_encode()` — new `unknown_value` parameter (default `None`). Non-canonical amino acids (e.g. `B`, `X`, `Z`) are now distinguishable from canonical AAs with NA source data: NA-valued canonical AAs remain `None`; non-canonical AAs receive `unknown_value`.
- `AAIndex1.to_dataframe()` — new `include_gap` parameter (default `True`). Set `include_gap=False` to drop the `'-'` gap-placeholder column, producing a clean 20-column DataFrame suitable for ML pipelines.
- `AAIndex1.plot_distribution(record_code)` — bar chart of amino acid values for a given AAindex1 record. Requires `matplotlib`; raises `ImportError` with install instructions if absent. Returns the `matplotlib.axes.Axes` object for further customisation.
- `_AAIndexMatrix.plot_heatmap(record_code)` — 20×20 heatmap of the pairwise matrix for a given AAindex2 or AAindex3 record. NA values are rendered as masked (white) cells. Requires `matplotlib`. Returns the `Axes` object.

### Fixed
- `AAIndex1.get_all_categories(category_file)` — the `category_file` parameter was previously accepted but silently ignored; it is now forwarded to `_parse_categories()`. Default changed from `"aaindex_categories.txt"` to `"aaindex_to_category.txt"` to match `_parse_categories()` and preserve existing behaviour.
- `_parse_aaindex()` on all three classes — added guard against unknown field codes. Previously, a line starting with an unrecognised character would raise `KeyError` when the parser tried to append to an unknown bucket. Such lines now emit a `UserWarning` and are skipped gracefully.
- Test expected value dicts for the gap placeholder `'-'` corrected from `int` `0` to `float` `0.0`. While `0 == 0.0` in Python, the fixtures are now type-consistent with `test_gap_placeholder_is_float` and the actual stored value.
- `deploy_test_pypi.yml` — removed duplicate `--index-url` / `-i` flags in the TestPyPI install-and-verify step.
- `deploy_pypi.yml` — the "Verify package version" step now reads `SOURCE_VERSION` from `pyproject.toml` via `tomllib` (Python 3.11+ stdlib) instead of importing the just-installed PyPI package. This makes the comparison meaningful: local source version vs. what was installed from PyPI.

### Changed
- `search()` on all three database classes now returns results **sorted alphabetically by accession code** instead of in database insertion order. This makes the output deterministic and reproducible across Python versions.
- Type annotations modernised to Python 3.10+ union syntax throughout `aaindex1.py` and `_aaindex_matrix.py`: `Optional[X]` → `X | None`, `List[X]` → `list[X]`, `Dict` → `dict`, `Union[X, Y]` → `X | Y`. The `from typing import` block is reduced to only `Any` and `Iterator`.
- `pyproject.toml` — `[tool.ruff] target-version` corrected from `"py38"` to `"py310"` to match `requires-python = ">=3.10"`. Ruff will now correctly apply Python 3.10+ lint and upgrade rules.

## [1.3.0]

### Added
- `search_fuzzy(query, n=10, cutoff=0.0)` on all three databases — fuzzy/partial search using `difflib.SequenceMatcher`. Results are ordered by descending similarity score (most relevant first) and limited to the top `n` matches. Accepts a `cutoff` threshold in [0, 1].
- `AAIndex1.valid_categories()` — returns a sorted list of all distinct category strings present in the database.
- `AAIndex1.batch_encode(sequences, record_code, gap_value=0.0)` — encodes a list of protein sequences as numeric feature vectors using a single index in one call. The record is looked up once and reused across all sequences, making this more efficient than calling `encode()` in a loop.
- `version` attribute on all three database class instances (`aaindex1.version`, `aaindex2.version`, `aaindex3.version`) — exposes the installed package version string without needing to import `importlib.metadata` separately.
- `AAIndex1.encode(sequence, record_code)` — encodes a protein sequence as a numeric feature vector using a named index. NA positions are returned as `None`; gaps use a configurable `gap_value` (default `0.0`).
- `AAIndex1.get_correlated_indices(record_code, min_correlation=0.8, depth=1)` — traverses the stored `correlation_coefficients` chains to return a dict of related index codes and their Pearson r values, up to `depth` hops from the seed.
- `AAIndex1.compare_indices(code1, code2)` — computes the Pearson correlation coefficient between two AAindex1 records using their amino acid values (NA positions excluded pairwise). Returns `None` when insufficient non-NA values are available.
- `to_dict(record_code=None)` on all three databases — exports one record or the full database as a plain Python dict.
- `to_json(record_code=None, indent=4)` on all three databases — serialises one record or the full database to a JSON string.
- `to_dataframe(record_code=None)` on all three databases — exports amino acid values (AAindex1) or pairwise matrix scores (AAindex2/3) as a `pandas.DataFrame`. pandas is optional; a clear `ImportError` is raised if it is not installed.

### Fixed
- `AAIndex1`: amino acids with no measured value (NA in source data) are now stored as `None` instead of `0.0`, preventing silent arithmetic errors. Arithmetic on `None` now raises `TypeError` as expected.
- `search()` on all three databases now searches all text fields — description, accession code, PMID, references, and notes — instead of description only.
- `AAIndex2`/`AAIndex3`: records with full 20×20 matrices (non-symmetric) are now parsed correctly. Previously, symmetric fill overwrote correct directional values. Each record now exposes an `is_symmetric` boolean field.
- `correlation_coefficients` values are now stored and returned as `float` (previously stored as raw strings from the source file). Existing JSON caches are normalised to float on load so no re-parse is required.
- Gap placeholder `'-'` in `AAIndex1` values dicts is now stored as `float` `0.0` instead of `int` `0`, making it type-consistent with all other numeric values. Existing caches are normalised on load.
- `_parse_categories()` no longer writes an intermediate `aaindex_categories.txt` file to the package data directory. Filtering of comment lines is now done in memory, fixing failures on read-only installs (system Python, Docker, virtualenvs without writable site-packages).
- `_parse_categories()` now rejects relative paths with directory components (e.g. `../../etc/passwd`) to prevent path traversal; only bare filenames or absolute paths are accepted.
- `data_dir` and `aaindex_filename` property setters on all three database classes now raise `ValueError` if the value contains path separators, preventing path traversal via these attributes.
- `AAIndex1._parse_aaindex()`: added guard against I-block with fewer than 3 lines; malformed records now emit a `UserWarning` and are skipped rather than raising `IndexError`.
- `_parse_aaindex()` on `AAIndex2`/`AAIndex3`: added bounds check on the rows/cols header (`len(parts) >= 6`); malformed headers emit a `UserWarning` instead of crashing with `IndexError`.
- `_parse_aaindex()` on `AAIndex2`/`AAIndex3`: incomplete matrix rows or columns now emit a `UserWarning` instead of silently producing a partial matrix with no indication of data loss.
- JSON cache write in `_parse_aaindex()` on all three classes is now wrapped in `try/except OSError` — a read-only data directory no longer prevents the parsed data from being returned and used.
- `_AAIndexMatrix.amino_acids()` result is now cached after the first call, matching the behaviour of `AAIndex1.amino_acids()` and avoiding repeated list allocation on repeated calls.
- `Map.__delitem__` now uses `dict.pop(key, None)` to remove the key from `__dict__`, preventing `KeyError` when a key is present in the underlying dict but was never routed through `__setitem__`.

### Changed
- `get_record_by_category()` now raises `ValueError` (with a did-you-mean hint) when an unknown category is passed, instead of silently returning an empty dict.
- `__getitem__` on all three databases now includes a did-you-mean hint in the `ValueError` message when an accession code is not found and a close match exists.
- Database parsing is now deferred (lazy loading): `AAIndex1`, `AAIndex2`, and `AAIndex3` no longer parse or load JSON on import. Data is loaded on first access, reducing cold-start time.
- `parse_aaindex()` and `parse_categories()` are now private methods (`_parse_aaindex()` and `_parse_categories()`). The old public names are retained as deprecated shims that emit `DeprecationWarning`.
- `Map.__repr__` now returns a structured `AAIndexRecord(...)` representation with metadata fields listed individually and large nested dicts (values, matrix) summarised as `{...} (N entries)`. `Map.__str__` is aliased to the same output.
- Verified: all public method names are already snake_case; no camelCase renames were required.
- `get_record_by_category()` is an `AAIndex1`-only method. `AAIndex2` and `AAIndex3` records have no category field and therefore do not support this method.

## [1.2.0]

### Fixed
- `Map.__getattr__` now raises `AttributeError` for missing keys instead of silently returning `None` — fixes `hasattr()` behaviour and prevents silent bugs.
- `.readthedocs.yml` Sphinx configuration path corrected from `docs/conf.py` to `docs/source/conf.py`.
- Removed unused `Optional` import from `aaindex1.py`.
- `MANIFEST.in` glob tightened from `graft aaindex/data*` to `graft aaindex/data`.

### Added
- `_AAIndexMatrix` shared base class in `aaindex/_aaindex_matrix.py` — eliminates code duplication between AAIndex2 and AAIndex3 by centralising parsing, search, lookup, `values()`, and dunder methods.
- Canonical `Map` class in `_aaindex_matrix.py` — single source for the dot-notation dict wrapper, imported by all three modules.
- `__all__` exported from every module (`aaindex1`, `aaindex2`, `aaindex3`, `_aaindex_matrix`).
- `values(record_code)` method on AAIndex2 and AAIndex3 for retrieving the full matrix dict.
- Type annotations on all public methods and properties across all modules.
- Google-style docstrings on all classes and methods.
- Python 3.11 and 3.12 classifiers in project metadata.
- `pyproject.toml` as the single build/metadata configuration file.
- Edge-case tests: non-canonical amino acid "Z" in `get()`, empty-string search, `TypeError` for non-string `__getitem__` input, whitespace-padded record codes.
- Explicit full-matrix symmetry tests for AAIndex2 and AAIndex3.
- `values()` tests for AAIndex2 and AAIndex3.
- Python 3.11 and 3.12 added to CI test matrix in `build_test.yml`.
- `CODECOV_TOKEN` secret reference added to codecov upload step.
- Single-source versioning via `importlib.metadata.version()` in `__init__.py`.
- `__all__` in `__init__.py` for explicit public API.
- `aaindex/py.typed` PEP 561 marker for type-checker compatibility.
- `[tool.ruff]` and `[tool.ruff.lint]` configuration in `pyproject.toml`.
- `docs/requirements.txt` for Read the Docs builds.
- `test_repr_format` test for AAIndex1.
- Changelog page added to Sphinx documentation.
- `.pre-commit-config.yaml` with ruff and pre-commit-hooks.
- `sphinx-rtd-theme` added to `[project.optional-dependencies.docs]`.
- `aaindex/data/README.md` updated to list `aaindex2.json` and `aaindex3.json`.
- JSON cache files (`aaindex/data/*.json`) added to `.gitignore`.

### Changed
- `AAIndex2` and `AAIndex3` now inherit from `_AAIndexMatrix`; module files reduced from ~500+ lines each to ~30 lines.
- `AAIndex1` imports `Map` from `_aaindex_matrix` instead of defining its own copy.
- All `.format()` string formatting replaced with f-strings across the codebase.
- `__init__.py` cleaned up: removed `__name__` override, camelCase attributes (`__authorEmail__`), and redundant metadata (`__download_url__`, `__status__`, `__keywords__`, `__test_suite__`). Only `__version__`, `__author__`, and `__license__` remain.
- README updated: removed stale "plans to include AAindex 2 & 3" text, added full AAIndex2 and AAIndex3 usage examples, updated install instructions from `setup.py` to `pip install .`.
- `num_records()` in `AAIndex1` now uses `len(self.aaindex_json)` instead of sorting all keys.
- Deploy workflows (`deploy_pypi.yml`, `deploy_test_pypi.yml`) upgraded from Python 3.10 to 3.12.
- `TODO.md` trimmed to only remaining actionable items.
- `aaindex/README.md` simplified to a one-liner pointing at the project README.
- CircleCI badge removed from root `README.md`.


### Removed
- `setup.cfg` and `setup.py` — replaced by `pyproject.toml`.
- Duplicate `Map` class definitions from `aaindex1.py`, `aaindex2.py`, and `aaindex3.py`.
- `.circleci/` directory (stale, referenced deleted `setup.py`).
- `old_stuff/` directory (legacy archived code).
- `aaindex.egg-info/` build artifact.

## [1.1.2] - 2023-11-16

### Added
- Unit tests for lowercase record index lookups on `AAIndex1`.
- Code coverage reporting via `codecov` integrated into the CI workflow.
- Method-level comments added to all `AAIndex1` class methods and test classes.

### Changed
- Updated `setup.cfg` with current package metadata and dependencies.
- Updated GitHub Actions workflows with improved coverage and security scanning steps.
- All comment underline separators standardised from `------` to `======`.
- All references to `AAIndex` renamed to `AAindex` for consistency.
- Removed camel-cased variable and function names in favour of lowercase with underscores (PEP 8).

### Fixed
- Various bug fixes identified during testing and code review.
- Removed `delayed==0.11.0b1` pinned pre-release dependency from requirements.
- Removed one-hot encoding logic and unused `scikit-learn` dependency.

## [1.1.0] - 2023-09-18

### Added
- Unit tests verifying the complete set of keys present in each `aaindex_json` record.
- Unit tests for `last_updated` attribute value.
- Unit tests for dot-notation record access.

### Changed
- Improved internal record structure to include all standard AAIndex field keys.
- Updated CI/CD workflows.

### Fixed
- Various bug fixes and code quality improvements.

## [1.0.0] - 2023-01-13

### Added
- Complete set of AAIndex record keys present on all parsed records, including `category`, `correlation_coefficients`, `notes`, `pmid`, and `references`.
- Category field (`category`) added to each `AAIndex1` record output object.
- Search-by-keyword functionality via the `search()` method.
- Download URL added to `setup.cfg`/`setup.py`.
- Keywords metadata added to `setup.cfg`/`setup.py`.

### Changed
- Restructured package directory layout.
- Improved record structure to include full bibliographic metadata.
- Updated GitHub Actions workflows.
- Removed `get_category_from_record` standalone function; category is now accessed directly via `aaindex1['XXXX']['category']`.
- Removed `AAINDEX_FILENAME` constant; filename is now a static class attribute.

### Fixed
- Bug fixes relating to directory structure changes.

## [0.0.1] - 2022-02-19

### Added
- Initial beta release of the `aaindex` Python package.
- `AAIndex1` class with parser for the AAIndex1 amino acid property index database.
- `parse_aaindex()` method to read and convert the flat `aaindex1` data file into JSON.
- `parse_categories()` and `get_all_categories()` methods for category lookups.
- `record_codes()`, `record_names()`, `num_records()`, and `amino_acids()` utility methods.
- `__getitem__` support for subscript-style record access (e.g. `aaindex1['CHOP780206']`).
- `__sizeof__` method returning the size of the underlying database file.
- `Map` helper class for dot-notation access to record dicts.
- Module-level `aaindex1` instance.
- `setup.py` and `setup.cfg` for PyPI packaging.
- MIT licence.
- Initial `README.md` with usage examples and database description.
- GitHub Actions workflow for automated build and test.

[1.3.1]: https://github.com/amckenna41/aaindex/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/amckenna41/aaindex/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/amckenna41/aaindex/compare/v1.1.2...v1.2.0
[1.1.2]: https://github.com/amckenna41/aaindex/compare/v1.1.2...v1.2.0
[1.1.2]: https://github.com/amckenna41/aaindex/compare/v1.1.0...v1.1.2
[1.1.0]: https://github.com/amckenna41/aaindex/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/amckenna41/aaindex/compare/v0.0.1...v1.0.0
[0.0.1]: https://github.com/amckenna41/aaindex/releases/tag/v0.0.1
