# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- `Map.__getattr__` now raises `AttributeError` for missing keys instead of silently returning `None` — fixes `hasattr()` behaviour and prevents silent bugs.
- `.readthedocs.yml` Sphinx configuration path corrected from `docs/conf.py` to `docs/source/conf.py`.
- Removed unused `Optional` import from `aaindex1.py`.
- `MANIFEST.in` glob tightened from `graft aaindex/data*` to `graft aaindex/data`.

### Added
- Python 3.11 and 3.12 added to CI test matrix in `build_test.yml`.
- `CODECOV_TOKEN` secret reference added to codecov upload step.
- Single-source versioning via `importlib.metadata.version()` in `__init__.py`.
- `__all__` in `__init__.py` for explicit public API.
- `aaindex/py.typed` PEP 561 marker for type-checker compatibility.
- `[tool.ruff]` and `[tool.ruff.lint]` configuration in `pyproject.toml`.
- `docs/requirements.txt` for Read the Docs builds.
- `tests/conftest.py` with shared `aai1`, `aai2`, `aai3` fixtures.
- `test_plot`, `test_plot_invalid_record`, and `test_repr_format` tests for AAIndex1.
- Changelog page added to Sphinx documentation.
- `.pre-commit-config.yaml` with ruff and pre-commit-hooks.
- `sphinx-rtd-theme` added to `[project.optional-dependencies.docs]`.
- `aaindex/data/README.md` updated to list `aaindex2.json` and `aaindex3.json`.
- JSON cache files (`aaindex/data/*.json`) added to `.gitignore`.

### Changed
- `num_records()` in `AAIndex1` now uses `len(self.aaindex_json)` instead of sorting all keys.
- Deploy workflows (`deploy_pypi.yml`, `deploy_test_pypi.yml`) upgraded from Python 3.10 to 3.12.
- `TODO.md` trimmed to only remaining actionable items.
- `aaindex/README.md` simplified to a one-liner pointing at the project README.
- CircleCI badge removed from root `README.md`.

### Removed
- `.circleci/` directory (stale, referenced deleted `setup.py`).
- `old_stuff/` directory (legacy archived code).
- `aaindex.egg-info/` build artifact.

## [1.2.0]

### Added
- `_AAIndexMatrix` shared base class in `aaindex/_aaindex_matrix.py` — eliminates code duplication between AAIndex2 and AAIndex3 by centralising parsing, search, lookup, `values()`, `plot()`, and dunder methods.
- Canonical `Map` class in `_aaindex_matrix.py` — single source for the dot-notation dict wrapper, imported by all three modules.
- `__all__` exported from every module (`aaindex1`, `aaindex2`, `aaindex3`, `_aaindex_matrix`).
- `values(record_code)` method on AAIndex2 and AAIndex3 for retrieving the full matrix dict.
- `plot(record_code)` method on AAIndex1 (bar chart), AAIndex2, and AAIndex3 (heatmap) for quick visual inspection.
- Type annotations on all public methods and properties across all modules.
- Google-style docstrings on all classes and methods.
- Python 3.11 and 3.12 classifiers in project metadata.
- `pyproject.toml` as the single build/metadata configuration file.
- Edge-case tests: non-canonical amino acid "Z" in `get()`, empty-string search, `TypeError` for non-string `__getitem__` input, whitespace-padded record codes.
- Explicit full-matrix symmetry tests for AAIndex2 and AAIndex3.
- `values()` tests for AAIndex2 and AAIndex3.

### Changed
- `AAIndex2` and `AAIndex3` now inherit from `_AAIndexMatrix`; module files reduced from ~500+ lines each to ~30 lines.
- `AAIndex1` imports `Map` from `_aaindex_matrix` instead of defining its own copy.
- All `.format()` string formatting replaced with f-strings across the codebase.
- `__init__.py` cleaned up: removed `__name__` override, camelCase attributes (`__authorEmail__`), and redundant metadata (`__download_url__`, `__status__`, `__keywords__`, `__test_suite__`). Only `__version__`, `__author__`, and `__license__` remain.
- README updated: removed stale "plans to include AAindex 2 & 3" text, added full AAIndex2 and AAIndex3 usage examples, updated install instructions from `setup.py` to `pip install .`.

### Removed
- `setup.cfg` and `setup.py` — replaced by `pyproject.toml`.
- Duplicate `Map` class definitions from `aaindex1.py`, `aaindex2.py`, and `aaindex3.py`.

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

[Unreleased]: https://github.com/amckenna41/aaindex/compare/v1.1.2...HEAD
[1.1.2]: https://github.com/amckenna41/aaindex/compare/v1.1.0...v1.1.2
[1.1.0]: https://github.com/amckenna41/aaindex/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/amckenna41/aaindex/compare/v0.0.1...v1.0.0
[0.0.1]: https://github.com/amckenna41/aaindex/releases/tag/v0.0.1
