Changelog
=========

The full changelog is maintained in the repository root:
`CHANGELOG.md <https://github.com/amckenna41/aaindex/blob/main/CHANGELOG.md>`_.

Latest: v1.3.0
---------------

* Added ``search_fuzzy(query, n=10, cutoff=0.0)`` on all three databases — fuzzy/partial search with ranked results.
* Added ``AAIndex1.valid_categories()`` — returns a sorted list of all distinct category strings.
* Added ``AAIndex1.batch_encode(sequences, record_code)`` — encodes a list of sequences in one call.
* Added ``version`` attribute on all three database instances — exposes the installed package version string directly on the object.
* Added ``AAIndex1.encode()``, ``get_correlated_indices()``, ``compare_indices()``, ``to_dict()``, ``to_json()``, ``to_dataframe()`` across all databases.
* Fixed: NA amino acid values now stored as ``None`` instead of ``0.0``; arithmetic on ``None`` raises ``TypeError`` as expected.
* Fixed: ``search()`` now searches all text fields — description, accession code, PMID, references, and notes.
* Fixed: Non-symmetric 20×20 matrices in ``AAIndex2``/``AAIndex3`` now parsed correctly; each record exposes ``is_symmetric``.
* Fixed: ``correlation_coefficients`` values are now ``float`` (previously raw strings); existing JSON caches are normalised on load.
* Fixed: Gap placeholder ``'-'`` is now ``float`` ``0.0`` instead of ``int`` ``0`` for type consistency.
* Fixed: ``_parse_categories()`` no longer writes an intermediate file, fixing failures on read-only installs; path traversal inputs are now rejected.
* Fixed: ``data_dir`` and ``aaindex_filename`` property setters now reject values with path separators.
* Fixed: Malformed I-blocks and rows/cols headers now emit ``UserWarning`` instead of crashing with ``IndexError``.
* Fixed: JSON cache write is now non-fatal on read-only data directories.
* Fixed: ``amino_acids()`` result cached in ``_AAIndexMatrix`` after first call; ``Map.__delitem__`` uses ``pop`` to avoid ``KeyError``.
* Changed: Database parsing is now deferred (lazy loading); ``parse_aaindex()``/``parse_categories()`` are now private with deprecated public shims; ``Map.__repr__`` returns a structured ``AAIndexRecord(...)`` string.
* See `CHANGELOG.md <https://github.com/amckenna41/aaindex/blob/main/CHANGELOG.md>`_ for the complete history.
