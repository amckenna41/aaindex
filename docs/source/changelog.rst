Changelog
=========

The full changelog is maintained in the repository root:
`CHANGELOG.md <https://github.com/amckenna41/aaindex/blob/main/CHANGELOG.md>`_.

Latest: v1.3.2
---------------

* Fixed: AAindex2 record ``DOSZ010101`` was parsed into an empty matrix. Its source header has a ``cols rows =`` typo, and reading the rows/cols amino acid strings by token position took ``'='`` as the entire column order — the record collapsed to 22×1, every ``get()`` returned ``None`` and it was wrongly flagged ``is_symmetric=True``. All 440 values are restored.
* Fixed: record accessors leaked references into the loaded database. ``__getitem__``, ``to_dict()``, ``values()``, ``search()``, ``search_fuzzy()``, ``get_record_by_category()`` and ``amino_acids()`` now return copies, so mutating a returned record no longer changes what later callers see.
* Fixed: a truncated JSON cache permanently broke the package. Cache writes are now atomic, and an unreadable cache falls back to parsing the raw database with a ``UserWarning``.
* Fixed: a malformed I-block raised ``IndexError`` and a record missing from the category file raised ``KeyError``; both now warn and degrade gracefully (uncategorised records get category ``'unknown'``).
* Fixed: ``code in db`` now normalises case and whitespace, so it agrees with ``db[code]``.
* Fixed: ``AAIndex1.search()`` also matches the ``category`` field, consistent with ``search_fuzzy()``.
* Changed: ``aaindex2.json`` and ``aaindex3.json`` are now shipped with the package alongside ``aaindex1.json``, so a fresh install never re-parses the raw databases or writes into ``site-packages``.
* Changed: removed the redundant ``aaindex_categories.txt`` data file (byte-identical to the comment-stripped ``aaindex_to_category.txt``).

v1.3.1
------

* Added ``AAIndex1.multi_encode(sequence, record_codes)`` — encode one sequence against multiple records in a single call.
* Added ``unknown_value`` parameter to ``AAIndex1.encode()`` and ``batch_encode()``, distinguishing non-canonical amino acids from canonical ones with NA source data.
* Added ``include_gap`` parameter to ``AAIndex1.to_dataframe()`` for a clean 20-column DataFrame.
* Added ``AAIndex1.plot_distribution()`` and ``_AAIndexMatrix.plot_heatmap()`` (require ``matplotlib``).
* Fixed: ``AAIndex1.get_all_categories(category_file)`` no longer ignores its argument.
* Fixed: unknown field codes in a raw database file emit a ``UserWarning`` instead of raising ``KeyError``.
* Changed: ``search()`` returns results sorted alphabetically by accession code.

v1.3.0
------

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
