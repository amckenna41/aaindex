Changelog
=========

The full changelog is maintained in the repository root:
`CHANGELOG.md <https://github.com/amckenna41/aaindex/blob/main/CHANGELOG.md>`_.

Latest: v1.2.0
---------------

* Added ``_AAIndexMatrix`` shared base class — eliminates code duplication between AAIndex2/3.
* Added ``Map`` class in ``_aaindex_matrix.py`` — single source for dot-notation dict wrapper.
* Added ``values()``, ``plot()``, type annotations, and Google-style docstrings across all modules.
* Migrated build configuration to ``pyproject.toml``; removed ``setup.cfg`` and ``setup.py``.
* AAIndex2 and AAIndex3 now inherit from ``_AAIndexMatrix`` (~30 lines each).
* Added edge-case, symmetry, and values tests.
* See `CHANGELOG.md <https://github.com/amckenna41/aaindex/blob/main/CHANGELOG.md>`_ for the complete history.
