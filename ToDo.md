- [X] Look into implementing AAIndex2 & 3
- [ ] Visualise function in class where you can do histogram/heat map etc of each amino acid and corresponding value.
- [ ] Implement functionality to work with AAIndex 2 and 3 (https://github.com/Pymol-Scripts/Pymol-script-repo/blob/master/aaindex.py , https://github.com/pycogent/pycogent/blob/master/cogent/parse/aaindex.py)
- [X] Move aaindex into __init__ file, change imports in folders and readme's.
- [X] Get rid of AAINDEX_FILENAME - make aaindex1 name static.
- [ ] Add Contributing section to readme (https://github.com/arc298/instagram-scraper).
- [X] Add references section to readme.
- [ ] Remove ds_store from all folders.
- [X] Ensure all files listed in readme of /data
- [ ] Go through flake8 and bandit output.
- [ ] Search by keyword, other ways to search?
- [ ] Split individual aaindex properties into individual python @property decorator, using a Record class.
- [ ] Use aaindex API URLs search_url = 'https://www.genome.jp/dbget-bin/www_bfind_sub?locale=en&serv=gn&keywords=charge&page=1&max_hit=0&dbkey=aaindex' , record_url = 'https://www.genome.jp/dbget-bin/www_bget?aaindex:KLEP840101'
- [X] Parse main URL to check last updated: "Last updated: February 13, 2017"
- [ ] Pull request to biopython for aaindex.
- [X] Add "python3 -m twine upload --repository testpypi dist/*" to test pypi workflow.
- [X] Add Keywords to setup.py & cfg
- [ ] Add demo using jupyter notebook.
- [X] Update build/test workflow, similar to pySAR.
- [ ] read the docs.
- [ ] Add readthedocs badge - [![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](http://ansicolortags.readthedocs.io/?badge=latest)
- [ ] Add emojis to readme
- [ ] Check variable naming conventions.
- [ ] Change "secrets.PY_PI..." to "secrets.PYPI...".
- [X] Add comments to .circle/workflow.
- [X] Add spacing in between individual references, see if it improves readability, revert if not.
- [ ] Add documentation section in readme.
- [X] Add path-ignore keywords to GitHub Action.
- [X] Reorder software metadata in setup.py to be in order of main func, create __description__ var.
- [X] Add download_url to setup.py - url of zipped package.