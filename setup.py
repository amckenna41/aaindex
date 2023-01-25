###############################################################################
#####   Setup.py - installs all the required packages and dependancies    #####
###############################################################################

import pathlib
from setuptools import setup, find_packages
import sys
import aaindex 

#ensure python version is greater than 3
if sys.version_info[0] < 3:
    sys.exit('Python 3 is the minimum version requirement')

#parse readme file
HERE = pathlib.Path(__file__).parent
README = (HERE / 'README.md').read_text()

setup(name=aaindex.__name__,
      version=aaindex.__version__,
      description=aaindex.__description__,
      long_description = README,
      long_description_content_type = "text/markdown",
      author=aaindex.__author__,
      author_email=aaindex.__authorEmail__,
      maintainer=aaindex.__maintainer__,
      license=aaindex.__license__,
      url=aaindex.__url__,
      download_url=aaindex.__download_url__,
      keywords=aaindex.__keywords__,
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
      install_requires=[
          'numpy>=1.16.6',
        #   'delayed',
          'requests>=2.25',
          'urllib3>=1.26'
      ],

     packages=find_packages(),
     include_package_data=True,
     zip_safe=True)