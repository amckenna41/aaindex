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

HERE = pathlib.Path(__file__).parent

README = (HERE / 'README.md').read_text()

setup(name='aaindex',
      version=aaindex.__version__,
      description='A Python package ',
      long_description = README,
      long_description_content_type = "text/markdown",
      author=aaindex.__license__,
      author_email=aaindex.__authorEmail__,
      license=aaindex.__license__,
      url=aaindex.__url__,
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Mathematics'
    ],

      install_requires=[
          'numpy>=1.16.6',
          'delayed',
          'scikit-learn==0.24.1',
          'requests>=2.25',
          'urllib3>=1.26'
      ],

     packages=find_packages(),
     include_package_data=True,
     zip_safe=True)
