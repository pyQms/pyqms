#!/usr/bin/env python

from setuptools import setup
import os

version_txt_path = os.path.join(
    os.path.dirname(__file__),
    'pyqms',
    'version.txt'
)
with open(version_txt_path, 'r') as version_file:
    pyqms_version = version_file.readline().strip()

setup(
    name             = 'pyqms',
    version          = pyqms_version,
    packages         = [ 'pyqms' ],
    package_dir      = { 'pyqms' : 'pyqms' },
    description      = 'pyQms - Python module for accurate quantification of mass spectrometry data',
    package_data     = {
        'pyqms' : [
            'version.txt',
            'kb/ext/unimod.xml',
        ]
    },
    install_requires         = ['pymzml', 'openpyxl'],
    long_description = "pyQms enables universal and accurate quantification of mass spectrometry data",
    author           = 'Johannes Leufken, Anna Niehues, L. Peter Sarin, Florian Wessels, Michael Hippler, Sebastian A. Leidel and Christian Fufezan',
    author_email     = 'christian@fufezan.net',
    url              = 'https://github.com/pyQms/pyqms',
    license          = 'The MIT License (MIT)',
    platforms        = 'any that supports python 3.4',
    classifiers      = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: SunOS/Solaris',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
