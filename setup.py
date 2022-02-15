#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

long_description = """PSDretrieval retrieves the particle size distribution from multi-frequency radar Doppler spectra.
"""

setup(name='PSDretrieval',
    description='PSD retrieval',
    author='Markus Karrer',
    author_email='karrer.markus@web.de',
    url='https://github.com/markuskarrer/PSDretrieval',
    packages= ['PSDretrieval'],
    package_data = {
        '': ['sample_data'],
    },
    long_description = long_description,
    license = 'GPL',
    include_package_data=True
     )


