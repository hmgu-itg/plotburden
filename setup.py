#!/usr/bin/env python3
import re

from setuptools import setup, find_packages


with open("plotburden/__init__.py") as f:
    version = re.search(r"__version__ = \'(.*?)\'", f.read()).group(1)

setup(
    name = "plotburden",
    version = version,
    install_requires = [
        'click',
        'numpy',
        'pandas',
        'bokeh',
        'requests',
        'seaborn'
    ],
    include_package_data = True,
    entry_points = {
        'console_scripts': [
            'calculate-plot = plotburden.calculate_plot:cli',
        ],
    },
# Metadata
    author = "Arthur Gilly",
    author_email = "arthur.gilly@helmholtz-muenchen.de",
    maintainer = 'Young-Chan Park',
    maintainer_email = 'young-chan.park@helmholtz-muenchen.de',
    packages = find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
)