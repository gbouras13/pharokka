#!/usr/bin/env python

"""The setup script."""

from sys import version_info
from setuptools import setup, find_packages
import os


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

version_file = open('VERSION', 'r')
v = version_file.read().strip()

setup(
    author="George Bouras",
    author_email='george.bouras@adelaide.edu.au',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="phage annotation program",
    entry_points={
        'console_scripts': [
            'phrokka=phrokka.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='phrokka',
    name='phrokka',
    packages=find_packages(include=['phrokka', 'phrokka.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/gbouras13/phrokka',
    version=v,
    zip_safe=False,
)
