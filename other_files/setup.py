#!/usr/bin/env python

import VERSION
"""The setup script."""

from setuptools import setup, find_packages


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

test_requirements = ['pytest>=3', ]

setup(
    author="George Bouras",
    author_email='george.bouras@adelaide.edu.au',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.9',
    ],
    description="phage annotation program",
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='phrokka',
    name='phrokka',
    packages=find_packages(include=['phrokka', 'phrokka.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/gbouras13/phrokka',
    version=VERSION.get_version(),
    zip_safe=False,
)
