#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

setup(
    name='heat1d',
    version='0.1.6',
    description="Thermal model for planetary science applications",
    long_description=readme + '\n\n' + history,
    author="Paul O. Hayne",
    author_email='paul.hayne@lasp.colorado.edu',
    url='https://github.com/phayne/heat1d',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'heat1d=heat1d.cli:main'
        ]
    },
    package_dir={'heat1d':
                 'heat1d'},
    include_package_data=True,
    install_requires=[
        'Click>=6.0',
        'numpy',
        'matplotlib',
        'planets',
    ],
    license="MIT license",
    zip_safe=False,
    keywords='heat1d',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
    test_suite='tests',
    tests_require=['pytest'],
    setup_requires=['pytest_runner'],
)
