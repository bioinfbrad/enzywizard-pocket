#!/usr/bin/env python
from setuptools import setup, find_packages
import os

# Read the version from version.py without importing the package
version_file = os.path.join(os.path.dirname(__file__), 'src', 'enzywizard_pocket', 'version.py')
with open(version_file) as f:
    exec(f.read())  # defines __version__

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="enzywizard-pocket",
    version=__version__,
    author="bioinfbrad",
    description="Detect and characterize binding pockets from a cleaned protein structure",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bioinfbrad/enzywizard-pocket",
    package_dir={"": "src"},                  # Code is in the 'src/' directory
    packages=find_packages(where="src"),
    python_requires=">=3.10",
    install_requires=[
        "bio-pyvol",                         # Key dependency for pocket detection
        "biopython",
        "numpy",
        "scipy",
        "pandas",
        "scikit-learn",
        "trimesh",
        "packaging",
    ],
    entry_points={
        "console_scripts": [
            "enzywizard-pocket = enzywizard_pocket.cli:main",   # Confirmed entry point[reference:4][reference:5]
        ],
    },
    include_package_data=True,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
