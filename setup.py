#!/usr/bin/env python
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="enzywizard-pocket",
    version="1.0.1",                         # From version.py[reference:2][reference:3]
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
)
