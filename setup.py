from setuptools import setup, find_packages
import os

setup(
    name="cbioportal_etl",
    version="2.1.0",
    author="Jessica Wong & Miguel Brown",
    description="A tool/ETL for converting data from CAVATICA and Data Warehouse to PedcBioportal format",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "wheel",
        "psycopg2-binary",
        "boto3",
        "pandas",
        "sevenbridges-python",
        "PyYAML",
        "dsnparse",
        "scipy",
        "bravado",
        "mysqlclient",
        "importer"
    ],
    entry_points={
        "console_scripts": [
            "cbioportal_etl=cbioportal_etl.main:main",
        ]
    },
    package_data={
    "cbioportal_etl": [
        "REFS/*",
        "scripts/*",
        "STUDY_CONFIGS/*",
        ],
    },
    python_requires=">=3.6",
)
