from setuptools import setup, find_packages
import os

setup(
    name="cbio-etl",
    version="2.0.0",
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
            "cbio-etl=cbioportal_etl.cli:main",
        ]
    },
    package_data={
        "cbioportal_etl": [
            "REFS/*",
            "scripts/*",
            "STUDY_CONFIGS/*",
        ],
    },
    python_requires=">=3.10",
)
