"""Python install requirements for cbio-etl."""

from setuptools import find_packages, setup

setup(
    name="cbio-etl",
    version="2.4.8",
    author="Jessica Wong & Miguel Brown",
    description="A tool/ETL for converting data from CAVATICA and Data Warehouse to PedcBioportal format",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "psycopg2-binary==2.9.10",
        "boto3==1.37.6",
        "pandas==2.2.3",
        "sevenbridges-python==2.11.2",
        "PyYAML==6.0.2",
        "dsnparse==0.2.1",
        "scipy==1.15.2",
        "importer==1.2",
        "truststore==0.10.4",
        "pybedtools==0.12.0",
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
