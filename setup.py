from setuptools import setup, find_packages

setup(
    name="cbio_etl_runner",
    version="1.0.0",
    author="Your Name",
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
            "cbio_etl_runner=cbio_etl_runner:main",
        ]
    },
    package_data={
        "": [
            "COLLABORATIONS/*", 
            "credentials/*",
            "docs/*",
            "PUBLICATIONS/*",
            "REFS/*",
            "scripts/*", 
            "STUDY_CONFIGS/*",
            "utilities/*",
            "LICENSE",
            "*.md",
            "python_requirements_file.txt"
            ]
    },
    python_requires=">=3.6",
)
