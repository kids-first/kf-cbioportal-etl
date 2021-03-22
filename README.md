# Outline on ETL for converting data from cavatica and data service to pedcbioportal format
In general, we are creating upload packages converting our data and metadata to satisfy the requirements outlined [here](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats)
## Software Prerequisites

+ `python3` v3.5.3+
+ `pip install -r requirements.txt`
+ `bedtools` (https://bedtools.readthedocs.io/en/latest/content/installation.html)


## Processing
` python main.py --kf_ids <kf ids> -c <cavatica_token> -k <kf_token>`

### help
`python main.py -h`
```
usage: main.py [-h] -i KF_IDS [KF_IDS ...] -c CAVATICA_TOKEN -k KF_TOKEN

optional arguments:
  -h, --help            show this help message and exit
  -i KF_IDS [KF_IDS ...], --kf_ids KF_IDS [KF_IDS ...]
                        List of kf study ids
  -c CAVATICA_TOKEN, --cavatica_token CAVATICA_TOKEN
                        Cavatica token
  -k KF_TOKEN, --kf_token KF_TOKEN
                        KF bearer token token
```