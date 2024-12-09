## Prerequisites
1. **System-Level Dependencies**:
  Not required, but depending on your system, might be needed:
    - `pkg-config`: May be required for building some Python libraries.
    - `libmysqlclient-dev`: May be required for `mysqlclient`.
    - `build-essential`: Provides `gcc` for compiling Python extensions.
    - Install these on Ubuntu/Debian:
      ```bash
      sudo apt update
      sudo apt install pkg-config libmysqlclient-dev build-essential
      ```

2. **db.ini File**:
    - Make a db.ini file and paste this (replacing `<user>` and `<password>` with your actual credentials):
      ```plaintext
      [postgresql]
      database=postgres
      host=d3b-warehouse-aurora-prd.d3b.io
      user=<user>
      password=<password>
      ```

3. **Seven Bridges Credentials**:
    - Set up credentials:
      ```bash
      mkdir -p ~/.sevenbridges
      vim ~/.sevenbridges/credentials
      ```
    - Paste this into the credentials file (replacing `<token>` with your actual token):
      ```plaintext
      [default]
      api_endpoint   = https://cavatica-api.sbgenomics.com/v2
      auth_token     = <token>
      advance_access = false
      ```

4. **Seven Bridges Tools**:
    - Install Seven Bridges CLI tools (not needed for ETL, but in case you want it installed):
      ```bash
      bash -c 'curl https://igor.sbgenomics.com/downloads/sb/install.sh -sSf | sudo -H sh'
      pip3 install pipx
      pipx ensurepath
      source ~/.bashrc
      pipx install sbpack
      ```