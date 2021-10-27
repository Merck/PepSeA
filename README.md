# Peptide SAR: MAFFT based API for multiple sequence alignment

API for the alignment of multiple peptide sequences represented in the [HELM notation](https://www.pistoiaalliance.org/helm-notation/). Input sequences can contain up to 256 non-natural amino acids.

## Deployment

1) Download and install MAFFT. Proceed to the MAFFT [page](https://mafft.cbrc.jp/alignment/software/). In the ___Download and Installation___ section choose your system and follow the installation instructions. 

    By default MAFFT is installed into the __/usr/local/bin/__ folder. The same path for MAFFT is specified in the __local.env__ file. Please, edit this file if you installed MAFFT into other folder.

2) Create Python 3.8 virtual environment and activate it:
    ```bash
    python3 -m venv ENV_NAME
    source ENV_NAME/bin/activate
    ```

3) Install dependencies specified in __requirements.txt__:
    ```bash
    python3 -m pip install -r requirements.txt
    ```
   
4) Run the API using uvicorn HTTP server. Execute the following command through the terminal, while you are in the root directory:
    ```bash
    uvicorn alignment.api:api --env-file local.env
    ```
   
5) Access the API through a web-browser. Copy the address specified on the last line of the terminal, after execution of the __uvicorn__ command(by default it is http://127.0.0.1:8000), or you can use any API testing tool (e.g. Postman: https://www.postman.com/api-platform/)


## Testing
Execute the following command to run all the tests:

    python -m pytest -v
    
Execute the following command to run a specific test:

    python -m pytest tests/unit_tests/*script_name* -v

## Links    
For more information about MAFFT follow this link https://mafft.cbrc.jp/alignment/software/
 
For more information on FastAPI deployment you can refer to: https://fastapi.tiangolo.com/deployment/
  
For more information on uvicorn package you can refer to: https://www.uvicorn.org/
