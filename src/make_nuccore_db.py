import subprocess
from Bio import Entrez
import pandas as pd
import os
from dotenv import load_dotenv
from tqdm import tqdm
from pathlib import Path
import argparse
import yaml
from email.utils import parsedate_to_datetime
from datetime import date
import urllib.request
import shutil

load_dotenv()

Entrez.email = os.getenv('EMAIL')
Entrez.api_key = os.getenv('NCBI_API_KEY')

def download_dataset(accessions_file: Path, dehydrated_file: Path, db_path: Path, max_workers = 10, api_key = None):
    print("Downloading Dehydrated Dataset...")
    subprocess.run([
        "./datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        accessions_file,
        "--dehydrated",
        "--filename",
        dehydrated_file,  
        "--include",
        "gbff",
        ] + (["--api-key", api_key] if api_key is not None else []),
        check = True
    )
    print()
    
    print("Unzipping Dehydrated Dataset...")
    subprocess.run([
        "unzip",
        dehydrated_file,
        "-d",
        db_path,
        ], 
        check = True
    )
    print()
    
    print("Rehydrating Dataset...")
    subprocess.run([
        "./datasets",
        "rehydrate",
        "--directory",
        db_path,
        "--max-workers",
        str(max_workers)
        ] + (["--api-key", api_key] if api_key is not None else [])
    )
    print()
    



load_dotenv()

REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
NCBI_API_KEY = os.getenv('NCBI_API_KEY')

parser = argparse.ArgumentParser()
parser.add_argument("dbPath")
parser.add_argument("-u", "--updatedb", action="store_true")
parser.add_argument("-m", "--max_workers", help="max workers for rehydration", default = 10)
args = parser.parse_args()

db_path = Path(args.dbPath)

conf_path = db_path / "config.yaml"


with open(conf_path, 'r') as f:
    conf = yaml.safe_load(f)
    f.close()

refseq_path = db_path / conf['refseq_path']
lineage_path = db_path / conf['lineage_path']
log_path = db_path / conf["log_path"]
genomes_dir = db_path / conf['genbank_path']

last_timestamp = date.fromisoformat(conf['last_updated'])

with urllib.request.urlopen(REFSEQ_URL) as r:
    new_timestamp = parsedate_to_datetime(r.headers["Last-Modified"]).date()
    
    
print("Harvesting Accessions...")
lineages = pd.read_csv(lineage_path, sep = '\t')
accessions = list(lineages['accession'])
with open('temp/accessions.txt', 'w') as f:
    f.write('\n'.join(accessions))
    f.close()

download_dataset(
    Path("temp/accessions.txt"), 
    Path("temp/genomes.zip"), 
    genomes_dir, 
    max_workers=int(30),
    api_key=NCBI_API_KEY
)

fasta_dir = genomes_dir / "ncbi_dataset" / "data"
file_count = sum(1 for x in fasta_dir.rglob('*') if x.is_file()) - 2
with open ("log.txt", 'w') as log:
    for file in tqdm(fasta_dir.rglob("*"), total=file_count):
        try:
            if file.is_file():
                dest_path = genomes_dir / file.name
                if file.suffix == ".gbff":
                    dest_path = genomes_dir / (file.parent.name + file.suffix)
                shutil.move(file, dest_path)
        except Exception as e:
            log.write(str(e)+ "\n")
    log.close()
shutil.rmtree(fasta_dir)
    

    
        
        
    
