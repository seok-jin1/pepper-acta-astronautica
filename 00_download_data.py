import os
import subprocess
import requests
import csv
import io

# Setup
PROJECT_ACCESSION = "PRJNA1145089"
OUTPUT_DIR = "external_data/PRJNA1145089"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"🚀 Starting download for {PROJECT_ACCESSION}...")

# Fetch Run Info
url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={PROJECT_ACCESSION}&result=read_run&fields=run_accession,fastq_ftp,sample_alias,sample_title&format=tsv&download=true"
response = requests.get(url)
reader = csv.DictReader(io.StringIO(response.text), delimiter='	')

# Filter and Download
count = 0
for row in reader:
    alias = row.get('sample_alias', '')
    if '16s' in alias or '16S' in alias:  # Filter for 16S samples
        run_acc = row['run_accession']
        ftp_links = row['fastq_ftp'].split(';')
        
        print(f"\n📦 Processing {alias} ({run_acc})...")
        
        for link in ftp_links:
            filename = os.path.basename(link)
            filepath = os.path.join(OUTPUT_DIR, filename)
            
            if os.path.exists(filepath):
                print(f"  ✅ {filename} exists. Skipping.")
                continue
                
            download_url = f"ftp://{link}"
            # Use wget for robust downloading
            cmd = ["wget", "-q", "--show-progress", "-O", filepath, download_url]
            
            try:
                print(f"  ⬇️  Downloading {filename}...")
                subprocess.run(cmd, check=True)
                count += 1
            except subprocess.CalledProcessError:
                print(f"  ❌ Failed to download {filename}")

print(f"\n✨ Download complete! Total files downloaded: {count}")
