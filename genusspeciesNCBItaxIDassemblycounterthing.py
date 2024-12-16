# -*- coding: utf-8 -*-
"""NCBItaxonomygenomecounterthing.py

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1Rs6JKKTOnhBwXJT6Oug-UjKXp8q0saNF
    
    takes genus/species list as nput (see end of script for example) and returns table of ncbi taxIDs, scientific name on file with NCBI, and a count of genome assemblies.
    
    makes sure you do this first:
    
    !pip install Bio
    !pip install ncbi-datasets-pylib

"""


from Bio import Entrez
import requests
import time

verbose = 1
max_retries = 5
naptime = 2

def get_taxon_id_and_genomes(genus_species_list, max_retries):
    Entrez.email = "jjacobs@atcc.org"  # Always provide your email
    results = []

    for name in genus_species_list:
        genus, species = name.split()[0], " ".join(name.split()[1:])
        taxon_id = None
        scientific_name = "Not found"
        assembly_count = 0
        retries = 0

        # NCBI is finicky - so let's try more than once up to max_retries to get the results, waiting naptime seconds between each attempt
        while retries < max_retries:
            try:
                # Get taxon ID and scientific name using Entrez
                handle = Entrez.esearch(db="taxonomy", term=name)
                record = Entrez.read(handle)
                handle.close()
                if record["IdList"]:
                    taxon_id = record["IdList"][0]
                    handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
                    records = Entrez.read(handle)
                    handle.close()
                    scientific_name = records[0]["ScientificName"]
                    if verbose: print("Found ",taxon_id,"\t",scientific_name)
                else:
                    # If the full name isn't found, search for the genus only
                    handle = Entrez.esearch(db="taxonomy", term=genus)
                    record = Entrez.read(handle)
                    handle.close()
                    if record["IdList"]:
                        taxon_id = record["IdList"][0]
                        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
                        records = Entrez.read(handle)
                        handle.close()
                        scientific_name = records[0]["ScientificName"]
                        if verbose: print("Found Genus",taxon_id,"\t",scientific_name)

                if verbose: print(record)

                if taxon_id:
                    # Use NCBI Datasets API v2.0 to get taxon report
                    url = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/dataset_report"
                    headers = {
                        'accept': 'application/json',
                        'content-type': 'application/json'
                    }
                    data = {
                        "taxons": [taxon_id]
                    }
                    response = requests.post(url, headers=headers, json=data)
                    if response.status_code == 200:
                        taxon_report = response.json()
                        report = taxon_report["reports"][0]["taxonomy"]

                        # Extract count of assemblies
                        for count in report["counts"]:
                            if count["type"] == "COUNT_TYPE_ASSEMBLY":
                                assembly_count = count["count"]
                                break

                results.append(f"{genus}\t{species}\t{name}\t{taxon_id}\t{scientific_name}\t{assembly_count}")
                break  # Exit the retry loop if successful
            except Exception as e:
                retries += 1
                if verbose: print(f"Error fetching data for {name}: {e}. Retrying ({retries}/{max_retries})...")
                time.sleep(naptime)  # Wait for a second before retrying
                if retries == max_retries:
                    results.append(f"{genus}\t{species}\t{name}\t{taxon_id}None\t{scientific_name}\tError")
                    if verbose: print(f"Failed to fetch data for {name} after {max_retries} retries.")

    return results

# Example usage
genus_species_list = ["Anabaena sp.", "Anoxystipes aquaeolei", "Aquamonas haywardensis", "Homo sapiens", "Something unknown"]
taxon_ids = get_taxon_id_and_genomes(genus_species_list, max_retries)
for line in taxon_ids:
    print(line)
