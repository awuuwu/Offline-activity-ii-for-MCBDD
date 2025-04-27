# pip install chembl-webresource-client requests
import time
import json
import requests
from collections import defaultdict
from chembl_webresource_client.new_client import new_client

print("Step 1: Retrieving approved drugs from ChEMBL...")
molecule_client = new_client.molecule
approved_iter = molecule_client.filter(max_phase=4).only(['molecule_chembl_id', 'pref_name', 'first_approval'])

approved_drugs = []
for idx, mol in enumerate(approved_iter, 1):
    year_raw = mol.get('first_approval')
    try:
        year = int(year_raw) if year_raw else None
    except ValueError:
        year = None
    approved_drugs.append({
        'chembl_id': mol['molecule_chembl_id'],
        'name': mol.get('pref_name') or '',
        'approval_year': year
    })
    if idx % 500 == 0:
        print(f"  Retrieved {idx} molecules...")

print(f"Total approved drugs retrieved: {len(approved_drugs)}")

print("Sorting drugs by approval year and name...")
approved_drugs_sorted = sorted(
    approved_drugs,
    key=lambda x: ((x['approval_year'] is None, x['approval_year'] or 0), x['name'])
)
print("Sorting complete.")

# Step 2: Fetch UniProt accessions for drugs approved since 2019
print("Step 2: Retrieving protein targets for drugs approved since 2019...")
mechanism_client = new_client.mechanism
target_client = new_client.target

drug_to_accessions = defaultdict(set)
recent_drugs = [d for d in approved_drugs_sorted if d['approval_year'] and d['approval_year'] >= 2019]
total_recent = len(recent_drugs)
print(f"Found {total_recent} drugs approved since 2019.")

for i, drug in enumerate(recent_drugs, 1):
    chembl_id = drug['chembl_id']
    print(f"  ({i}/{total_recent}) Processing drug {chembl_id} ({drug['name']})...")
    mechs = mechanism_client.filter(molecule_chembl_id=chembl_id)
    for mech in mechs:
        target_chembl_id = mech.get('target_chembl_id')
        if not target_chembl_id:
            continue
        try:
            target = target_client.get(target_chembl_id)
        except Exception:
            continue
        if target.get('target_type') != 'SINGLE PROTEIN':
            continue
        for comp in target.get('target_components', []):
            accession = comp.get('accession')
            if accession:
                drug_to_accessions[chembl_id].add(accession)
    time.sleep(0.5)
print("Protein target retrieval complete.")


print("Step 3: Retrieving UniProt keywords for each protein accession...")
UNIPROT_BASE = 'https://rest.uniprot.org/uniprotkb/'

def fetch_uniprot_keywords(accession):
    url = f'{UNIPROT_BASE}{accession}.json'
    resp = requests.get(url, headers={'Accept': 'application/json'})
    resp.raise_for_status()
    data = resp.json()
    return [kw.get('name') for kw in data.get('keywords', [])]

all_accessions = [(drug, acc) for drug, accs in drug_to_accessions.items() for acc in accs]
total_acc = len(all_accessions)
print(f"Total protein accessions to process: {total_acc}")
drug_protein_keywords = defaultdict(dict)

for idx, (chembl_id, acc) in enumerate(all_accessions, 1):
    print(f"  ({idx}/{total_acc}) Fetching keywords for {acc} (drug {chembl_id})...")
    try:
        keywords = fetch_uniprot_keywords(acc)
    except Exception as e:
        print(f"    Warning: Failed to fetch {acc}: {e}")
        keywords = []
    drug_protein_keywords[chembl_id][acc] = keywords
    time.sleep(0.5)

print("Writing results to chembl_uniprot_keywords.json...")
output = {
    'approved_drugs_sorted': approved_drugs_sorted,
    'drug_to_accessions': {k: list(v) for k, v in drug_to_accessions.items()},
    'drug_protein_keywords': drug_protein_keywords
}
with open('chembl_uniprot_keywords.json', 'w') as out_f:
    json.dump(output, out_f, indent=2)
print('Data retrieval complete.')
