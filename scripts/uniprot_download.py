import requests
import csv
import time
from datetime import date


#UniProt ACCESSION NUMARALARI
ACCESSIONS = [
    "O08762",   
    "A4D126",
    "P21673",
    "R9RZK8",
    "O43312",
    "P32458",
    "Q96SW2",
    "Q96BZ9",
    "P15514",
    "D2IYS2",
    "P02774",
    "O00300",
    "P11916",
    "Q62735",
    "P18858",
    "Q96G21",
    "Q9SJQ9",
    "Q9NUB1",
    "P0ABD3",
    "Q6C2E3",
]

# BioMarker_ID ataması 
BIOMARKER_IDS = [
    "BM0001",
    "BM0002",
    "BM0003",
    "BM0004",
    "BM0005",
    "BM0006",
    "BM0007",
    "BM0008",
    "BM0009",
    "BM0010",
    "BM0011",
    "BM0012",
    "BM0013",
    "BM0014",
    "BM0015",
    "BM0016",
    "BM0017",
    "BM0018",
    "BM0019",
    "BM0020",
]

RETRIEVED_BY = "Şevvalnur ALTINOK"  



#  UniProt'tan bir proteinin verisini çekiyor

def fetch_uniprot_json(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    response = requests.get(url, timeout=15)
    response.raise_for_status()  
    return response.json()

# Gelen JSON dosyasından gerekenleri çekiyor

def parse_record(accession, data):
    
    protein_name = (
        data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "NA")
    )

    organism = data.get("organism", {}).get("scientificName", "NA")
    tax_id = data.get("organism", {}).get("taxonId", "NA")

    sequence = data.get("sequence", {}).get("value", "")
    length = data.get("sequence", {}).get("length", 0)

    entry_type = data.get("entryType", "")
    reviewed = "Reviewed" if "reviewed" in entry_type.lower() else "Unreviewed"

    return {
        "UniProt_Accession": accession,
        "Protein_Name": protein_name,
        "Organism": organism,
        "Taxonomy_ID": tax_id,
        "Sequence": sequence,
        "Sequence_Length": length,
        "Reviewed_Status": reviewed,
        "Source_URL": f"https://www.uniprot.org/uniprotkb/{accession}",
        "Retrieval_Date": date.today().isoformat(),
    }



# Tüm listeyi indiriyor

def download_batch(accession_list):
    
    results = []
    errors = []

    for acc in accession_list:
        try:
            data = fetch_uniprot_json(acc)
            record = parse_record(acc, data)
            results.append(record)
            print(f"OK: {acc} indirildi.")

        except requests.exceptions.Timeout:
            errors.append([acc, "DOWNLOAD_ERROR", "Connection timeout", date.today().isoformat()])
            print(f"HATA: {acc} zaman aşımına uğradı.")

        except requests.exceptions.HTTPError as e:
            errors.append([acc, "DOWNLOAD_ERROR", f"HTTP error: {e}", date.today().isoformat()])
            print(f"HATA: {acc} için HTTP hatası: {e}")

        except Exception as e:
            errors.append([acc, "DOWNLOAD_ERROR", str(e), date.today().isoformat()])
            print(f"HATA: {acc} için beklenmeyen hata: {e}")

        time.sleep(0.5)   

    return results, errors



# CSV tablosu yazılıyor

def write_csv(results, bm_id_map, output_path="biomarker_proteins_20.csv"):
    if not results:
        print("Yazılacak sonuç yok, CSV oluşturulmadı.")
        return

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["BioMarker_ID"] + list(results[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for r in results:
            row = {"BioMarker_ID": bm_id_map[r["UniProt_Accession"]]}
            row.update(r)
            writer.writerow(row)

    print(f"CSV yazıldı: {output_path}")



# FASTA dosyasını yazan ana fonksıyon

def write_fasta(results, bm_id_map, output_path="biomarker_sequences_20.fasta"):
    with open(output_path, "w", encoding="utf-8") as f:
        for r in results:
            bm_id = bm_id_map[r["UniProt_Accession"]]
            protein_clean = r["Protein_Name"].replace(" ", "_")
            organism_clean = r["Organism"].replace(" ", "_")

            header = f">{bm_id}|{r['UniProt_Accession']}|{protein_clean}|{organism_clean}"
            f.write(header + "\n")
            f.write(r["Sequence"] + "\n")

    print(f"FASTA yazıldı: {output_path}")



# source_manifest.csv dosyasını yazan fonksiyon

def write_manifest(results, bm_id_map, retrieved_by, output_path="source_manifest.csv"):
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "BioMarker_ID", "UniProt_Accession", "Source_Database",
            "Source_URL", "Retrieval_Method", "Retrieval_Date",
            "Retrieved_By", "File_Format", "Notes"
        ])

        for r in results:
            bm_id = bm_id_map[r["UniProt_Accession"]]
            writer.writerow([
                bm_id,
                r["UniProt_Accession"],
                "UniProtKB/Swiss-Prot",
                r["Source_URL"],
                "REST API",
                r["Retrieval_Date"],
                retrieved_by,
                "JSON/FASTA",
                ""
            ])

    print(f"Manifest yazıldı: {output_path}")



#  Hata dosyasını yazan fonksiyon

def write_errors(errors, output_path="download_errors.csv"):
    if not errors:
        print("Hiç hata yok.")
        return

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Accession", "Error_Type", "Message", "Date"])
        writer.writerows(errors)

    print(f"Hata kaydı yazıldı: {output_path}")



# Ana çalıştırma bloğu

if __name__ == "__main__":
    if len(ACCESSIONS) != len(BIOMARKER_IDS):
        raise ValueError("ACCESSIONS ve BIOMARKER_IDS listeleri aynı uzunlukta olmalı!")

    # accession -> BioMarker_ID eşleşme sözlüğü
    bm_id_map = dict(zip(ACCESSIONS, BIOMARKER_IDS))

    print(f"{len(ACCESSIONS)} protein indirilecek...\n")

    results, errors = download_batch(ACCESSIONS)

    write_csv(results, bm_id_map)
    write_fasta(results, bm_id_map)
    write_manifest(results, bm_id_map, RETRIEVED_BY)
    write_errors(errors)

    print(f"\nTamamlandı. Başarılı: {len(results)}, Hatalı: {len(errors)}")
