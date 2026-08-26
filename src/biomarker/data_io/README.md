# data_io

**Durum:** 🔴 Planlama aşamasında (Hafta 1-2)

Kullanıcının protein dizisini sisteme sokması ve doğrulanması ile ilgili modüller:
- FASTA dosyası yükleme
- Manuel dizi girişi
- UniProt/NCBI ID ile çevrimiçi veri çekme
- Dizi geçerlilik kontrolü

## Planlanan dosyalar
- `fasta_loader.py` — FASTA parse ve doğrulama
- `fetch_uniprot.py` — UniProt API üzerinden veri çekme (pilot veri setinde kullanılan script buraya taşınacak)
