# analysis

**Durum:** 🔴 Planlama aşamasında (Hafta 3-12)

Protein dizisi üzerinde çalışan tüm analiz algoritmaları burada yer alır.
Her modül birbirinden bağımsız çalışabilir, `data_io` üzerinden aldığı
diziyle beslenir.

## Planlanan dosyalar
- `physicochemical.py` — Moleküler ağırlık, izoelektrik nokta, amino asit bileşimi, instabilite indeksi
- `hydropathy.py` — Kyte-Doolittle hidropati profili (kaydırma penceresi metodu)
- `secondary_structure.py` — Chou-Fasman algoritması ile ikincil yapı tahmini
- `motif_finder.py` — PROSITE benzeri motiflerin regex ile taranması
- `homology_blast.py` — NCBI BLAST API entegrasyonu
