from Bio import Entrez

# Ganti dengan accession number yang ingin dicari
accession_number = "MZ099333"

# Gunakan email untuk Entrez
Entrez.email = "raffaelsiahaan@gmail.com"  # Pastikan mengganti dengan email yang valid

# Mendapatkan data dari Entrez
handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
record = handle.read()

# Cari nama spesies dalam data genetik (biasanya dalam baris 'SOURCE')
import re
match = re.search(r"SOURCE\s+(.*)", record)

if match:
    organism = match.group(1).strip()
    print(f"Nama spesies: {organism}")
else:
    print("Nama spesies tidak ditemukan.")
