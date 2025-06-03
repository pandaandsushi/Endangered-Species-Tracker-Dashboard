import requests
import json

data = [
    "EX_E", "LR_cd", "LR_lc", "LR_nt"
]

api_key = '7eJk5U1XUUz7ASTVtDVG2U4ymDUwsMkh6soG'

# Iterate through the list of red list categories
for category in data:
    if "_" in category:
        kategori_baru = category.replace("_", "/")
        print(kategori_baru)
        url = f'https://api.iucnredlist.org/api/v4/red_list_categories/{kategori_baru}'
    else:
        url = f'https://api.iucnredlist.org/api/v4/red_list_categories/{category}'

    params = {
        'page': 1,  # Halaman pertama
    }

    headers = {
        'Authorization': f'Bearer {api_key}'
    }

    # Menulis hasil langsung ke file JSON berdasarkan kategori
    filename = f'{category}.json'  # Membuat nama file berdasarkan kategori
    
    with open(filename, 'w') as f:
        first_chunk = True  # Menandai apakah ini chunk pertama untuk menambahkan koma setelahnya
        
        while True:
            response = requests.get(url, headers=headers, params=params)

            if response.status_code == 200:
                data = response.json()
                
                if 'assessments' in data:
                    # Menulis setiap data 'assessments' sebagai objek JSON
                    for assessment in data['assessments']:
                        if not first_chunk:
                            f.write(",\n")  # Menambahkan koma setelah objek sebelumnya
                        else:
                            first_chunk = False  # Untuk chunk pertama, tidak perlu koma
                        
                        # Menulis objek 'assessment' ke dalam file
                        json.dump(assessment, f, ensure_ascii=False, indent=4)

                params['page'] += 1 

                # Jika tidak ada data lebih lanjut, berhenti
                if not data['assessments']:
                    break
            else:
                print(f"Terjadi kesalahan: {response.status_code}")
                break

    print(f"Data untuk kategori {category} telah ditulis ke file '{filename}'.")
