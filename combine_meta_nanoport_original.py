# programmino per aggiungere la colonna cluster del file metadati di Michele al file metadati originale, facendo check su meandose
import pandas as pd
import numpy as np

filea = "metadata_paper_release_JI_COMBINED_metadata_microlearner_prostate_base_05_24_MODIFIEDYUN.csv"
fileb = "metadata_cluster_nanoport_michele1june26.csv"
output = "metadata_paper_release_JI_COMBINED_metadata_microlearner_prostate_base_05_24_MODIFIEDYUNWITHCLUSTER.csv"

col_id = None  # usa la prima colonna di entrambi i file
col_dose_a = "mean dose to the rectum (Gy)"
col_dose_b = "Dmean_retto_EQD2Gy_Tcorr"
col_cluster = "cluster"

# Carica i file
a = pd.read_csv(filea)
b = pd.read_csv(fileb)

# Prima colonna = ID
id_a = a.columns[0]
id_b = b.columns[0]

# Controlli colonne
for col in [col_dose_a]:
    if col not in a.columns:
        raise ValueError(f"Colonna mancante in filea: {col}")

for col in [col_dose_b, col_cluster]:
    if col not in b.columns:
        raise ValueError(f"Colonna mancante in fileb: {col}")

# Normalizza gli ID come stringhe
a[id_a] = a[id_a].astype(str).str.strip()
b[id_b] = b[id_b].astype(str).str.strip()

# Controlla che tutti gli ID di b siano in a
ids_a = set(a[id_a])
ids_b = set(b[id_b])

missing = ids_b - ids_a
if missing:
    raise ValueError(f"Questi ID sono presenti in fileb ma non in filea: {sorted(missing)}")

# Controlla valori cluster
def convert_cluster(x):
    if pd.isna(x) or str(x).strip().upper() in ["N/A", "NA", ""]:
        return -2

    try:
        val = int(float(x))
    except ValueError:
        raise ValueError(f"Valore cluster non valido: {x}")

    if val not in [1, 2, 3]:
        raise ValueError(f"Valore cluster diverso da 1, 2 o 3: {x}")

    return val

b[col_cluster] = b[col_cluster].apply(convert_cluster)

# Merge temporaneo per controllare le dosi
check = b[[id_b, col_dose_b]].merge(
    a[[id_a, col_dose_a]],
    left_on=id_b,
    right_on=id_a,
    how="left"
)

# Converte le dosi a numeri
check[col_dose_b] = pd.to_numeric(check[col_dose_b], errors="coerce")
check[col_dose_a] = pd.to_numeric(check[col_dose_a], errors="coerce")

# Controllo corrispondenza dose
tol = 1e-6
diff = np.abs(check[col_dose_b] - check[col_dose_a])

bad = check[diff > tol]

if len(bad) > 0:
    print("ATTENZIONE: dosi non corrispondenti:")
    print(bad[[id_b, col_dose_b, col_dose_a]])
    raise ValueError("Trovate dosi non corrispondenti tra i due file.")

# Crea mappa ID -> cluster
cluster_map = dict(zip(b[id_b], b[col_cluster]))

# Aggiunge cluster ad a
a[col_cluster] = a[id_a].map(cluster_map).fillna(-1).astype(int)

# Salva output
a.to_csv(output, index=False)

print(f"File salvato: {output}")