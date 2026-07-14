#!/usr/bin/env python3

import argparse
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Aggiunge la colonna cluster2 del secondo CSV al primo CSV, "
            "controllando che ID e cluster coincidano."
        )
    )

    parser.add_argument("file1", help="Primo file CSV")
    parser.add_argument("file2", help="Secondo file CSV")
    parser.add_argument("output", help="File CSV di uscita")

    args = parser.parse_args()

    # Lettura dei file.
    # dtype=str evita che ID come ML-0001 vengano modificati.
    df1 = pd.read_csv(args.file1, dtype=str)
    df2 = pd.read_csv(args.file2, dtype=str)

    # La prima colonna di ciascun file viene considerata la colonna ID.
    id_col1 = df1.columns[0]
    id_col2 = df2.columns[0]

    # Controllo delle colonne necessarie.
    if "cluster" not in df1.columns:
        sys.exit(
            "Errore: nel primo file non esiste una colonna chiamata 'cluster'."
        )

    required_columns_file2 = {"cluster", "cluster2"}
    missing_columns = required_columns_file2 - set(df2.columns)

    if missing_columns:
        sys.exit(
            "Errore: nel secondo file mancano le colonne: "
            + ", ".join(sorted(missing_columns))
        )

    # Rimuove eventuali spazi dagli ID.
    df1[id_col1] = df1[id_col1].str.strip()
    df2[id_col2] = df2[id_col2].str.strip()

    # Controllo degli ID duplicati.
    duplicated_file1 = df1[id_col1].duplicated(keep=False)
    duplicated_file2 = df2[id_col2].duplicated(keep=False)

    if duplicated_file1.any():
        duplicated_ids = df1.loc[duplicated_file1, id_col1].unique()

        sys.exit(
            "Errore: nel primo file sono presenti ID duplicati:\n"
            + "\n".join(map(str, duplicated_ids))
        )

    if duplicated_file2.any():
        duplicated_ids = df2.loc[duplicated_file2, id_col2].unique()

        sys.exit(
            "Errore: nel secondo file sono presenti ID duplicati:\n"
            + "\n".join(map(str, duplicated_ids))
        )

    # Tiene dal secondo file solo le colonne necessarie.
    # Le colonne vengono rinominate per evitare conflitti durante il merge.
    df2_selected = df2[
        [id_col2, "cluster", "cluster2"]
    ].rename(
        columns={
            id_col2: id_col1,
            "cluster": "cluster_file2",
            "cluster2": "cluster2_file2",
        }
    )

    # Left join:
    # tutte le righe del primo file vengono mantenute.
    merged = df1.merge(
        df2_selected,
        on=id_col1,
        how="left",
        indicator=True,
        validate="one_to_one",
    )

    # True quando l'ID del primo file è presente anche nel secondo file.
    id_trovato = merged["_merge"] == "both"

    # Converte cluster in numerico.
    # In questo modo 1 e 1.0 vengono considerati uguali.
    cluster_file1_numeric = pd.to_numeric(
        merged["cluster"],
        errors="coerce",
    )

    cluster_file2_numeric = pd.to_numeric(
        merged["cluster_file2"],
        errors="coerce",
    )

    # ID presente in entrambi i file e cluster uguale.
    cluster_corretto = (
        id_trovato
        & cluster_file1_numeric.notna()
        & cluster_file2_numeric.notna()
        & (cluster_file1_numeric == cluster_file2_numeric)
    )

    # Costruzione della nuova colonna cluster2:
    #
    # -1: ID presente nel primo file ma non nel secondo
    # -2: ID presente in entrambi, ma cluster diverso
    # valore di cluster2 del secondo file: ID e cluster coincidono
    merged["cluster2_output"] = "-1"

    merged.loc[
        id_trovato & ~cluster_corretto,
        "cluster2_output",
    ] = "-2"

    merged.loc[
        cluster_corretto,
        "cluster2_output",
    ] = merged.loc[
        cluster_corretto,
        "cluster2_file2",
    ]

    # Il file finale contiene tutte le colonne originali del primo file
    # più la nuova colonna cluster2.
    output_df = df1.copy()
    output_df["cluster2"] = merged["cluster2_output"]

    output_df.to_csv(args.output, index=False)

    # Statistiche sulle righe del primo file.
    ids_non_trovati = merged.loc[
        ~id_trovato,
        id_col1,
    ].tolist()

    ids_cluster_diverso = merged.loc[
        id_trovato & ~cluster_corretto,
        id_col1,
    ].tolist()

    ids_aggiunti = merged.loc[
        cluster_corretto,
        id_col1,
    ].tolist()

    n_file1 = len(df1)
    n_file2 = len(df2)
    n_aggiunti = len(ids_aggiunti)
    n_non_trovati = len(ids_non_trovati)
    n_cluster_diverso = len(ids_cluster_diverso)
    n_non_aggiunti = n_non_trovati + n_cluster_diverso

    print(f"Righe totali del primo file: {n_file1}")
    print(f"Righe totali del secondo file: {n_file2}")
    print(f"Righe correttamente associate: {n_aggiunti}")
    print(f"Righe con ID non trovato: {n_non_trovati}")
    print(f"Righe con cluster diverso: {n_cluster_diverso}")
    print(f"Righe totali non associate: {n_non_aggiunti}")

    print("\nID presenti nel primo file ma non nel secondo:")
    if ids_non_trovati:
        for patient_id in ids_non_trovati:
            print(patient_id)
    else:
        print("Nessuno")

    print("\nID presenti in entrambi i file ma con cluster diverso:")
    if ids_cluster_diverso:
        for patient_id in ids_cluster_diverso:
            cluster1 = merged.loc[
                merged[id_col1] == patient_id,
                "cluster",
            ].iloc[0]

            cluster2 = merged.loc[
                merged[id_col1] == patient_id,
                "cluster_file2",
            ].iloc[0]

            print(
                f"{patient_id}: "
                f"cluster primo file = {cluster1}, "
                f"cluster secondo file = {cluster2}"
            )
    else:
        print("Nessuno")

    print(f"\nFile di uscita scritto in: {args.output}")


if __name__ == "__main__":
    main()