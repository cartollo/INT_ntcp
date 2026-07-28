#!/bin/sh

N="$1"
K="$2"
OUTPUT_DIR="$3"
ROOT_PREFIX="$4"

if [ -z "$N" ] || [ -z "$K" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$ROOT_PREFIX" ]; then
    echo "Usage: sh bootstrap.sh <N> <K> <output_dir> <root_prefix>"
    exit 1
fi

if [ "$N" -le 0 ] || [ "$K" -le 0 ]; then
    echo "N e K devono essere maggiori di zero."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

worker=1

while [ "$worker" -le "$K" ]; do
    (
        btfile="${OUTPUT_DIR}/bootstrapout_${worker}.txt"
        index="$worker"

        echo "Avvio ciclo $worker, output append: $btfile"

        while [ "$index" -le "$N" ]; do
            rootfile="${OUTPUT_DIR}/${ROOT_PREFIX}_${index}.root"

            echo "Ciclo $worker: bootstrap $index -> $rootfile"

            ./myntcpanalysis \
                -out "$rootfile" \
                -seed "$index" \
                -btfilename "$btfile"

            index=$((index + K))
        done

        echo "Ciclo $worker completato."
    ) &

    worker=$((worker + 1))
done

wait

echo "Tutti i bootstrap completati."