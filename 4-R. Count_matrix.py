from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
COUNT_DIR = SCRIPT_DIR / "rna_data" / "counts"

INPUT_PATH = COUNT_DIR / "featurecounts_raw.txt"
OUTPUT_PATH = COUNT_DIR / "gene_counts.csv"

if not INPUT_PATH.exists():
    print("featurecounts_raw.txt 파일이 없습니다.")
    raise SystemExit(1)

try:
    counts = pd.read_csv(INPUT_PATH, sep="\t", comment="#")
    sample_columns = counts.columns[6:]
    rename_map = {}

    for column in sample_columns:
        sample_name = (Path(column).name.replace(".sorted.bam", ""))
        rename_map[column] = sample_name

    result = counts[["Geneid"] + list(sample_columns)].copy()
    result = result.rename(columns={"Geneid": "Gene", **rename_map, })
    result.to_csv(OUTPUT_PATH, index=False)

    print(f"count matrix 저장 완료: {OUTPUT_PATH}")

except Exception as error:
    print(f"count matrix 저장 실패: {error}")
    raise SystemExit(1)