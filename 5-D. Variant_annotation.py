from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
VARIANT_DIR = DATA_DIR / "variants"
FILTERED_DIR = VARIANT_DIR / "filtered"
ANNOTATION_DIR = VARIANT_DIR / "annotated"

ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)

ANNOTATION_PATH = (VARIANT_DIR / "annotation_table.csv")

if not ANNOTATION_PATH.exists():

    print("annotation_table.csv가 없습니다.")
    raise SystemExit(1)

pass_files = sorted(FILTERED_DIR.glob("*.pass_variants.csv"))

if not pass_files:
    print("PASS 변이 파일을 찾지 못했습니다.")
    raise SystemExit(1)

annotation = pd.read_csv(ANNOTATION_PATH)

success_count = 0
fail_count = 0

for pass_path in pass_files:
    sample = pass_path.name.replace(".pass_variants.csv", "")
    try:
        variants = pd.read_csv(pass_path)

        if variants.empty:
            print("PASS 변이가 없습니다.")
            success_count += 1
            continue

        merged = variants.merge(
            annotation,
            on=[
                "CHROM",
                "POS",
                "REF",
                "ALT"
            ],
            how="left"
        )
 
        merged["Variant"] = (
            merged["CHROM"].astype(str)
            + ":"
            + merged["POS"].astype(str)
            + " "
            + merged["REF"].astype(str)
            + ">"
            + merged["ALT"].astype(str)
        )

        merged["Gene"] = (merged["Gene"].fillna("Not_annotated"))
        merged["Effect"] = (merged["Effect"].fillna("Not_annotated"))
        merged["ProteinChange"] = (merged["ProteinChange"].fillna("Not_annotated"))
        merged["ClinicalSignificance"] = (merged["ClinicalSignificance"].fillna("Unknown"))

        columns =[
            "Variant",
            "Gene",
            "Effect",
            "ProteinChange",
            "ClinicalSignificance",
            "QUAL",
            "DP",
            "AF",
            "VARIANT_TYPE",
            "SAMPLE"
        ]

        result = (merged[columns].sort_values(["AF", "DP"], ascending=[False, False]))
        output_path = (ANNOTATION_DIR / f"{sample}.annotated_variants.csv")
        result.to_csv(output_path, index=False)
        success_count += 1

    except Exception as error:
        print(f"{sample} annotation 실패: {error}")
        fail_count += 1

if fail_count == 0:
    print(f"전체 annotation 성공: {success_count}개 샘플")

else:
    print(f"전체 annotation 종료: 성공 {success_count}개, 실패 {fail_count}개")