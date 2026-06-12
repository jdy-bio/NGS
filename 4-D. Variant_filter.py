from pathlib import Path
import gzip
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
VARIANT_DIR = DATA_DIR / "variants"
FILTERED_DIR = VARIANT_DIR / "filtered"
FILTERED_DIR.mkdir(parents=True, exist_ok=True)

# 필터 기준
MIN_QUAL = 30
MIN_DP = 20
MIN_AF = 0.05

# INFO 필드 딕셔너리 변환 함수
def parse_info(info_text):
    result = {}
    for item in str(info_text).split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result

# DP 추출 함수
def get_depth(row):
    info = parse_info(row["INFO"])
    if "DP" in info:
        return int(str(info["DP"]).split(",")[0])

    sample_values = str(row["SAMPLE"]).split(":")
    format_keys = str(row["FORMAT"]).split(":")

    if "DP" in format_keys:
        dp_index = format_keys.index("DP")
        if dp_index < len(sample_values):
            return int(sample_values[dp_index])
    return 0


# AF 추출 함수
def get_allele_frequency(row):
    info = parse_info(row["INFO"])

    if "AF" in info:
        return float(str(info["AF"]).split(",")[0])

    format_keys = str(row["FORMAT"]).split(":")
    sample_values = str(row["SAMPLE"]).split(":")

    if ("AD" in format_keys and "DP" in format_keys):

        ad_index = format_keys.index("AD")
        dp_index = format_keys.index("DP")

        if (ad_index < len(sample_values) and dp_index < len(sample_values)):

            ad_values = sample_values[ad_index].split(",")
            depth = int(sample_values[dp_index])

            if (len(ad_values) >= 2 and depth > 0):

                alt_depth = int(ad_values[1])
                return alt_depth / depth
    return 0.0

# 변이 분류 함수

def classify_variant(row):
    reasons = []
    if row["QUAL"] < MIN_QUAL:
        reasons.append("LOW_QUAL")

    if row["DP"] < MIN_DP:
        reasons.append("LOW_DP")

    if row["AF"] < MIN_AF:
        reasons.append("LOW_AF")

    if not reasons:
        return "PASS"

    return ";".join(reasons)

vcf_files = sorted(list(VARIANT_DIR.glob("*.raw.vcf.gz")) + list(VARIANT_DIR.glob("*.raw.vcf")))

if not vcf_files:
    print("raw VCF 파일을 찾지 못했습니다.")
    raise SystemExit(1)

success_count = 0
fail_count = 0


for vcf_path in vcf_files:

    if vcf_path.name.endswith( ".raw.vcf.gz"):
        sample = vcf_path.name.replace(".raw.vcf.gz", "")

    else:
        sample = vcf_path.name.replace(".raw.vcf", "")

    try:
        df = pd.read_csv(
            vcf_path,
            comment="#",
            sep="\t",
            header=None,
            names=[
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                "SAMPLE",
            ],
            compression="infer"
        )

        if df.empty:
            print(f"{sample} 검출된 변이가 없습니다.")
            success_count += 1
            continue

        df["QUAL"] = pd.to_numeric(df["QUAL"], errors="coerce").fillna(0)
        df["DP"] = df.apply(get_depth, axis=1)
        df["AF"] = df.apply(get_allele_frequency, axis=1)
        df["VARIANT_TYPE"] = df.apply(lambda row: ("SNV"if(len(str(row["REF"])) == 1 and len(str(row["ALT"])) == 1) else "INDEL"), axis=1)
        df["QC_RESULT"] = df.apply(classify_variant, axis=1)

        all_result_path = (FILTERED_DIR / f"{sample}.variant_qc.csv")
        pass_result_path = (FILTERED_DIR / f"{sample}.pass_variants.csv")
        
        df.to_csv(all_result_path, index=False)

        pass_df = df[df["QC_RESULT"] == "PASS"]
        pass_df.to_csv(pass_result_path, index=False)
        success_count += 1

    except Exception as error:
        print(f"{sample} 변이 필터링 실패:" f"{error}")
        fail_count += 1

if fail_count == 0:
    print(f"전체 변이 필터링 완료: {success_count}개 샘플")
else:
    print(f"전체 필터링 종료: 성공 {success_count}개, 실패 {fail_count}개")