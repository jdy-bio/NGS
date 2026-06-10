from pathlib import Path
import pandas as pd
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
BAM_DIR = DATA_DIR / "bam"
QC_DIR = DATA_DIR / "qc"

QC_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH = (QC_DIR / "bam_qc_summary.csv")

bam_files = sorted(BAM_DIR.glob("*.bam"))

if not bam_files:
    print(" BAM 파일을 찾지 못했습니다.")
    raise SystemExit(1)

rows = []

for bam_path in bam_files:

    total_reads = 0
    mapped_reads = 0
    unmapped_reads = 0
    duplicate_reads = 0

    mapq_values = []

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):

                total_reads += 1

                # mapping 여부 확인
                if read.is_unmapped:
                    unmapped_reads += 1

                else:
                    mapped_reads += 1
                    mapq_values.append(read.mapping_quality)
                    
                # duplicate 여부 확인
                if read.is_duplicate:
                    duplicate_reads += 1

        # 비율 계산
        mapping_rate = (mapped_reads / total_reads if total_reads > 0 else 0)
        duplicate_rate = (duplicate_reads / total_reads if total_reads > 0 else 0)
        mean_mapq = (sum(mapq_values) / len(mapq_values) if mapq_values else 0)

        sample = (bam_path.name.replace(".sorted.bam", "").replace(".bam", ""))

        rows.append({
            "sample": sample,
            "bam_file": bam_path.name,
            "total_reads": total_reads,
            "mapped_reads": mapped_reads,
            "unmapped_reads": unmapped_reads,
            "mapping_rate": round(mapping_rate, 4),
            "mean_mapq": round(mean_mapq, 2),
            "duplicate_reads": duplicate_reads,
            "duplicate_rate": round(duplicate_rate, 4),
        })

        print(f"{sample} BAM QC 완료")

    except Exception as error:
        print(f"{bam_path.name} BAM QC 실패: "f"{error}")

result_df = pd.DataFrame(rows)

try:

    if result_df.empty:
        print("저장할 결과가 없습니다.")

    else:
        result_df.to_csv(OUTPUT_PATH, index=False)
        print(f"저장 완료: {OUTPUT_PATH}")

except Exception as error:
    print(f"저장 실패: {error}")