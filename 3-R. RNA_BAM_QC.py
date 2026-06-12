from pathlib import Path
import pandas as pd
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
BAM_DIR = SCRIPT_DIR / "rna_data" / "bam"
QC_DIR = SCRIPT_DIR / "rna_data" / "qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

bam_files = sorted(BAM_DIR.glob("*.sorted.bam"))

if not bam_files:
    print(f"BAM 파일을 찾지 못했습니다: {BAM_DIR}")
    raise SystemExit(1)

rows = []

for bam_path in bam_files:
    sample = bam_path.name.replace(".sorted.bam", "")
    total_reads = 0
    mapped_reads = 0
    uniquely_mapped_reads = 0
    mapq_values = []

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:

            for read in bam.fetch(until_eof=True):
                total_reads += 1
                
                if read.is_unmapped:
                    continue
                mapped_reads += 1
                mapq_values.append(read.mapping_quality)
                
                if read.mapping_quality == 255:
                    uniquely_mapped_reads += 1

        mapping_rate = (mapped_reads / total_reads if total_reads > 0 else 0)
        unique_mapping_rate = (uniquely_mapped_reads / mapped_reads if mapped_reads > 0 else 0)
        mean_mapq = (sum(mapq_values) / len(mapq_values) if mapq_values else 0)

        rows.append({
            "sample": sample,
            "total_reads": total_reads,
            "mapped_reads": mapped_reads,
            "mapping_rate": round(mapping_rate, 4),
            "uniquely_mapped_reads": uniquely_mapped_reads,
            "unique_mapping_rate": round(unique_mapping_rate, 4),
            "mean_mapq": round(mean_mapq, 2),
        })
        print(f"{sample} BAM QC 완료")

    except Exception as error:
        print(f"{sample} BAM QC 실패: {error}")

result_df = pd.DataFrame(rows)
output_path = QC_DIR / "rna_bam_qc_summary.csv"

if result_df.empty:
    print("BAM QC 결과가 없습니다.")

else:
    try:
        result_df.to_csv(output_path, index=False)

        print()
        print(result_df)
        print(f"\nBAM QC 저장 완료: {output_path}")

    except Exception as error:
        print(f"BAM QC 결과 저장 실패: {error}")
