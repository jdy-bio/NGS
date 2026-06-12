from pathlib import Path
import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / "Analyzed_RNA"
OUTPUT_DIR = SCRIPT_DIR / "rna_data" / "qc"

INPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 품질 기준 설정
MIN_AVG_Q = 30
MIN_LENGTH = 20

r1_files = sorted(list(INPUT_DIR.glob("*_R1.fastq.gz")) + list(INPUT_DIR.glob("*_R1.fastq")))

rows = []
success_count = 0
fail_count = 0

for r1_path in r1_files:
    if r1_path.name.endswith("_R1.fastq.gz"):
        sample = r1_path.name.replace("_R1.fastq.gz", "")
        r2_path = INPUT_DIR / f"{sample}_R2.fastq.gz"
    else:
        sample = r1_path.name.replace("_R1.fastq", "")
        r2_path = INPUT_DIR / f"{sample}_R2.fastq"

    if not r2_path.exists():
        print(f"{sample} R2 파일이 없습니다.")
        fail_count += 1
        continue

    try:
        for mate, fastq_path in {"R1": r1_path, "R2": r2_path}.items():
            read_count = 0
            pass_count = 0
            lengths = []
            mean_q_values = []
            gc_values = []

            for record in SeqIO.parse(fastq_path, "fastq"):
                seq = str(record.seq).upper()
                q = record.letter_annotations["phred_quality"]

                mean_q = sum(q) / len(q)
                gc = (seq.count("G") + seq.count("C")) / len(seq) * 100

                read_count += 1
                lengths.append(len(seq))
                mean_q_values.append(mean_q)
                gc_values.append(gc)

                if mean_q >= MIN_AVG_Q and len(seq) >= MIN_LENGTH:
                    pass_count += 1

            rows.append({
                "sample": sample,
                "mate": mate,
                "reads": read_count,
                "pass_reads": pass_count,
                "pass_rate": round(pass_count / read_count, 4) if read_count else 0,
                "mean_length": round(sum(lengths) / len(lengths), 2) if lengths else 0,
                "mean_q": round(sum(mean_q_values) / len(mean_q_values), 2) if mean_q_values else 0,
                "mean_gc_percent": round(sum(gc_values) / len(gc_values), 2) if gc_values else 0,
            })

        success_count += 1

    except Exception as error:
        print(f"{sample} QC 실패: {error}")
        fail_count += 1

result_df = pd.DataFrame(rows)
output_path = OUTPUT_DIR / "rna_fastq_qc_summary.csv"

try:
    if result_df.empty:
        print("분석된 RNA FASTQ가 없습니다.")
    else:
        result_df.to_csv(output_path, index=False)
        print(f"FASTQ QC 저장 완료: {output_path}")
except Exception as error:
    print(f"저장 실패: {error}")

print(f"FASTQ QC 종료: 성공 {success_count}개, 실패 {fail_count}개")
