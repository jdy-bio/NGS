from pathlib import Path
import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYZED_DIR = SCRIPT_DIR / "Analyzed"
OUTPUT_DIR = SCRIPT_DIR / "data"

ANALYZED_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 품질 기준 설정
MIN_AVG_Q = 30
MIN_LENGTH = 20

rows = []
r1_files_gz = list(ANALYZED_DIR.glob("*_R1.fastq.gz"))
r1_files_fastq = list(ANALYZED_DIR.glob("*_R1.fastq"))
r1_files = sorted(r1_files_gz + r1_files_fastq)

for r1_path in r1_files:

    if r1_path.name.endswith("_R1.fastq.gz"):
        sample = r1_path.name.replace("_R1.fastq.gz","")
        r2_path = (ANALYZED_DIR / f"{sample}_R2.fastq.gz")

    else:

        sample = r1_path.name.replace("_R1.fastq", "")
        r2_path = (ANALYZED_DIR / f"{sample}_R2.fastq")


    if not r2_path.exists():
        continue

    # R1과 R2 파일을 함께 처리하기 위한 딕셔너리 생성
    paired_files = {"R1": r1_path, "R2": r2_path}

    for mate, fastq_path in paired_files.items():

        read_count = 0
        pass_count = 0

        lengths = []
        average_qualities = []
        gc_values = []

        for record in SeqIO.parse(fastq_path, "fastq"):

            sequence = str(record.seq).upper()
            qualities = (record.letter_annotations["phred_quality"])

            average_q = (sum(qualities) / len(qualities))

            gc_percent = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
            read_count += 1
            lengths.append(len(sequence))

            average_qualities.append(average_q)

            gc_values.append(gc_percent)

            if (average_q >= MIN_AVG_Q and len(sequence) >= MIN_LENGTH):
                pass_count += 1

        rows.append({
            "sample": sample,
            "mate": mate,
            "file_name": fastq_path.name,
            "reads": read_count,
            "pass_reads": pass_count,

            "pass_rate": (round(pass_count / read_count, 4) 
                          if read_count > 0
                          else 0),
            
            "mean_length": (round(sum(lengths) / len(lengths), 2)
                if lengths
                else 0),

            "mean_q": (round(sum(average_qualities) / len(average_qualities), 2)
                if average_qualities
                else 0),

            "mean_gc_percent": (round(sum(gc_values) / len(gc_values), 2)
                if gc_values
                else 0)
        })

result_df = pd.DataFrame(rows)

output_path = (OUTPUT_DIR / "fastq_qc_summary.csv")

try:

    if result_df.empty:
        print("Analyzed 폴더에서 분석 가능한 R1/R2 FASTQ 파일이 없습니다.")

    else:
        result_df.to_csv(output_path, index=False)
        print(f"저장 완료: {output_path}")

except Exception as error:
    print(f"저장 실패: {error}")