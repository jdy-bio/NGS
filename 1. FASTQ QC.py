from pathlib import Path
import pandas as pd
from Bio import SeqIO

# 코드 파일이 있는 위치를 기준으로 폴더 설정
SCRIPT_DIR = Path(__file__).resolve().parent

# 분석할 FASTQ 파일을 넣는 폴더
ANALYZED_DIR = SCRIPT_DIR / "Analyzed"

# 분석 결과를 저장하는 폴더
OUTPUT_DIR = SCRIPT_DIR / "data"


# 폴더가 없으면 자동 생성
ANALYZED_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# 3. 품질 기준 설정
MIN_AVG_Q = 30
MIN_LENGTH = 20

rows = []


# 5. Analyzed 폴더에서 R1 파일 자동 검색

# 압축된 FASTQ 검색
r1_files_gz = list(ANALYZED_DIR.glob("*_R1.fastq.gz"))

# 압축되지 않은 FASTQ 검색
r1_files_fastq = list(ANALYZED_DIR.glob("*_R1.fastq"))

# 두 종류를 합친 후 파일명 순서대로 정렬
r1_files = sorted(r1_files_gz + r1_files_fastq)

# R1 파일을 기준으로 샘플별 반복 분석

for r1_path in r1_files:

    # FASTQ 확장자에 따라 샘플 이름 추출

    if r1_path.name.endswith("_R1.fastq.gz"):

        sample = r1_path.name.replace("_R1.fastq.gz","")

        r2_path = (ANALYZED_DIR / f"{sample}_R2.fastq.gz")

    else:

        sample = r1_path.name.replace("_R1.fastq", "")

        r2_path = (ANALYZED_DIR / f"{sample}_R2.fastq")


    # R2 파일이 없으면 현재 샘플 건너뛰기

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

        # FASTQ 파일을 read 단위로 읽기

        for record in SeqIO.parse(fastq_path, "fastq"):

            sequence = str(record.seq).upper()
            qualities = (record.letter_annotations["phred_quality"])

            # read 평균 품질 계산
            average_q = (sum(qualities) / len(qualities))

            # read GC 비율 계산
            gc_percent = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
            read_count += 1
            lengths.append(len(sequence))

            average_qualities.append(average_q)

            gc_values.append(gc_percent)

            # 품질 기준 통과 여부
            if (average_q >= MIN_AVG_Q and len(sequence) >= MIN_LENGTH):
                pass_count += 1


        # FASTQ 파일의 분석 결과 저장

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

# 결과 저장

output_path = (OUTPUT_DIR / "fastq_qc_summary.csv")

try:

    if result_df.empty:
        print(
            "저장 실패: Analyzed 폴더에서 "
            "분석 가능한 R1/R2 FASTQ 파일을 찾지 못했습니다."
        )

    else:
        result_df.to_csv(
            output_path,
            index=False
        )

        print(
            f"저장 완료: {output_path}"
        )

except Exception as error:

    print(
        f"저장 실패: {error}"
    )