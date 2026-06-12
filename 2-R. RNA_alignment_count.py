from pathlib import Path
import platform
import shutil
import subprocess

if platform.system() != "Linux":
    print("Linux 또는 WSL Ubuntu에서 실행하세요.")
    raise SystemExit(1)

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / "Analyzed_RNA"
RNA_DATA_DIR = SCRIPT_DIR / "rna_data"

TRIM_DIR = RNA_DATA_DIR / "trimmed_fastq"
STAR_INDEX_DIR = RNA_DATA_DIR / "star_index"
BAM_DIR = RNA_DATA_DIR / "bam"
COUNT_DIR = RNA_DATA_DIR / "counts"
LOG_DIR = RNA_DATA_DIR / "logs"

for directory in [
    INPUT_DIR,
    TRIM_DIR,
    STAR_INDEX_DIR,
    BAM_DIR,
    COUNT_DIR,
    LOG_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)

# 실행 설정

THREADS = 4

REFERENCE = INPUT_DIR / "reference.fa"
GTF = INPUT_DIR / "annotation.gtf"

required_tools = [
    "fastp",
    "STAR",
    "samtools",
    "featureCounts",
]

missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]

if missing_tools:
    print("필요한 외부 프로그램이 없습니다.")
    print(", ".join(missing_tools))
    raise SystemExit(1)

if not REFERENCE.exists():
    print("Analyzed_RNA 폴더에 reference.fa가 없습니다.")
    raise SystemExit(1)

if not GTF.exists():
    print("Analyzed_RNA 폴더에 annotation.gtf가 없습니다.")
    raise SystemExit(1)

star_index_marker = STAR_INDEX_DIR / "Genome"

if not star_index_marker.exists():
    try:
        subprocess.run(
            [
                "STAR",
                "--runMode", "genomeGenerate",
                "--runThreadN", str(THREADS),
                "--genomeDir", str(STAR_INDEX_DIR),
                "--genomeFastaFiles", str(REFERENCE),
                "--sjdbGTFfile", str(GTF),
                "--sjdbOverhang", "149",
            ],
            check=True
        )
    except subprocess.CalledProcessError as error:
        print(f"STAR index 생성 실패: {error}")
        raise SystemExit(1)

r1_files = sorted(
    list(INPUT_DIR.glob("*_R1.fastq.gz"))
    + list(INPUT_DIR.glob("*_R1.fastq"))
)

if not r1_files:
    print("R1 FASTQ를 찾지 못했습니다.")
    raise SystemExit(1)

success_count = 0
fail_count = 0
bam_files = []

for r1_path in r1_files:
    gzipped = r1_path.name.endswith("_R1.fastq.gz")

    if gzipped:
        sample = r1_path.name.replace("_R1.fastq.gz", "")
        r2_path = INPUT_DIR / f"{sample}_R2.fastq.gz"
    else:
        sample = r1_path.name.replace("_R1.fastq", "")
        r2_path = INPUT_DIR / f"{sample}_R2.fastq"

    if not r2_path.exists():
        print(f"{sample} 실패: R2 파일이 없습니다.")
        fail_count += 1
        continue

    trim_r1 = TRIM_DIR / f"{sample}_R1.trim.fastq.gz"
    trim_r2 = TRIM_DIR / f"{sample}_R2.trim.fastq.gz"
    star_prefix = BAM_DIR / f"{sample}."
    star_bam = BAM_DIR / f"{sample}.Aligned.sortedByCoord.out.bam"
    final_bam = BAM_DIR / f"{sample}.sorted.bam"

    try:
        subprocess.run(
            [
                "fastp",
                "-i", str(r1_path),
                "-I", str(r2_path),
                "-o", str(trim_r1),
                "-O", str(trim_r2),
                "-w", str(THREADS),
                "-h", str(LOG_DIR / f"{sample}.fastp.html"),
                "-j", str(LOG_DIR / f"{sample}.fastp.json"),
            ],
            check=True
        )

        subprocess.run(
            [
                "STAR",
                "--runThreadN", str(THREADS),
                "--genomeDir", str(STAR_INDEX_DIR),
                "--readFilesIn", str(trim_r1), str(trim_r2),
                "--readFilesCommand", "zcat",
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outFileNamePrefix", str(star_prefix),
            ],
            check=True
        )

        star_bam.replace(final_bam)

        subprocess.run(["samtools", "index", "-@", str(THREADS),str(final_bam),], check=True)

        bam_files.append(final_bam)
        success_count += 1
        print(f"{sample} RNA 정렬 완료")

    except subprocess.CalledProcessError as error:
        print(f"{sample} RNA 정렬 실패: {error}")
        fail_count += 1

if not bam_files:
    print("BAM 파일이 없습니다.")
    raise SystemExit(1)

count_output = COUNT_DIR / "featurecounts_raw.txt"

try:
    subprocess.run(
        [
            "featureCounts",
            "-T", str(THREADS),
            "-p",
            "-a", str(GTF),
            "-o", str(count_output),
            *[str(path) for path in bam_files],
        ],
        check=True
    )

    print(f"유전자별 count 생성 완료: {count_output}")

except subprocess.CalledProcessError as error:
    print(f"featureCounts 실패: {error}")
    raise SystemExit(1)

print(f"RNA 정렬 종료: 성공 {success_count}개, 실패 {fail_count}개")
