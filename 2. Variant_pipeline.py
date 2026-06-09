from pathlib import Path
import platform
import shutil
import subprocess

if platform.system() != "Linux":
    print("실행 실패: Linux 또는 WSL Ubuntu에서 실행하세요.")
    raise SystemExit(1)

SCRIPT_DIR = Path(__file__).resolve().parent

ANALYZED_DIR = SCRIPT_DIR / "Analyzed"
DATA_DIR = SCRIPT_DIR / "data"
TRIM_DIR = DATA_DIR / "trimmed_fastq"
BAM_DIR = DATA_DIR / "bam"
VARIANT_DIR = DATA_DIR / "variants"
LOG_DIR = DATA_DIR / "logs"

for directory in [
    ANALYZED_DIR,
    DATA_DIR,
    TRIM_DIR,
    BAM_DIR,
    VARIANT_DIR,
    LOG_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)

# 실행 설정

THREADS = 4
REFERENCE = ANALYZED_DIR / "reference.fa"

if shutil.which("bwa-mem2"):
    aligner = "bwa-mem2"

elif shutil.which("bwa"):
    aligner = "bwa"

else:
    aligner = None


required_tools = ["fastp", "samtools", "bcftools",]
missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]

if aligner is None:
    missing_tools.append("bwa-mem2 또는 bwa")

if missing_tools:
    print("실행 실패: 필요한 프로그램이 설치되어 있지 않습니다.")
    print("설치되지 않은 프로그램: " + ", ".join(missing_tools))
    print()
    raise SystemExit(1)

if not REFERENCE.exists():
    print("실행 실패: Analyzed 폴더에 reference.fa가 없습니다.")
    raise SystemExit(1)

if aligner == "bwa-mem2":
    reference_index_exists = any(ANALYZED_DIR.glob("reference.fa.*"))

else:
    bwa_index_files = [
        Path(str(REFERENCE) + ".amb"),
        Path(str(REFERENCE) + ".ann"),
        Path(str(REFERENCE) + ".bwt"),
        Path(str(REFERENCE) + ".pac"),
        Path(str(REFERENCE) + ".sa"),
    ]
    reference_index_exists = all(index_file.exists() for index_file in bwa_index_files)

if not reference_index_exists:
    try:
        subprocess.run([aligner, "index", str(REFERENCE),], check=True)
        print("reference index 생성 완료")

    except subprocess.CalledProcessError as error:
        print(f"reference index 생성 실패: {error}")
        raise SystemExit(1)

r1_files = sorted(list(ANALYZED_DIR.glob("*_R1.fastq.gz")) +list(ANALYZED_DIR.glob("*_R1.fastq")))

if not r1_files:
    print("실행 실패: Analyzed 폴더에서 R1 FASTQ를 찾지 못했습니다.")
    raise SystemExit(1)

success_count = 0
fail_count = 0

for r1_path in r1_files:

    if r1_path.name.endswith("_R1.fastq.gz"):
        sample = r1_path.name.replace("_R1.fastq.gz", "")
        r2_path = (ANALYZED_DIR / f"{sample}_R2.fastq.gz")

    else:
        sample = r1_path.name.replace("_R1.fastq", "")
        r2_path = (ANALYZED_DIR / f"{sample}_R2.fastq")

    if not r2_path.exists():
        print(f"{sample} 분석 실패: R2 파일이 없습니다.")
        fail_count += 1
        continue

    trim_r1 = (TRIM_DIR / f"{sample}_R1.trim.fastq.gz")
    trim_r2 = (TRIM_DIR / f"{sample}_R2.trim.fastq.gz")
    bam_path = (BAM_DIR / f"{sample}.sorted.bam")
    raw_vcf = (VARIANT_DIR / f"{sample}.raw.vcf.gz")
    fastp_html = (LOG_DIR / f"{sample}.fastp.html")
    fastp_json = (LOG_DIR / f"{sample}.fastp.json")
    fastp_log = (LOG_DIR / f"{sample}.fastp.log")
    bwa_log = (LOG_DIR / f"{sample}.bwa.log")
    variant_log = (LOG_DIR / f"{sample}.variant.log")

    try:
        # fastp 전처리
        with open(fastp_log, "w") as log_file:
            subprocess.run(
                [
                    "fastp",

                    "-i",
                    str(r1_path),

                    "-I",
                    str(r2_path),

                    "-o",
                    str(trim_r1),

                    "-O",
                    str(trim_r2),

                    "-w",
                    str(THREADS),

                    "-h",
                    str(fastp_html),

                    "-j",
                    str(fastp_json),
                ],

                stdout=log_file,
                stderr=log_file,
                check=True
            )

        # BWA 정렬 + sorted BAM 생성
        align_command = [aligner, "mem", "-t", str(THREADS), str(REFERENCE), str(trim_r1), str(trim_r2),]

        with open(bwa_log, "w") as bwa_log_file:

            with subprocess.Popen(align_command, stdout=subprocess.PIPE, stderr=bwa_log_file) as bwa_process:
                subprocess.run(["samtools", "sort", "-@", str(THREADS), "-o", str(bam_path),], stdin=bwa_process.stdout, check=True)

                if bwa_process.stdout is not None:
                    bwa_process.stdout.close()

                bwa_return_code = (bwa_process.wait())

                if bwa_return_code != 0:
                    raise subprocess.CalledProcessError(bwa_return_code, align_command)

        # BAM index 생성
        subprocess.run(["samtools", "index", "-@", str(THREADS), str(bam_path), ], check=True)


        # 변이 검출
        with open(variant_log, "w") as variant_log_file:

            with subprocess.Popen(["bcftools", "mpileup", "-f", str(REFERENCE), str(bam_path), "-Ou",], stdout=subprocess.PIPE, stderr=variant_log_file) as mpileup_process:
                subprocess.run(["bcftools", "call", "-m", "-v", "-Oz", "-o", str(raw_vcf),], stdin=mpileup_process.stdout, stderr=variant_log_file, check=True)

                if mpileup_process.stdout is not None:
                    mpileup_process.stdout.close()

                mpileup_return_code = (mpileup_process.wait())

                if mpileup_return_code != 0:
                    raise subprocess.CalledProcessError(mpileup_return_code,["bcftools", "mpileup"])

        # VCF index 생성
        subprocess.run(["bcftools", "index", "-t", str(raw_vcf),], check=True)
        print(f"{sample} 분석 및 저장 완료")
        success_count += 1

    except subprocess.CalledProcessError as error:
        print(f"{sample} 분석 실패: " f"{error}")
        print(f"로그 확인: {LOG_DIR}")
        fail_count += 1

if fail_count == 0:
    print(f"전체 분석 완료: "f"{success_count}개 샘플")

else:
    print(f"전체 분석 종료: " f"성공 {success_count}개, "f"실패 {fail_count}개")