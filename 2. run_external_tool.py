from pathlib import Path
import shutil
import subprocess

# 코드가 있는 폴더 기준으로 경로 설정
SCRIPT_DIR = Path(__file__).resolve().parent

# 분석할 FASTQ와 reference.fa를 넣는 폴더
ANALYZED_DIR = SCRIPT_DIR / "Analyzed"

# 결과를 저장하는 폴더
DATA_DIR = SCRIPT_DIR / "data"

TRIM_DIR = DATA_DIR / "trimmed_fastq"
BAM_DIR = DATA_DIR / "bam"
VARIANT_DIR = DATA_DIR / "variants"
LOG_DIR = DATA_DIR / "logs"

# 폴더 자동 생성

ANALYZED_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR.mkdir(parents=True, exist_ok=True)

TRIM_DIR.mkdir(parents=True, exist_ok=True)

BAM_DIR.mkdir(parents=True, exist_ok=True)

VARIANT_DIR.mkdir(parents=True, exist_ok=True)

LOG_DIR.mkdir(parents=True, exist_ok=True)


# 실행 설정

RUN_EXTERNAL_TOOLS = False
THREADS = 4

REFERENCE = ANALYZED_DIR / "reference.fa"


# 외부 프로그램 실행 여부

if not RUN_EXTERNAL_TOOLS:
    print(
        "실행 중지: "
        "RUN_EXTERNAL_TOOLS를 True로 변경하세요."
    )
    raise SystemExit(0)


# =========================================================
# 5. 필요한 외부 프로그램 확인
# =========================================================

if shutil.which("bwa-mem2"):
    aligner = "bwa-mem2"

elif shutil.which("bwa"):
    aligner = "bwa"

else:
    aligner = None


required_tools = [
    "fastp",
    "samtools",
    "bcftools"
]

missing_tools = [
    tool
    for tool in required_tools
    if shutil.which(tool) is None
]

if aligner is None:
    missing_tools.append(
        "bwa-mem2 또는 bwa"
    )


if missing_tools:
    print(
        "실행 실패: 필요한 프로그램이 없습니다."
    )

    print(
        ", ".join(missing_tools)
    )

    raise SystemExit(1)


# =========================================================
# 6. reference 파일 확인
# =========================================================

if not REFERENCE.exists():
    print(
        "실행 실패: "
        "Analyzed 폴더에 reference.fa가 없습니다."
    )

    raise SystemExit(1)


# =========================================================
# 7. R1 FASTQ 파일 자동 검색
# =========================================================

r1_files = sorted(
    list(
        ANALYZED_DIR.glob(
            "*_R1.fastq.gz"
        )
    )
    +
    list(
        ANALYZED_DIR.glob(
            "*_R1.fastq"
        )
    )
)


if not r1_files:
    print(
        "실행 실패: "
        "Analyzed 폴더에서 R1 FASTQ를 찾지 못했습니다."
    )

    raise SystemExit(1)


# =========================================================
# 8. 샘플별 분석
# =========================================================

success_count = 0
fail_count = 0


for r1_path in r1_files:

    # -----------------------------------------------------
    # 샘플 이름과 R2 경로 생성
    # -----------------------------------------------------

    if r1_path.name.endswith(
        "_R1.fastq.gz"
    ):
        sample = r1_path.name.replace(
            "_R1.fastq.gz",
            ""
        )

        r2_path = (
            ANALYZED_DIR
            / f"{sample}_R2.fastq.gz"
        )

    else:
        sample = r1_path.name.replace(
            "_R1.fastq",
            ""
        )

        r2_path = (
            ANALYZED_DIR
            / f"{sample}_R2.fastq"
        )


    # -----------------------------------------------------
    # R2 존재 여부 확인
    # -----------------------------------------------------

    if not r2_path.exists():
        print(
            f"{sample} 분석 실패: "
            "R2 파일이 없습니다."
        )

        fail_count += 1
        continue


    # -----------------------------------------------------
    # 결과 파일 경로 설정
    # -----------------------------------------------------

    trim_r1 = (
        TRIM_DIR
        / f"{sample}_R1.trim.fastq.gz"
    )

    trim_r2 = (
        TRIM_DIR
        / f"{sample}_R2.trim.fastq.gz"
    )

    bam_path = (
        BAM_DIR
        / f"{sample}.sorted.bam"
    )

    raw_vcf = (
        VARIANT_DIR
        / f"{sample}.raw.vcf.gz"
    )

    fastp_html = (
        LOG_DIR
        / f"{sample}.fastp.html"
    )

    fastp_log = (
        LOG_DIR
        / f"{sample}.fastp.log"
    )


    try:

        # =================================================
        # 9. fastp 전처리
        # =================================================

        with open(
            fastp_log,
            "w"
        ) as log_file:

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
                    str(fastp_html)
                ],

                stdout=log_file,
                stderr=log_file,
                check=True
            )


        # =================================================
        # 10. BWA 정렬 + sorted BAM 생성
        # =================================================

        align_command = [
            aligner,
            "mem",

            "-t",
            str(THREADS),

            str(REFERENCE),
            str(trim_r1),
            str(trim_r2)
        ]


        with subprocess.Popen(
            align_command,
            stdout=subprocess.PIPE
        ) as bwa_process:

            subprocess.run(
                [
                    "samtools",
                    "sort",

                    "-@",
                    str(THREADS),

                    "-o",
                    str(bam_path)
                ],

                stdin=bwa_process.stdout,
                check=True
            )


        # =================================================
        # 11. BAM index 생성
        # =================================================

        subprocess.run(
            [
                "samtools",
                "index",
                str(bam_path)
            ],

            check=True
        )


        # =================================================
        # 12. 변이 검출
        # =================================================

        with subprocess.Popen(
            [
                "bcftools",
                "mpileup",

                "-f",
                str(REFERENCE),

                str(bam_path),

                "-Ou"
            ],

            stdout=subprocess.PIPE
        ) as mpileup_process:

            subprocess.run(
                [
                    "bcftools",
                    "call",

                    "-m",
                    "-v",

                    "-Oz",

                    "-o",
                    str(raw_vcf)
                ],

                stdin=mpileup_process.stdout,
                check=True
            )


        # =================================================
        # 13. VCF index 생성
        # =================================================

        subprocess.run(
            [
                "bcftools",
                "index",

                "-t",
                str(raw_vcf)
            ],

            check=True
        )


        print(
            f"{sample} 분석 및 저장 완료"
        )

        success_count += 1


    except subprocess.CalledProcessError as error:

        print(
            f"{sample} 분석 실패: "
            f"{error}"
        )

        fail_count += 1


# =========================================================
# 14. 전체 완료 메시지
# =========================================================

if fail_count == 0:
    print(
        f"전체 분석 완료: "
        f"{success_count}개 샘플"
    )

else:
    print(
        f"전체 분석 종료: "
        f"성공 {success_count}개, "
        f"실패 {fail_count}개"
    )