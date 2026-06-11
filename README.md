# NGS DNA 변이 분석 파이프라인

Linux 또는 Windows WSL2 + Ubuntu 환경에서 paired-end FASTQ 파일을 입력으로 받아
품질 확인, 전처리, 정렬, BAM 생성, 변이 검출, 변이 필터링, annotation,
NCBI ClinVar 확인, 최종 리포트 생성을 수행하는 학습·포트폴리오용 프로젝트입니다.

## 전체 분석 흐름

```text
FASTQ
→ FASTQ QC
→ fastp 전처리
→ BWA 정렬
→ sorted BAM 생성
→ BAM QC
→ bcftools 변이 검출
→ raw VCF 생성
→ QUAL·DP·AF 필터링
→ 변이 annotation
→ NCBI ClinVar 조회
→ 최종 변이 리포트

## 필요한 Python 패키지

### 현재 파이프라인 실행에 필요한 패키지

- `pandas`: CSV와 결과 테이블 처리
- `numpy`: 수치 계산
- `biopython`: FASTQ·FASTA 처리 및 NCBI Entrez 조회
- `pysam`: BAM 파일 읽기와 QC
- `matplotlib`: 그래프 생성
- `scipy`: 통계 계산
- `scikit-learn`: PCA와 머신러닝 분석 확장
- `statsmodels`: 다중검정 보정과 통계 분석
- `tqdm`: 진행률 표시
- `PyYAML`: 설정 파일 사용 시 활용
- `openpyxl`: Excel 파일 저장
- `jupyterlab`: 분석 결과 확인과 실습

---

## 필요한 외부 NGS 프로그램

- `fastp`: FASTQ 품질 관리와 trimming
- `BWA` 또는 `BWA-MEM2`: FASTQ read 정렬
- `samtools`: BAM 정렬, index, QC
- `bcftools`: SNP·INDEL 검출과 VCF 생성

---

## 한 번에 설치하기

Ubuntu 또는 WSL Ubuntu 터미널에 아래 내용을 그대로 복사해 실행합니다.

```bash
sudo apt update && \
sudo apt install -y \
    python3 \
    python3-pip \
    python3-venv \
    fastp \
    bwa \
    samtools \
    bcftools && \
python3 -m venv .venv && \
source .venv/bin/activate && \
python -m pip install --upgrade pip && \
python -m pip install \
    pandas \
    numpy \
    biopython \
    pysam \
    matplotlib \
    scipy \
    scikit-learn \
    statsmodels \
    tqdm \
    pyyaml \
    openpyxl \
    jupyterlab

## 입력 파일 준비

`Analyzed` 폴더에 다음 파일을 넣습니다.

Analyzed/
├── reference.fa
├── SAMPLE01_R1.fastq.gz
├── SAMPLE01_R2.fastq.gz
├── SAMPLE02_R1.fastq.gz
└── SAMPLE02_R2.fastq.gz

압축되지 않은 FASTQ도 사용할 수 있습니다.

SAMPLE01_R1.fastq
SAMPLE01_R2.fastq

FASTQ 파일명은 다음 규칙을 따라야 합니다.

샘플명_R1.fastq.gz
샘플명_R2.fastq.gz

또는:

샘플명_R1.fastq
샘플명_R2.fastq

샘플이 많아도 R1 파일을 자동으로 검색하고 대응하는 R2 파일을 연결합니다.

기준서열 파일명은 다음과 같아야 합니다.

reference.fa