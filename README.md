바이오인포매틱스 NGS 입니다.

이 프로젝트는 Linux 또는 Windows WSL2 + Ubuntu 에서 실행하도록 권장드립니다.
이 프로젝트는 또한 python 으로 제작 되었습니다.

# 필요한 패키지
pandas
numpy
biopython
pysam
matplotlib
scipy
scikit-learn
statsmodels
tqdm
PyYAML
openpyxl
jupyterlab

# 외부 NGS 프로그램
- fastp
- BWA 또는 BWA-MEM2
- samtools
- bcftools

## 설치 복사/붙여넣기 용

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
