from pathlib import Path
import subprocess
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PIPELINE_SCRIPTS = [
    "1-R. FASTQ_QC.py",
    "2-R. RNA_alignment_count.py",
    "3-R. RNA_BAM_QC.py",
    "4-R. Count_matrix.py",
    "5-R. DEG_analysis.py",
    "6-R. GO_KEGG_GSEA.py",
]

missing = [name for name in PIPELINE_SCRIPTS if not (SCRIPT_DIR / name).exists()]

if missing:
    print("실행 실패: 다음 파일이 없습니다.")
    for name in missing:
        print(f"- {name}")
    raise SystemExit(1)

for step, name in enumerate(PIPELINE_SCRIPTS, start=1):
    print()
    print("=" * 70)
    print(f"{step}단계 실행: {name}")
    print("=" * 70)

    result = subprocess.run(
        [sys.executable, str(SCRIPT_DIR / name)],
        cwd=SCRIPT_DIR
    )

    if result.returncode != 0:
        print(f"파이프라인 중단: {name} 실행 실패")
        raise SystemExit(result.returncode)

print()
print("Bulk RNA-seq 전체 분석 완료")
