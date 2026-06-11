from pathlib import Path
import subprocess
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PIPELINE_SCRIPTS = [
    "1. FASTQ QC.py",
    "2. Variant_pipeline.py",
    "3. Bam QC.py",
    "4. Variant_filter.py",
    "5. Variant_annotation.py",
    "6. NCBI_variant.py",
    "7. Final_variant.py",
]

missing_scripts = [script_name for script_name in PIPELINE_SCRIPTS if not (SCRIPT_DIR / script_name).exists()]

if missing_scripts:
    print("다음 코드 파일이 없습니다.")

    for script_name in missing_scripts:
        print(f"- {script_name}")

    raise SystemExit(1)

success_count = 0

for step_number, script_name in enumerate(PIPELINE_SCRIPTS, start=1):
    script_path = SCRIPT_DIR / script_name

    print()
    print("=" * 70)
    print(f"{step_number}단계 실행: {script_name}")
    print("=" * 70)

    result = subprocess.run([sys.executable, str(script_path)], cwd=SCRIPT_DIR)

    if result.returncode != 0:
        print()
        print(f"{script_name} 실행 실패")
        
        raise SystemExit(result.returncode)
    success_count += 1
    print(f"{script_name} 완료")

print()
print("=" * 70)
print(f"DNA 변이 분석 전체 완료: {success_count}개 단계")