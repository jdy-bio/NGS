from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
VARIANT_DIR = DATA_DIR / "variants"
REVIEWED_DIR = VARIANT_DIR / "reviewed"
REPORT_DIR = VARIANT_DIR / "final_reports"

REPORT_DIR.mkdir(parents=True, exist_ok=True)

def classify_final_result(row):
    clinvar_significance = str(row.get("ClinVar_Significance", "")).strip()
    local_significance = str(row.get("ClinicalSignificance", "")).strip()

    significance = (
        clinvar_significance
        if clinvar_significance and clinvar_significance.lower() != "nan"
        else local_significance
    )

    normalized = significance.lower()

    if "pathogenic" in normalized and "conflicting" not in normalized:
        if "likely" in normalized:
            return "Likely_pathogenic_candidate"
        return "Pathogenic_candidate"

    if ("uncertain" in normalized or "vus" in normalized or "conflicting" in normalized):
        return "Review_required"

    if "benign" in normalized:
        return "Benign_candidate"

    return "Unclassified"

reviewed_files = sorted(REVIEWED_DIR.glob("*.ncbi_reviewed.csv"))

if not reviewed_files:
    print("NCBI 검토 결과를 찾지 못했습니다.")
    raise SystemExit(1)

success_count = 0
fail_count = 0
all_sample_summaries = []

for reviewed_path in reviewed_files:
    sample = reviewed_path.name.replace(".ncbi_reviewed.csv", "")

    try:
        variants = pd.read_csv(reviewed_path)

        if variants.empty:
            summary = {
                "Sample": sample,
                "Total_variants": 0,
                "Pathogenic_candidates": 0,
                "Likely_pathogenic_candidates": 0,
                "Review_required": 0,
                "Benign_candidates": 0,
                "Final_sample_result": "No_reportable_variant",
            }

            all_sample_summaries.append(summary)
            variants.to_csv(REPORT_DIR / f"{sample}.final_variant_report.csv", index=False,)

            print(f"{sample} 최종 리포트 완료")
            success_count += 1
            continue

        variants["FinalInterpretation"] = variants.apply(classify_final_result, axis=1,)

        interpretation_priority = {
            "Pathogenic_candidate": 1,
            "Likely_pathogenic_candidate": 2,
            "Review_required": 3,
            "Unclassified": 4,
            "Benign_candidate": 5,
        }

        variants["Priority"] = (variants["FinalInterpretation"].map(interpretation_priority).fillna(99))

        sort_columns = ["Priority"]
        ascending_values = [True]

        if "AF" in variants.columns:
            sort_columns.append("AF")
            ascending_values.append(False)

        if "DP" in variants.columns:
            sort_columns.append("DP")
            ascending_values.append(False)

        variants = variants.sort_values(sort_columns,ascending=ascending_values,)

        report_columns = [
            column
            for column in [
                "Variant",
                "Gene",
                "Effect",
                "ProteinChange",
                "QUAL",
                "DP",
                "AF",
                "VARIANT_TYPE",
                "ClinicalSignificance",
                "ClinVar_Significance",
                "ClinVar_ReviewStatus",
                "ClinVar_Conditions",
                "ClinVar_Accession",
                "ClinVar_URL",
                "FinalInterpretation",
            ]
            if column in variants.columns
        ]

        final_report = variants[report_columns].copy()
        csv_path = REPORT_DIR / f"{sample}.final_variant_report.csv"
        final_report.to_csv(csv_path, index=False)

        pathogenic_count = int((variants["FinalInterpretation"] == "Pathogenic_candidate").sum())
        likely_pathogenic_count = int((variants["FinalInterpretation"] == "Likely_pathogenic_candidate").sum())
        review_count = int((variants["FinalInterpretation"] == "Review_required").sum())
        benign_count = int((variants["FinalInterpretation"] == "Benign_candidate").sum())

        if pathogenic_count > 0:
            sample_result = "Pathogenic_variant_candidate_detected"
        elif likely_pathogenic_count > 0:
            sample_result = "Likely_pathogenic_variant_candidate_detected"
        elif review_count > 0:
            sample_result = "Manual_review_required"
        else:
            sample_result = "No_pathogenic_candidate_detected"

        summary = {
            "Sample": sample,
            "Total_variants": int(len(variants)),
            "Pathogenic_candidates": pathogenic_count,
            "Likely_pathogenic_candidates": likely_pathogenic_count,
            "Review_required": review_count,
            "Benign_candidates": benign_count,
            "Final_sample_result": sample_result,
        }

        all_sample_summaries.append(summary)

        txt_lines = [
            f"Sample: {sample}",
            "=" * 60,
            f"Total variants: {len(variants)}",
            f"Pathogenic candidates: {pathogenic_count}",
            f"Likely pathogenic candidates: {likely_pathogenic_count}",
            f"Review required: {review_count}",
            f"Benign candidates: {benign_count}",
            "",
            f"Final sample result: {sample_result}",
            "",
            "Variant details",
            "-" * 60,
        ]

        for _, row in variants.iterrows():
            txt_lines.extend(
                [
                    f"Variant: {row.get('Variant', '')}",
                    f"Gene: {row.get('Gene', '')}",
                    f"Protein change: {row.get('ProteinChange', '')}",
                    f"Depth: {row.get('DP', '')}",
                    f"VAF: {row.get('AF', '')}",
                    f"Local annotation: {row.get('ClinicalSignificance', '')}",
                    f"ClinVar significance: {row.get('ClinVar_Significance', '')}",
                    f"ClinVar review status: {row.get('ClinVar_ReviewStatus', '')}",
                    f"Final interpretation: {row.get('FinalInterpretation', '')}",
                    f"NCBI URL: {row.get('ClinVar_URL', '')}",
                    "",
                ]
            )

        txt_path = REPORT_DIR / f"{sample}.final_variant_report.txt"
        txt_path.write_text("\n".join(txt_lines), encoding="utf-8")
        success_count += 1

    except Exception as error:
        print(f"{sample} 최종 리포트 생성 실패: {error}")
        fail_count += 1

summary_df = pd.DataFrame(all_sample_summaries)
summary_path = REPORT_DIR / "all_samples_variant_summary.csv"
summary_df.to_csv(summary_path, index=False)

if fail_count == 0:
    print(f"전체 최종 리포트 완료: {success_count}개 샘플")
else:
    print(f"전체 최종 리포트 종료: 성공 {success_count}개, 실패 {fail_count}개")
