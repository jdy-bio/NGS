from pathlib import Path
import re
import pandas as pd
import gseapy as gp

SCRIPT_DIR = Path(__file__).resolve().parent
DEG_DIR = SCRIPT_DIR / "rna_data" / "deg"
PATHWAY_DIR = SCRIPT_DIR / "rna_data" / "pathway"

PATHWAY_DIR.mkdir(parents=True, exist_ok=True)
ORGANISM = "Human"
GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2021_Human",]

ENRICHR_CUTOFF = 0.05
GSEA_PERMUTATIONS = 100
GSEA_MIN_SIZE = 5
GSEA_MAX_SIZE = 500

def safe_filename(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")

deg_files = sorted(DEG_DIR.glob("*_deg_results.csv"))

if not deg_files:
    print(
        f"비교별 DEG 결과 파일을 찾지 못했습니다: "
        f"{DEG_DIR}"
    )
    raise SystemExit(1)

print("확인된 DEG 비교 파일:")

for path in deg_files:
    print(f"- {path.name}")

print()

summary_rows = []

success_count = 0
fail_count = 0

for deg_path in deg_files:

    comparison = deg_path.name.removesuffix("_deg_results.csv")
    comparison_safe = safe_filename(comparison)
    significant_path = (DEG_DIR / f"{comparison}_significant_deg.csv")
    comparison_output_dir = (PATHWAY_DIR / comparison_safe)
    comparison_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"기능 분석 시작: {comparison}")

    if not significant_path.exists():
        print(
            f"유의 DEG 파일이 없습니다: "
            f"{significant_path.name}"
        )

        summary_rows.append({
            "comparison": comparison,
            "go_kegg_status": "SKIPPED",
            "gsea_status": "SKIPPED",
            "significant_gene_count": 0,
            "ranked_gene_count": 0,
            "status": "FAILED",
            "error": "significant DEG file missing",
        })

        fail_count += 1
        print()
        continue

    try:
        deg = pd.read_csv(deg_path)
        significant = pd.read_csv(significant_path)

    except Exception as error:
        print(f"DEG 파일 읽기 실패: {error}")

        summary_rows.append({
            "comparison": comparison,
            "go_kegg_status": "FAILED",
            "gsea_status": "FAILED",
            "significant_gene_count": 0,
            "ranked_gene_count": 0,
            "status": "FAILED",
            "error": str(error),
        })

        fail_count += 1
        print()
        continue

    required_deg_columns = {"Gene", "log2FoldChange",}

    if not required_deg_columns.issubset(deg.columns):
        missing_columns = (required_deg_columns - set(deg.columns))

        error_message = (
            "DEG 결과에 필요한 컬럼이 없습니다: "
            + ", ".join(sorted(missing_columns))
        )

        print(error_message)

        summary_rows.append({
            "comparison": comparison,
            "go_kegg_status": "FAILED",
            "gsea_status": "FAILED",
            "significant_gene_count": 0,
            "ranked_gene_count": 0,
            "status": "FAILED",
            "error": error_message,
        })

        fail_count += 1
        print()
        continue

    if "Gene" not in significant.columns:
        error_message = ("significant DEG 파일에 Gene 컬럼이 없습니다.")
        print(error_message)
        summary_rows.append({
            "comparison": comparison,
            "go_kegg_status": "FAILED",
            "gsea_status": "FAILED",
            "significant_gene_count": 0,
            "ranked_gene_count": 0,
            "status": "FAILED",
            "error": error_message,
        })

        fail_count += 1
        print()
        continue

    gene_list = (significant["Gene"].dropna().astype(str).str.strip())
    gene_list = gene_list[gene_list != ""]
    gene_list = (gene_list.drop_duplicates().tolist())
    go_kegg_status = "SKIPPED"
    go_kegg_error = ""

    if gene_list:

        try:
            enrichr_output_dir = (comparison_output_dir / "enrichr")
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=GENE_SETS,
                organism=ORGANISM,
                outdir=str(enrichr_output_dir),
                cutoff=ENRICHR_CUTOFF,
            )

            go_kegg_result_path = (comparison_output_dir / f"{comparison_safe}_go_kegg_results.csv")
            enr.results.to_csv(go_kegg_result_path, index=False)
            go_kegg_status = "SUCCESS"

            print(f"GO·KEGG 완료: 유전자 {len(gene_list)}개")

        except Exception as error:
            go_kegg_status = "FAILED"
            go_kegg_error = str(error)

            print(f"GO·KEGG 분석 실패: {error}")

    else:
        print("GO·KEGG 분석에 사용할 유의한 DEG가 없습니다.")

    ranked = (deg[["Gene", "log2FoldChange", ]].dropna().copy())
    ranked["Gene"] = (ranked["Gene"].astype(str).str.strip())
    ranked = ranked[ranked["Gene"] != ""]
    ranked["log2FoldChange"] = pd.to_numeric(ranked["log2FoldChange"], errors="coerce")
    ranked = (ranked.dropna(subset=["log2FoldChange"]).drop_duplicates(subset="Gene", keep="first").sort_values( "log2FoldChange", ascending=False))

    gsea_status = "SKIPPED"
    gsea_error = ""

    if not ranked.empty:

        try:
            gsea_output_dir = (comparison_output_dir / "gsea")

            pre_res = gp.prerank(
                rnk=ranked,
                gene_sets=GENE_SETS,
                outdir=str(gsea_output_dir),
                permutation_num=GSEA_PERMUTATIONS,
                seed=42,
                min_size=GSEA_MIN_SIZE,
                max_size=GSEA_MAX_SIZE,
            )

            gsea_result_path = (comparison_output_dir / f"{comparison_safe}_gsea_results.csv")
            pre_res.res2d.to_csv(gsea_result_path, index=False)
            ranked_output_path = (comparison_output_dir / f"{comparison_safe}_gene_ranking.csv")
            ranked.to_csv(ranked_output_path, index=False)
            gsea_status = "SUCCESS"

            print(f"GSEA 완료: ranking 유전자 {len(ranked)}개")

        except Exception as error:
            gsea_status = "FAILED"
            gsea_error = str(error)

            print(f"GSEA 분석 실패: {error}")

    else:
        print("GSEA에 사용할 ranking 데이터가 없습니다.")

    if (go_kegg_status == "SUCCESS" or gsea_status == "SUCCESS"):
        overall_status = "SUCCESS"
        success_count += 1

    else:
        overall_status = "FAILED"
        fail_count += 1

    error_messages = []

    if go_kegg_error:
        error_messages.append(f"GO_KEGG: {go_kegg_error}")

    if gsea_error:
        error_messages.append(f"GSEA: {gsea_error}")

    summary_rows.append({
        "comparison": comparison,
        "go_kegg_status": go_kegg_status,
        "gsea_status": gsea_status,
        "significant_gene_count": len(gene_list),
        "ranked_gene_count": len(ranked),
        "status": overall_status,
        "error": " | ".join(error_messages),
    })

    print(f"{comparison} 기능 분석 종료")
    print()

summary_df = pd.DataFrame(summary_rows)
summary_path = (PATHWAY_DIR / "pathway_analysis_summary.csv")
summary_df.to_csv(summary_path, index=False)

print("전체 GO·KEGG·GSEA 분석 종료")
print(f"성공 비교: {success_count}개")
print(f"실패 비교: {fail_count}개")
print(f"요약 저장: {summary_path}")

if success_count == 0:
    raise SystemExit(1)