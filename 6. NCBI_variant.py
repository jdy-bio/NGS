from pathlib import Path
import time
import pandas as pd
from Bio import Entrez

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
VARIANT_DIR = DATA_DIR / "variants"
ANNOTATED_DIR = VARIANT_DIR / "annotated"
REVIEWED_DIR = VARIANT_DIR / "reviewed"

REVIEWED_DIR.mkdir(parents=True, exist_ok=True)

# 이메일 주소
Entrez.email = "eodud0215@gmail.com"

NCBI_API_KEY = "" # NCBI API key가 있다면 입력
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

MAX_RESULTS = 5
REQUEST_INTERVAL = 0.12 if NCBI_API_KEY else 0.40


def get_nested_value(data, keys, default=""):
    current = data
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def build_clinvar_query(row):
    rs_id = str(row.get("dbSNP", "")).strip()
    if rs_id and rs_id.lower() not in {"nan", "not_annotated"}:
        return f"{rs_id}[All Fields]"

    gene = str(row.get("Gene", "")).strip()
    protein_change = str(row.get("ProteinChange", "")).strip()

    if (
        gene
        and gene.lower() != "nan"
        and gene != "Not_annotated"
        and protein_change
        and protein_change.lower() != "nan"
        and protein_change != "Not_annotated"
    ):
        return f"{gene}[Gene Name] AND {protein_change}[All Fields]"

    variant = str(row.get("Variant", "")).strip()
    if variant:
        return f"{variant}[All Fields]"

    return ""


def search_clinvar(query):
    empty_result = {
        "NCBI_Query": query,
        "ClinVar_Found": False,
        "ClinVar_ID": "",
        "ClinVar_Title": "",
        "ClinVar_Accession": "",
        "ClinVar_Significance": "",
        "ClinVar_ReviewStatus": "",
        "ClinVar_Conditions": "",
        "ClinVar_URL": "",
    }

    if not query:
        return empty_result

    with Entrez.esearch(db="clinvar", term=query, retmax=MAX_RESULTS) as handle:
        search_record = Entrez.read(handle)

    id_list = search_record.get("IdList", [])
    if not id_list:
        return empty_result

    clinvar_id = id_list[0]

    with Entrez.esummary(db="clinvar", id=clinvar_id, retmode="xml") as handle:
        summary_record = Entrez.read(handle)

    summaries = summary_record.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    if not summaries:
        return empty_result

    summary = summaries[0]

    title = str(summary.get("title", summary.get("Title", "")))
    accession = str(summary.get("accession", summary.get("Accession", "")))

    significance = (
        get_nested_value(summary, ["germline_classification", "description"])
        or get_nested_value(summary, ["clinical_significance", "description"])
        or str(summary.get("clinical_significance", ""))
    )

    review_status = (
        get_nested_value(summary, ["germline_classification", "review_status"])
        or get_nested_value(summary, ["clinical_significance", "review_status"])
    )

    conditions = []
    trait_set = summary.get("trait_set", [])
    if isinstance(trait_set, list):
        for trait in trait_set:
            if isinstance(trait, dict):
                name = trait.get("trait_name") or trait.get("name") or trait.get("TraitName")
                if name:
                    conditions.append(str(name))

    return {
        "NCBI_Query": query,
        "ClinVar_Found": True,
        "ClinVar_ID": str(clinvar_id),
        "ClinVar_Title": title,
        "ClinVar_Accession": accession,
        "ClinVar_Significance": str(significance),
        "ClinVar_ReviewStatus": str(review_status),
        "ClinVar_Conditions": "; ".join(conditions),
        "ClinVar_URL": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/",
    }

annotated_files = sorted(ANNOTATED_DIR.glob("*.annotated_variants.csv"))

if not annotated_files:
    print("annotation 결과 파일을 찾지 못했습니다.")
    raise SystemExit(1)

success_count = 0
fail_count = 0

for annotated_path in annotated_files:
    sample = annotated_path.name.replace(".annotated_variants.csv", "")

    try:
        variants = pd.read_csv(annotated_path)

        if variants.empty:
            output_path = REVIEWED_DIR / f"{sample}.ncbi_reviewed.csv"
            variants.to_csv(output_path, index=False)
            print(f"{sample} NCBI 검토 완료: annotation 변이가 없습니다.")
            success_count += 1
            continue

        ncbi_results = []

        for _, row in variants.iterrows():
            query = build_clinvar_query(row)

            try:
                result = search_clinvar(query)
            except Exception as error:
                result = {
                    "NCBI_Query": query,
                    "ClinVar_Found": False,
                    "ClinVar_ID": "",
                    "ClinVar_Title": "",
                    "ClinVar_Accession": "",
                    "ClinVar_Significance": "",
                    "ClinVar_ReviewStatus": "",
                    "ClinVar_Conditions": "",
                    "ClinVar_URL": "",
                    "NCBI_Error": str(error),
                }

            ncbi_results.append(result)
            time.sleep(REQUEST_INTERVAL)

        reviewed = pd.concat(
            [
                variants.reset_index(drop=True),
                pd.DataFrame(ncbi_results).reset_index(drop=True),
            ],
            axis=1,
        )

        output_path = REVIEWED_DIR / f"{sample}.ncbi_reviewed.csv"
        reviewed.to_csv(output_path, index=False)

        found_count = int(reviewed["ClinVar_Found"].fillna(False).sum())
        print(f"{sample} NCBI 검토 및 저장 완료: ClinVar 검색 {found_count}개")
        success_count += 1

    except Exception as error:
        print(f"{sample} NCBI 검토 실패: {error}")
        fail_count += 1

if fail_count == 0:
    print(f"전체 NCBI 검토 완료: {success_count}개 샘플")
else:
    print(f"전체 NCBI 검토 종료: 성공 {success_count}개, 실패 {fail_count}개")