from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / "Analyzed_RNA"
COUNT_DIR = SCRIPT_DIR / "rna_data" / "counts"
DEG_DIR = SCRIPT_DIR / "rna_data" / "deg"
COUNT_PATH = COUNT_DIR / "gene_counts.csv"
GROUP_PATH = INPUT_DIR / "sample_groups.csv"

FDR_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0

CONTROL_GROUP = None

DEG_DIR.mkdir(parents=True, exist_ok=True)

def safe_filename(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")


def choose_control_group(groups, requested_control=None):
    if requested_control is not None:
        if requested_control not in groups:
            raise ValueError(f"지정한 기준군 '{requested_control}'이 sample_groups.csv에 없습니다.")
        return requested_control

    preferred_names = [
        "Control", "control", "CONTROL",
        "Ctrl", "ctrl", "CTRL",
        "C", "Vehicle", "vehicle", "Untreated", "untreated",
    ]

    for name in preferred_names:
        if name in groups:
            return name

    print(f"'{groups[0]}'을 기준군으로 사용합니다.")
    return groups[0]


def add_direction_column(deg_df):
    result = deg_df.copy()

    result["Direction"] = np.select(
        [
            ((result["padj"] <= FDR_CUTOFF) & (result["log2FoldChange"] >= LOG2FC_CUTOFF)),
            ((result["padj"] <= FDR_CUTOFF) & (result["log2FoldChange"] <= -LOG2FC_CUTOFF)),
        ],
        ["Up", "Down"],
        default="Not_significant",
    )

    return result.sort_values(["padj", "pvalue"], na_position="last",)

def save_volcano_plot(deg_df, comparison_name):
    plot_df = deg_df.copy()
    plot_df["safe_padj"] = (plot_df["padj"].fillna(1).clip(lower=1e-300))

    plt.figure()

    for label in ["Not_significant", "Up", "Down"]:
        part = plot_df[plot_df["Direction"] == label]
        plt.scatter(part["log2FoldChange"], -np.log10(part["safe_padj"]), label=label, alpha=0.7,)

    plt.axvline(LOG2FC_CUTOFF, linestyle="--")
    plt.axvline(-LOG2FC_CUTOFF, linestyle="--")
    plt.axhline(-np.log10(FDR_CUTOFF), linestyle="--")

    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(comparison_name)
    plt.legend()
    plt.tight_layout()
    plt.savefig(DEG_DIR / f"{safe_filename(comparison_name)}_volcano.png", dpi=160,)
    plt.close()

if not COUNT_PATH.exists():
    print(f"gene_counts.csv가 없습니다: {COUNT_PATH}")
    raise SystemExit(1)

if not GROUP_PATH.exists():
    print(f"sample_groups.csv가 없습니다: {GROUP_PATH}")
    raise SystemExit(1)

counts = pd.read_csv(COUNT_PATH)
metadata = pd.read_csv(GROUP_PATH)

if "Gene" not in counts.columns:
    print("gene_counts.csv에는 Gene 열이 필요합니다.")
    raise SystemExit(1)

required_metadata_columns = {"sample", "group"}

if not required_metadata_columns.issubset(metadata.columns):
    print("sample_groups.csv에는 sample, group 열이 필요합니다.")
    raise SystemExit(1)

if metadata["sample"].duplicated().any():
    duplicates = metadata.loc[metadata["sample"].duplicated(), "sample"].tolist()
    print("sample_groups.csv에 중복 샘플이 있습니다:")
    print(", ".join(map(str, duplicates)))
    raise SystemExit(1)

if metadata["group"].isna().any():
    print("sample_groups.csv의 group 열에 빈 값이 있습니다.")
    raise SystemExit(1)

samples = metadata["sample"].astype(str).tolist()
missing_samples = [sample for sample in samples if sample not in counts.columns]

if missing_samples:
    print("count matrix에 다음 샘플이 없습니다.")
    print(", ".join(missing_samples))
    raise SystemExit(1)

count_matrix = counts.set_index("Gene")[samples].T

try:
    count_matrix = count_matrix.astype(int)
except ValueError as error:
    print(f"count matrix를 정수로 변환할 수 없습니다: {error}")
    raise SystemExit(1)

if (count_matrix < 0).any().any():
    print("count matrix에는 음수 값이 들어갈 수 없습니다.")
    raise SystemExit(1)

metadata = metadata.set_index("sample").loc[samples]
metadata["group"] = metadata["group"].astype(str)

group_counts = metadata["group"].value_counts()
print("확인된 그룹 및 샘플 수:")
print(group_counts)
print()

if len(group_counts) < 2:
    print("DEG 분석에는 최소 2개 그룹이 필요합니다.")
    raise SystemExit(1)

small_groups = group_counts[group_counts < 2]
if not small_groups.empty:
    print("경고: 반복 샘플이 1개뿐인 그룹이 있습니다.")
    print(small_groups)
    print("통계 결과의 신뢰도가 낮을 수 있습니다.\n")

dds = DeseqDataSet(counts=count_matrix, metadata=metadata, design_factors="group", refit_cooks=True,)
dds.deseq2()

groups = metadata["group"].dropna().unique().tolist()

try:
    control_group = choose_control_group(groups, CONTROL_GROUP)
except ValueError as error:
    print(error)
    raise SystemExit(1)

comparison_groups = [group for group in groups if group != control_group]

print(f"기준군: {control_group}")
print("비교 대상:", ", ".join(comparison_groups))
print()

summary_rows = []

for treatment_group in comparison_groups:
    comparison_name = f"{treatment_group}_vs_{control_group}"
    comparison_file = safe_filename(comparison_name)

    print("=" * 60)
    print(f"DEG 분석 중: {comparison_name}")

    try:
        stats = DeseqStats(dds, contrast=["group", treatment_group, control_group,], )

        stats.summary()
        deg_df = stats.results_df.reset_index()
        deg_df = deg_df.rename(columns={"index": "Gene"})
        deg_df = add_direction_column(deg_df)

        result_path = DEG_DIR / f"{comparison_file}_deg_results.csv"
        deg_df.to_csv(result_path, index=False)

        significant = deg_df[deg_df["Direction"] != "Not_significant"].copy()

        significant_path = (DEG_DIR / f"{comparison_file}_significant_deg.csv")
        significant.to_csv(significant_path, index=False)

        up_count = int((significant["Direction"] == "Up").sum())
        down_count = int((significant["Direction"] == "Down").sum())

        save_volcano_plot(deg_df, comparison_name)

        summary_rows.append({
            "comparison": comparison_name,
            "control": control_group,
            "treatment": treatment_group,
            "total_genes": len(deg_df),
            "significant_genes": len(significant),
            "up_genes": up_count,
            "down_genes": down_count,
            "status": "SUCCESS",
            "error": "",
        })

        print(
            f"{comparison_name} 완료: "
            f"Up {up_count}개, Down {down_count}개"
        )

    except Exception as error:
        summary_rows.append({
            "comparison": comparison_name,
            "control": control_group,
            "treatment": treatment_group,
            "total_genes": 0,
            "significant_genes": 0,
            "up_genes": 0,
            "down_genes": 0,
            "status": "FAILED",
            "error": str(error),
        })
        print(f"{comparison_name} 분석 실패: {error}")

normalized = dds.layers["normed_counts"]
log_norm = np.log2(normalized + 1)

scaled = StandardScaler().fit_transform(log_norm)
pca = PCA(n_components=2)
coordinates = pca.fit_transform(scaled)

pca_df = pd.DataFrame({
    "Sample": samples,
    "Group": metadata["group"].values,
    "PC1": coordinates[:, 0],
    "PC2": coordinates[:, 1],
})

pca_df.to_csv(DEG_DIR / "pca_coordinates.csv", index=False,)
plt.figure()

for group in pca_df["Group"].unique():
    part = pca_df[pca_df["Group"] == group]
    plt.scatter(part["PC1"], part["PC2"], label=group,)

for _, row in pca_df.iterrows():
    plt.text(row["PC1"], row["PC2"], row["Sample"], fontsize=8,)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
plt.title("RNA-seq PCA")
plt.legend()
plt.tight_layout()
plt.savefig(DEG_DIR / "pca_plot.png", dpi=160)
plt.close()

summary_df = pd.DataFrame(summary_rows)
summary_path = DEG_DIR / "deg_comparison_summary.csv"
summary_df.to_csv(summary_path, index=False)

success_count = int((summary_df["status"] == "SUCCESS").sum())
fail_count = int((summary_df["status"] == "FAILED").sum())

print("전체 DEG 분석 종료")
print(f"성공 {success_count}개 비교, 실패 {fail_count}개 비교")
print(f"요약 저장: {summary_path}")