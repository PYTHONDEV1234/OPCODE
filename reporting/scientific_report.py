import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.platypus import Table, TableStyle


def generate_scientific_report(adata, output_dir):

    os.makedirs(output_dir, exist_ok=True)

    scores = adata.obs["V5_OPC_score"].values

    # ============================
    # BASIC STATS
    # ============================

    summary_stats = {
        "Total Cells": len(scores),
        "Mean Score": float(np.mean(scores)),
        "Std Dev": float(np.std(scores)),
        "Max Score": float(np.max(scores)),
        "Min Score": float(np.min(scores)),
        "Non-zero Cells": int(np.sum(scores != 0))
    }

    summary_df = pd.DataFrame.from_dict(summary_stats, orient="index", columns=["Value"])
    summary_csv_path = os.path.join(output_dir, "summary_statistics.csv")
    summary_df.to_csv(summary_csv_path)

    # ============================
    # CLASS MEANS
    # ============================

    class_df = None
    if "class" in adata.obs.columns:
        class_means = (
            adata.obs.groupby("class")["V5_OPC_score"]
            .mean()
            .sort_values(ascending=False)
        )
        class_df = class_means.reset_index()
        class_df.columns = ["Class", "Mean_V5_Score"]
        class_csv_path = os.path.join(output_dir, "class_means.csv")
        class_df.to_csv(class_csv_path, index=False)

    # ============================
    # HISTOGRAM FIGURE
    # ============================

    fig_path = os.path.join(output_dir, "score_distribution.png")
    plt.figure(figsize=(6,4))
    plt.hist(scores, bins=50)
    plt.title("V5 OPC Score Distribution")
    plt.xlabel("V5 OPC Score")
    plt.ylabel("Cell Count")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()

    # ============================
    # PDF REPORT
    # ============================

    pdf_path = os.path.join(output_dir, "V5_OPC_Scientific_Report.pdf")
    doc = SimpleDocTemplate(pdf_path)
    elements = []
    styles = getSampleStyleSheet()

    elements.append(Paragraph("V5 OPC Detection Scientific Report", styles["Title"]))
    elements.append(Spacer(1, 12))

    for key, value in summary_stats.items():
        elements.append(Paragraph(f"{key}: {value}", styles["Normal"]))
        elements.append(Spacer(1, 6))

    elements.append(Spacer(1, 12))

    elements.append(Image(fig_path, width=400, height=300))

    doc.build(elements)

    return {
        "summary_csv": summary_csv_path,
        "class_csv": class_csv_path if class_df is not None else None,
        "figure": fig_path,
        "pdf": pdf_path
    }