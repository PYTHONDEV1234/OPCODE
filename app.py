import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

from scoring_engine.v5_opc_scoring_engine import score_v5_opc
from scoring_engine.metadata_utils import merge_metadata
from analysis.effect_size import cohens_d
from analysis.validation_metrics import (
    compute_validation_metrics,
    compute_validation_integrity_score
)
from reporting.pdf_report import generate_pdf_report


# =====================================================
# APP CONFIG
# =====================================================

st.set_page_config(layout="wide")
st.title("OPC Detection Tool — V5 Production Build")
st.markdown("Upload an `.h5ad` dataset to compute V5 OPC scores and generate validation metrics.")


# =====================================================
# FILE UPLOAD
# =====================================================

uploaded_file = st.file_uploader("Upload .h5ad file", type=["h5ad"])
metadata_file = st.file_uploader("Optional: Upload Metadata CSV", type=["csv"])


# =====================================================
# MAIN PIPELINE
# =====================================================

if uploaded_file:

    st.subheader("Loading Dataset")
    adata = sc.read_h5ad(uploaded_file)
    st.write(adata)

    # -------------------------------------------------
    # SCORING
    # -------------------------------------------------

    st.subheader("Running V5 Scoring Engine")

    try:
        adata = score_v5_opc(adata)
    except Exception as e:
        st.error(f"Scoring failed: {e}")
        st.stop()

    # -------------------------------------------------
    # METADATA MERGE (OPTIONAL)
    # -------------------------------------------------

    if metadata_file:
        st.subheader("Merging Metadata")
        try:
            adata = merge_metadata(adata, metadata_file)
        except Exception as e:
            st.warning(f"Metadata merge failed: {e}")

    # =====================================================
    # VALIDATION DASHBOARD
    # =====================================================

    st.header("Validation Dashboard")

    scores = adata.obs["V5_OPC_score"].values

    col1, col2, col3 = st.columns(3)
    col1.metric("Total Cells", len(scores))
    col2.metric("Mean Score", round(np.mean(scores), 3))
    col3.metric("Max Score", round(np.max(scores), 3))

    # -------------------------------------------------
    # SCORE DISTRIBUTION
    # -------------------------------------------------

    st.subheader("V5 Score Distribution")

    fig, ax = plt.subplots()
    ax.hist(scores, bins=50)
    ax.set_xlabel("V5 OPC Score")
    ax.set_ylabel("Cell Count")
    ax.set_title("V5 OPC Score Distribution")
    st.pyplot(fig)

    # -------------------------------------------------
    # CLASS VALIDATION
    # -------------------------------------------------

    validation_metrics = None
    vis_score = None

    if "class" in adata.obs.columns:

        st.subheader("Class-Level Validation")

        class_means = (
            adata.obs.groupby("class")["V5_OPC_score"]
            .mean()
            .sort_values(ascending=False)
        )

        st.dataframe(class_means)

        # Compute automated validation metrics
        validation_metrics = compute_validation_metrics(adata)

        vis_score = compute_validation_integrity_score(validation_metrics)

        st.subheader("Validation Metrics")

        metrics_df = pd.DataFrame([validation_metrics])
        st.dataframe(metrics_df)

        st.subheader("Validation Integrity Score (VIS)")
        st.metric("VIS (0–1 scale)", round(vis_score, 3))

        # Bar chart of validation metrics
        fig2, ax2 = plt.subplots()
        ax2.bar(validation_metrics.keys(), validation_metrics.values())
        ax2.set_title("Validation Metric Breakdown")
        ax2.set_ylabel("Metric Value")
        ax2.set_xticklabels(validation_metrics.keys(), rotation=45)
        st.pyplot(fig2)

    # -------------------------------------------------
    # CLUSTER VALIDATION
    # -------------------------------------------------

    if "clusters" in adata.obs.columns:

        st.subheader("Cluster Ranking")

        cluster_means = (
            adata.obs.groupby("clusters")["V5_OPC_score"]
            .mean()
            .sort_values(ascending=False)
        )

        st.dataframe(cluster_means)

        fig3, ax3 = plt.subplots()
        cluster_means.plot(kind="bar", ax=ax3)
        ax3.set_ylabel("Mean V5 OPC Score")
        ax3.set_title("Cluster Mean Scores")
        st.pyplot(fig3)

    if "cluster_alias" in adata.obs.columns:

        st.subheader("Cluster Alias Ranking")

        cluster_means_alias = (
            adata.obs.groupby("cluster_alias")["V5_OPC_score"]
            .mean()
            .sort_values(ascending=False)
        )

        st.dataframe(cluster_means_alias.head(50))

    # -------------------------------------------------
    # DOWNLOAD CSV
    # -------------------------------------------------

    st.subheader("Download Scored Data")

    csv_data = adata.obs.to_csv(index=False).encode("utf-8")

    st.download_button(
        label="Download Scored Dataset CSV",
        data=csv_data,
        file_name="v5_scored_output.csv",
        mime="text/csv",
    )

    # -------------------------------------------------
    # GENERATE SCIENTIFIC REPORT
    # -------------------------------------------------

    st.subheader("Generate Scientific PDF Report")

    if st.button("Generate Report"):

        try:
            os.makedirs("outputs", exist_ok=True)

            report_path = generate_pdf_report(
                adata=adata,
                output_path="outputs/V5_OPC_Report.pdf",
                project_name="OPC Detection Tool",
                author="Ansel",
                year="2026",
                validation_metrics=validation_metrics,
                vis_score=vis_score
            )

            with open(report_path, "rb") as f:
                st.download_button(
                    label="Download Scientific Report PDF",
                    data=f,
                    file_name="V5_OPC_Report.pdf",
                    mime="application/pdf",
                )

            st.success("Report generated successfully.")

        except Exception as e:
            st.error(f"Report generation failed: {e}")

    st.success("V5 Scoring & Validation Complete.")