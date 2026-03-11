import matplotlib.pyplot as plt
import seaborn as sns


def plot_distribution(scores):
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(scores, bins=50, kde=False, ax=ax)
    ax.set_title("V5 OPC Score Distribution")
    ax.set_xlabel("V5 OPC Score")
    ax.set_ylabel("Cell Count")
    return fig


def plot_class_comparison(df, score_col="V5_OPC_score"):
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(data=df, x="class", y=score_col, ax=ax)
    ax.set_title("V5 OPC Score by Cell Class")
    return fig