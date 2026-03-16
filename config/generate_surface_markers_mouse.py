# config/generate_surface_markers_mouse.py
"""
Generate config/surface_markers_mouse.csv

Goal:
- Provide a curated, extensible set of *mouse* cell-surface candidates that are plausibly usable
  for FACS/MACS enrichment of OPC / oligodendrocyte-lineage populations.
- Keep the CSV simple and stable: columns = gene, protein, notes

How to run (from project root):
  python config/generate_surface_markers_mouse.py

This script overwrites:
  config/surface_markers_mouse.csv
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import List


@dataclass(frozen=True)
class SurfaceMarker:
    gene: str
    protein: str
    notes: str


def _markers_mouse() -> List[SurfaceMarker]:
    """
    Notes philosophy:
    - "Classic OPC surface marker": commonly used or widely accepted for OPC enrichment.
    - "Oligo-lineage / progenitor surface candidate": plausible but context-dependent; validate.
    - "Broad / not specific": useful as secondary marker only.
    - "Caveat": marker is shared with other populations (pericytes, endothelium, immune, etc).
    """
    return [
        # -------------------------
        # Classic / most-used OPC enrichment anchors
        # -------------------------
        SurfaceMarker(
            gene="Pdgfra",
            protein="PDGFRα (CD140a)",
            notes="Classic OPC surface marker; commonly used in mouse OPC enrichment panels.",
        ),
        SurfaceMarker(
            gene="Cspg4",
            protein="NG2 / CSPG4",
            notes="Classic OPC/pericyte-associated surface marker; not specific—combine with other gates.",
        ),
        SurfaceMarker(
            gene="Ptprz1",
            protein="RPTPζ / PTPRZ1",
            notes="OPC/oligo-lineage surface candidate; often enriched in progenitors (validate per dataset).",
        ),
        SurfaceMarker(
            gene="Itga6",
            protein="Integrin α6 (CD49f)",
            notes="Often used in glial progenitor sorting panels; combine with PDGFRA/other markers.",
        ),
        SurfaceMarker(
            gene="Itgb1",
            protein="Integrin β1 (CD29)",
            notes="Broadly expressed; use as secondary marker with more specific gates.",
        ),
        SurfaceMarker(
            gene="Vcam1",
            protein="VCAM1 (CD106)",
            notes="Reported on subsets of progenitors/reactive populations; use with caution and validation.",
        ),

        # -------------------------
        # High-utility “surface-enriched” candidates often seen in OPC-like clusters
        # (useful for cross-dataset signal; not necessarily OPC-specific)
        # -------------------------
        SurfaceMarker(
            gene="Cd9",
            protein="CD9",
            notes="Tetraspanin; frequently enriched in glial populations; useful as a surface-positive component.",
        ),
        SurfaceMarker(
            gene="Cd81",
            protein="CD81",
            notes="Tetraspanin; frequently enriched in glial populations; useful as a surface-positive component.",
        ),
        SurfaceMarker(
            gene="Cd63",
            protein="CD63",
            notes="Tetraspanin/lysosomal membrane marker; surface availability depends on context; validate.",
        ),

        # -------------------------
        # Oligo-lineage / progenitor surface candidates (validate)
        # -------------------------
        SurfaceMarker(
            gene="Lingo1",
            protein="LINGO1",
            notes="Oligo-lineage surface candidate; implicated in myelination biology; antibody availability varies.",
        ),
        SurfaceMarker(
            gene="Gpr17",
            protein="GPR17",
            notes="Oligo-lineage surface candidate (GPCR); expression is stage-dependent (often pre-OL).",
        ),
        SurfaceMarker(
            gene="Epha4",
            protein="EPHA4",
            notes="Receptor tyrosine kinase; neural/glial surface candidate; validate for OPC specificity.",
        ),
        SurfaceMarker(
            gene="Epha5",
            protein="EPHA5",
            notes="Receptor tyrosine kinase; neural/glial surface candidate; validate.",
        ),
        SurfaceMarker(
            gene="Ntrk2",
            protein="TrkB (NTRK2)",
            notes="Neural receptor; not OPC-specific; include only if it repeatedly ranks as enriched in highconf.",
        ),
        SurfaceMarker(
            gene="Fgfr3",
            protein="FGFR3",
            notes="Glial receptor; can be astrocyte-associated in some contexts; use with caution.",
        ),
        SurfaceMarker(
            gene="Notch1",
            protein="NOTCH1",
            notes="Stem/progenitor signaling receptor; broad; may help as secondary marker in specific datasets.",
        ),
        SurfaceMarker(
            gene="Notch2",
            protein="NOTCH2",
            notes="Stem/progenitor receptor; broad; validate.",
        ),

        # -------------------------
        # Adhesion / ECM receptors (common in glial/progenitor programs; validate)
        # -------------------------
        SurfaceMarker(
            gene="Sdc1",
            protein="Syndecan-1 (SDC1 / CD138)",
            notes="Heparan sulfate proteoglycan; not OPC-specific; may appear in progenitor-like states; validate.",
        ),
        SurfaceMarker(
            gene="Sdc4",
            protein="Syndecan-4 (SDC4)",
            notes="Adhesion/ECM; broad; secondary marker only.",
        ),
        SurfaceMarker(
            gene="Ncam1",
            protein="NCAM1 (CD56)",
            notes="Neural cell adhesion; broad neural marker; may help exclude non-neural contaminants rather than identify OPCs.",
        ),
        SurfaceMarker(
            gene="L1cam",
            protein="L1CAM (CD171)",
            notes="Often neuronal; useful as a NEGATIVE marker if present (exclude neurons).",
        ),
        SurfaceMarker(
            gene="Robo1",
            protein="ROBO1",
            notes="Guidance receptor; not OPC-specific; validate.",
        ),

        # -------------------------
        # Lipid / transport / membrane-associated candidates (validate)
        # -------------------------
        SurfaceMarker(
            gene="Lrp1",
            protein="LRP1 (CD91)",
            notes="Broad endocytic receptor; can be expressed in glia; not specific; secondary marker only.",
        ),
        SurfaceMarker(
            gene="Slc2a1",
            protein="GLUT1 (SLC2A1)",
            notes="Often endothelial; useful as a NEGATIVE marker if present (exclude endothelium).",
        ),
        SurfaceMarker(
            gene="Abca1",
            protein="ABCA1",
            notes="Transporter; broad; surface accessibility varies; validate.",
        ),

        # -------------------------
        # Common NEGATIVE depletion markers (these are not 'OPC positive' markers)
        # Put them here so the pipeline can reference them if present.
        # -------------------------
        SurfaceMarker(
            gene="Ptprc",
            protein="CD45 (PTPRC)",
            notes="IMMUNE depletion marker (negative gate).",
        ),
        SurfaceMarker(
            gene="Pecam1",
            protein="CD31 (PECAM1)",
            notes="ENDOTHELIAL depletion marker (negative gate).",
        ),
        SurfaceMarker(
            gene="Epcam",
            protein="EPCAM",
            notes="EPITHELIAL depletion marker (negative gate).",
        ),
        SurfaceMarker(
            gene="Aif1",
            protein="IBA1 (AIF1)",
            notes="Microglia-associated; not a clean surface target; treat as RNA negative marker / contamination flag.",
        ),

        # -------------------------
        # Mature oligodendrocyte “avoid” markers (usually not surface-sort targets, but useful flags)
        # (We include them as caution entries; your pipeline can later use them for exclusion logic.)
        # -------------------------
        SurfaceMarker(
            gene="Mog",
            protein="MOG",
            notes="Mature oligodendrocyte marker; if present strongly, consider excluding for OPC enrichment aims.",
        ),
        SurfaceMarker(
            gene="Mag",
            protein="MAG",
            notes="Mature oligodendrocyte marker; useful for maturity exclusion checks.",
        ),
    ]


def _validate(markers: List[SurfaceMarker]) -> List[SurfaceMarker]:
    cleaned: List[SurfaceMarker] = []
    seen = set()

    for m in markers:
        gene = (m.gene or "").strip()
        if not gene:
            continue
        # Keep original capitalization (mouse genes are usually TitleCase),
        # but deduplicate case-insensitively.
        key = gene.lower()
        if key in seen:
            continue
        seen.add(key)

        cleaned.append(
            SurfaceMarker(
                gene=gene,
                protein=(m.protein or "").strip(),
                notes=(m.notes or "").strip(),
            )
        )

    # Stable sort: keep intended order but ensure "classic anchors" appear early.
    # (We already ordered them first.)
    return cleaned


def write_surface_markers_mouse_csv(out_path: str):
    markers = _validate(_markers_mouse())

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene", "protein", "notes"])
        for m in markers:
            w.writerow([m.gene, m.protein, m.notes])

    print(f"Wrote {len(markers)} markers -> {out_path}")


def main():
    out_path = os.path.join(os.path.dirname(__file__), "surface_markers_mouse.csv")
    write_surface_markers_mouse_csv(out_path)


if __name__ == "__main__":
    main()