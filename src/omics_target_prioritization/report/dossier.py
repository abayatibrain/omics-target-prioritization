"""Jinja2 HTML one-pager dossier for a prioritized target.

The dossier renders, for a single :class:`~omics_target_prioritization.models.TargetScore`:

- the headline score and calibrated confidence label,
- a per-layer colocalization table (layer, tissue, PP.H4),
- a full provenance table (source, dataset, method, timestamp, parameters),

so a reviewer can see, on one page, *why* the target scored as it did.
"""

from __future__ import annotations

from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from omics_target_prioritization.evidence.confidence import CONFIDENCE_RULE_DOC
from omics_target_prioritization.models import TargetScore

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_TEMPLATE_NAME = "dossier.html.j2"


def _environment() -> Environment:
    """Build the Jinja2 environment bound to the package template directory."""
    return Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html", "xml"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )


def render_dossier(target: TargetScore, *, locus_id: str = "") -> str:
    """Render a target dossier to an HTML string.

    Parameters
    ----------
    target
        The prioritized target to document.
    locus_id
        Optional locus identifier shown in the header.

    Returns
    -------
    str
        The rendered HTML document.
    """
    env = _environment()
    template = env.get_template(_TEMPLATE_NAME)

    coloc_rows = [
        {
            "source": item.source,
            "layer": item.layer or "-",
            "tissue": item.tissue or "-",
            "score": f"{item.score:.3f}",
        }
        for item in target.evidence
        if item.source.endswith("_coloc")
    ]
    provenance_rows = [
        {
            "source": item.source,
            "dataset": item.provenance.dataset,
            "method": item.provenance.method,
            "timestamp": item.provenance.timestamp.isoformat(),
            "score": f"{item.score:.3f}",
        }
        for item in target.evidence
    ]

    return template.render(
        target=target,
        locus_id=locus_id,
        coloc_rows=coloc_rows,
        provenance_rows=provenance_rows,
        breakdown=sorted(target.breakdown.items()),
        confidence_rule=CONFIDENCE_RULE_DOC,
    )


def write_dossier(target: TargetScore, out_path: Path, *, locus_id: str = "") -> Path:
    """Render a dossier and write it to ``out_path``.

    Parameters
    ----------
    target
        The prioritized target to document.
    out_path
        Destination HTML file path. Parent directories are created.
    locus_id
        Optional locus identifier shown in the header.

    Returns
    -------
    pathlib.Path
        The path written (echoes ``out_path``).
    """
    html = render_dossier(target, locus_id=locus_id)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    return out_path
