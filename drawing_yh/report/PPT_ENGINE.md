# Shared report-draft PPT engine

`drawing_yh.report.ppt_engine` is the single renderer for project web-report drafts.
The Markdown report remains the source of truth; project code only supplies paths and optional overrides.

## Required project files

```text
sub_reports/
├── REPORT.md
├── _report_config.py
└── report_figs/
```

Each draft chunk in `REPORT.md` declares exactly one `PPT 页标题` and `PPT 副标题`.
`_report_config.py` supplies `SLIDES` and `load_chunks(fix_img_for_pages=False)`.

## Project entry point

```python
from pathlib import Path
from drawing_yh.report.ppt_engine import (
    PptBuildConfig,
    build_draft_ppt,
    default_template_candidates,
    load_draft_specs,
)
import _report_config as report_config

config = PptBuildConfig(
    output=Path("D:/Projects/My_Project/report/draft.pptx"),
    figure_dir=Path(__file__).parent / "report_figs",
    template_candidates=default_template_candidates(),
    title="Project title",
    subtitle="汇报草稿",
)
build_draft_ppt(load_draft_specs(report_config), config)
```

## Shared behavior

- content title 20 pt, full-slide centered;
- subtitle 16 pt, next row, left aligned;
- Microsoft YaHei explanatory text with fixed readable sizes;
- candidate layouts: `side`, `wide_bottom`, `bottom_columns`, `paired_summary`;
- score uses image readability, text fit, space utilization, semantic mapping, and overflow penalty;
- explicit image width overrides and project-specific dense-slide slugs;
- speaker notes from `主信息` and `讲者备注`;
- optional native PowerPoint sections;
- close only the target presentation, back up unsaved edits, atomically replace output, and reopen it.

## Unified launcher

From the `projects-yh` repository root:

```powershell
py -3.13 tools/build_report_ppt.py bone_marrow
py -3.13 tools/build_report_ppt.py b_cell
py -3.13 tools/build_report_ppt.py thymus
py -3.13 tools/build_report_ppt.py maodie
py -3.13 tools/build_report_ppt.py cellchat
py -3.13 tools/build_report_ppt.py all
```

Aliases `bm`, `bone`, and `bcell` are accepted. Existing project-local commands remain compatible.

## Migrated projects

- `Bone_Marrow_Aging/project_report/md_to_pptx.py`
- `B_Cell_Aging/report/sub_reports/build_ppt_draft.py`
- `Thymus_Aging/mouse_clock_bulk_reanalysis/report/md_to_pptx.py`
- `Bone_Marrow_Aging/analysis/maodie_blood_validation/md_to_pptx.py`
- `Bone_Marrow_Aging/analysis/cellchat_x_BMIF/md_to_pptx.py`
- existing `node build_ppt_draft.mjs` remains as a compatibility launcher.

Reports with `PPT 页标题` / `PPT 副标题` metadata use `load_draft_specs()` directly. Historical reports without those fields use `drawing_yh.report.legacy_ppt.load_legacy_specs()` to preserve project title maps, split multi-image pages, and move excess visible bullets into speaker notes; rendering still goes through the same engine.
