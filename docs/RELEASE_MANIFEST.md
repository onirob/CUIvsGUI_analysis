# Public Release Manifest

## Notebook handling

- `01_descriptive_analyses.ipynb`: 42 cells; 77 saved output objects visible; 12 sensitive output objects replaced in place by notices.
- `02_hypothesis_tests.ipynb`: 90 cells; 134 saved output objects visible; 48 sensitive output objects replaced in place by notices.

The visible outputs were recovered from the original notebooks without executing any cell. Aggregate model, diagnostic, and comparison outputs were restored and translated to match the English cell text. Article-reporting labels were aligned without rerunning code. Output-object redaction is limited to participant-level records, direct identifiers, or conversation-level records; absolute paths in otherwise safe saved text outputs are separately normalized to relative module, notebook, or `analysis_data/` paths without removing warnings or results. Three small-cell descriptive sections embedded within otherwise safe H3 output streams were withheld in place; the surrounding model tests, coefficients, contrasts, diagnostics, and conclusions remain visible.

## Included aggregate artifacts

549 result files passed the release filters. `docs/included_files.csv` records their source, public destination, and SHA-256 hash.

## Excluded categories

- All raw source tables at the root of the original `analysis_data` directory.
- Participant-level and cleaned-analysis datasets.
- Participant-authored conversation output and identifying notebook previews.
- Descriptive/count tables with a cell count below 10.
- Verbose working interpretation logs that duplicated safer structured result files.

The build recorded 58 excluded generated artifacts in the private build log; that path-level log is not included publicly. This total includes two H3 subgroup-count CSVs removed during the final small-cell audit.
