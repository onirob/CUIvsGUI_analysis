# Analysis and Reproducibility

This repository contains the analysis code and results for Figliè et al. (2026), *Neither Replacement nor Panacea: Comparing LLM-Based Conversational and Graphical Decision Support in Industrial Tasks*.

## Reporting authority

The article is the official source for the study design, primary analyses, and reported findings. The accompanying Supplementary Material is the source for the task decision reference solutions (S1), robustness analyses (S2), industry-role confounding checks (S3), and copy–paste use checks (S4). If a notebook label, exploratory output, or generated file differs from those reporting sources, the article and Supplementary Material take precedence in their respective scopes. The public repository follows the article's terminology: CUI denotes the chatbot condition and GUI denotes the dashboard condition.

The final sample comprised 134 participants in a 2 × 3 mixed design. Interface was manipulated between participants, while every participant completed low-, medium-, and high-complexity tasks in that fixed order.

The notebooks retain their saved outputs. No analysis was rerun and no numerical result or plot was regenerated while preparing this release. Outputs containing direct identifiers, participant-level rows, timestamps, platform data, or participant-authored conversation text were replaced in place by explicit privacy notices.

## Repository structure

- `code/python/`: preserved analysis notebooks with saved outputs and article-aligned section labels.
- `code/r/`: R scripts used for cumulative-link mixed models and related analyses.
- `results/primary/`: available aggregate outputs for the article's primary analyses; some reported results exist only in the notebooks.
- `results/sensitivity/`: industry-role-adjusted confounding checks corresponding to Supplementary Material S3 (Table S5).
- `results/robustness/`: alternative specifications and secondary checks corresponding to Supplementary Material S2 (Tables S1–S4), including the H5 interface-by-complexity tests.
- `figures/`: selected publication-safe figures available from the working project.
- `docs/`: concise methods, private-input schema, data-availability statement, and release documentation.

## Analysis map

| Hypotheses | Outcome | Primary model |
|---|---|---|
| H1A–H4A | Mental workload | Probit cumulative-link mixed models with participant random intercepts |
| H1B–H4B | Decision accuracy | Fractional-logit generalized estimating equations clustered by participant |
| H1C–H4C | Completion time | Gamma generalized estimating equations with a log link, clustered by participant |
| H5 | Intended reliance (`r1`–`r3`) | Item-wise probit cumulative-link mixed models; interface effects adjusted for task complexity |

H5 is evaluated primarily through the item-wise interface effects. The interface × task-complexity models are secondary robustness analyses from Supplementary Material S2.4 (Table S4), not the primary H5 test. The manuscript-aligned findings are summarized in [`RESULTS.md`](RESULTS.md), and the operational analysis definitions are in [`docs/METHODS.md`](docs/METHODS.md).

## Reproduction boundary

Participant-level source data are not included. Authorized researchers using an approved private copy can inspect the expected fields and derivations in [`docs/INPUT_SCHEMA.md`](docs/INPUT_SCHEMA.md). The notebooks document the original ordered workflow and retain their saved output streams, but they are not a turnkey executable package without the private inputs. See [`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md) for the access boundary.

## Software

The Python environment recorded during the analysis is in `requirements.txt`; required R packages are listed in `code/r/packages.txt`. These files document the environment and were not used to rerun the analyses for this release.

## Licensing

Code is licensed under the MIT License. Documentation, figures, and aggregate results are licensed under Creative Commons Attribution 4.0 International as described in `LICENSE-MATERIALS.md`.
