# Private Input Schema

The analysis code expects authorized, non-public source tables. This schema documents the computational interface and article-defined derivations without publishing participant records or example values.

## Source fields

| Input category | Principal private fields used |
|---|---|
| Group assignments | `participant_id`, `interface_cond` |
| Task responses | `participant_id`, `task_code`, task-specific answer/score fields, `duration_ms`, status fields, and task-event timestamps where applicable |
| Questionnaire responses | `participant_id`, `questionnaire_name`, `task_code`, `item_key`, `value_numeric`, and `value_text` where applicable |
| Platform demographics | Participant identifier and approved demographic variables used in descriptive or confounding checks |
| Chat records | Participant/conversation linkage and message text used only for the Supplementary Material S4 copy–paste use checks |

## Article-defined analysis fields

| Analysis field | Definition |
|---|---|
| `interface_cond` | `chatbot` = CUI; `dashboard` = GUI |
| `task_code` / complexity | `T1` = low/easy, `T2` = medium/mid, `T3` = high/hard; administered in this fixed order |
| `tlx_mean` | Mean of the six Raw NASA-TLX subscales |
| Ordered workload outcome | `tlx_mean` rounded to five-point increments and treated as ordered for probit CLMMs |
| `accuracy_score` | Task-specific raw score transformed by POMP to [0, 1]; task decision reference solutions are in Supplementary Material S1 |
| Completion time | Elapsed time from “I read the task” to answer submission, represented in milliseconds or converted to seconds for modeling |
| `bdli_mean` | Mean of the 12 seven-point BDLI items |
| `dl_ord` | Low [1, 5], medium (5, 6], high (6, 7] |
| Reliance items | `r1` sole reliance; `r2` reliance after colleague verification; `r3` reliance after verification with another analytical tool |

The original notebooks look for these inputs below an `analysis_data/` directory. That directory and all participant-level records are intentionally absent from the public repository.
