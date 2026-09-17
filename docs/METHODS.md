# Analysis Methods

The article is the official source for the design and analysis definitions. This document is a concise map from the article to the repository.

## Design and measures

- **Sample and design:** 134 participants; interface was between participants (CUI/chatbot vs GUI/dashboard), and task complexity was within participants (low, medium, high) in a fixed low-to-high order.
- **Mental workload:** the Raw NASA-TLX score was the mean of six subscales. `tlx_mean` was rounded to five-point increments and modeled as an ordered outcome.
- **Decision accuracy:** task-specific raw scores were converted to percentage-of-maximum-possible (POMP) scores on [0, 1]; task decision reference solutions are in Supplementary Material S1.
- **Completion time:** elapsed time from the participant's “I read the task” action to answer submission.
- **Data literacy:** the mean of the 12-item BDLI, with seven-point responses, categorized as low [1, 5], medium (5, 6], and high (6, 7].
- **Intended reliance:** three ordered items analyzed separately: sole reliance (`r1`), reliance after colleague verification (`r2`), and reliance after verification with another analytical tool (`r3`).

## Primary models

| Outcome | Model family | Repeated-measures handling |
|---|---|---|
| Mental workload | Probit cumulative-link mixed model | Participant random intercept |
| Accuracy | Fractional-logit GEE | Participant clustering |
| Completion time | Gamma GEE with log link | Participant clustering |
| Reliance | Item-wise probit cumulative-link mixed model, interface effect adjusted for complexity | Participant random intercept |

H2 and H3 use the interface interactions specified in the article; H4 tests the task-complexity effect while controlling for interface. Marginal effects and contrasts for the GEE analyses use participant-cluster bootstrap inference where reported, and follow-up families use Holm adjustment. For H5, the item-wise interface main effects are primary; interface × complexity models are Supplementary Material S2.4 (Table S4) robustness checks.

Industry-role-adjusted models are Supplementary Material S3 (Table S5) confounding checks. Other alternative specifications are Supplementary Material S2 (Tables S1–S4) robustness analyses. The supplementary copy–paste use checks are in S4 (Tables S6–S7). The complete-triple decision-accuracy robustness sample is 131 participants (Table S2); the notebook's 134-to-131 diagnostic is an optional audit, not an input to any fitted model.

No model was fitted and no result was recomputed while preparing the public repository.
