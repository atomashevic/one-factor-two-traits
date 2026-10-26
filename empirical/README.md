# Published PHQ-ADS values

`kroenke2016/` holds values transcribed from the Supplemental Digital Content of
Kroenke et al. (2016), *Psychosomatic Medicine*, 78(6), 716–727 (PMC4927366). Section 7 and Table 5 use them.

| File | Contents |
|---|---|
| `table_s2_correlations.csv` | Table S2: the 16 × 16 item correlations, all trials pooled |
| `table_s2_asymmetries.csv` | Item pairs whose upper and lower triangles in Table S2 differ, and the value used |
| `table_s3_loadings.csv` | Table S3: standardized one-factor, two-factor and bifactor loadings for each trial |
| `extract_tables.py` | Extracts the three files above from the supplement PDF |
| `table_s3.R` | Population one-factor approximation of each trial from the published two-factor loadings |

**The supplement PDF is not redistributed.** To rerun the extraction, download it from PubMed Central (PMC4927366) and save it as `empirical/data/raw/kroenke2016/kroenke2016_supplement.pdf`.

**Reproducing Table 5.** `checks/phq_table5.R` recomputes every Table 5 cell and each number in Section 7 from `table_s3_loadings.csv`.

**Licence.** These values come from the original article, not from these materials, so the licences of this repository do not cover them. Cite Kroenke et al. (2016) when using them.
