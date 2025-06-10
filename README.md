# Assessing the impact of the withdrawal of Zinc Oxide from weaner piglets' diet on antimicrobial resistance in _Escherichia coli_


<a href="https://doi.org/10.5281/zenodo.15625229"><img src="https://zenodo.org/badge/998958166.svg" alt="DOI"></a>


## Abstract
Antimicrobial resistance (AMR) is a global concern for both human and animal health. Antimicrobial consumption is particularly high in intensive livestock systems, especially in the pig sector [1]. In the UK, efforts to reduce antimicrobial use include stewardship practices, assurance schemes, and regulatory measures.
AMR arises not only from the use of antimicrobials but also from exposure to heavy metals such as zinc and copper, as well as biocides, which can co-select for resistance through shared cellular mechanisms [2]. Resistance may also be acquired via direct genetic linkage between resistance genes carried on mobile genetic elements.

This study was conducted as part of a broader, multidisciplinary project funded by the BBSRC, aiming to understand the impact of the UK ban on pharmaceutical zinc oxide (ZnO) supplementation in piglet feed on AMR levels and post-weaning diarrhoea (PWD) outcomes [3,4].

The following scripts have been specifically developed to investigate the interplay between zinc oxide use in weaner piglets’ diets and the development of antimicrobial resistance in _Escherichia coli_. These scripts operate using Antimicrobial Susceptibility Testing (AST) data obtained through disc diffusion assays, to determine the phenotypic resistance profiles of _E. coli_ isolates.


## Using the scripts
All the scripts have been designed to work with a .csv file organised as follows:

For the rows:
- One row = one isolate/bacteria,

And for the columns:

_E. coli_ isolates were obtained from weaner piglets faecal samples, and were distinguised by farm of origin, sampling time point, ZnO exposure, individual sample identifier within the sampling time point, agar type (i.e. isolated from selective or unselective agar), and individual isolate number within the sample (as multiple isolates were collected from each sample).
- Column 1 = Farm = identifier of the farm from which the isolate comes from,
- Column 2 = Visit = sampling time point (V1, V2 or V4 in our case),
- Column 3 = Group = pre-withdrawal or post-withdrawal (respectively, received or did not receive ZnO in our case),
- Column 4 = Sample = identifier of the sample within each sampling time point (e.g. W01, W02, etc. in our case),
- Column 5 = Agar = to distinguish isolates collected on Unselective or Selective agar,
- Column 6 = Pick = individual number of the isolate within each sample (e.g. P1, P2, etc in our case),

- Columns 7 to 20 = AST values (inhibition zone diameter in millimeters), one column per antibiotic,

Note: in the scripts, the antibiotics are in the following order (unless otherwise specified): Tetracycline (TET), Ampicillin (AMP), Amoxicillin_Clavulanate (AMC), Apramycin (APR), Gentamicin (GEN), Streptomycin (STR), Spectinomycin (SPE), Sulfamethoxazole_Trimethoprim (SUT), Florfenicol (FLO), Chloramphenicol (CHL), Ceftiofur (CFT), Cefotaxime (CTX), Enrofloxacin (ENR), Ciprofloxacin (CIP)

- Columns 21 to 23 = results of toxins PCR analysis (Yes if toxin detected by PCR, No otherwise) (e.g. STa, LTa, Stx2e),
- Columns 24 to 29 = results of fimbriae PCR analysis (Yes if fimbriae detected by PCR, No otherwise) (e.g. F5, F41, F4, F6, F18, and F17).


## Scripts (Rstudio)
- <code>1_AST_CombinedHistogramsATB.R</code>: displaying the distribution of the isolates based on their Antimicrobial Susceptibility values (i.e. zone diameter in millimeters),
- <code>2_AST_Pheatmap4</code>: displaying the phenotypic resistance profiles of the isolates using heatmaps,
- <code>3_AST_R%ComparisonPerVisit</code>: comparing the resistance profiles in the pre- vs post-withdrawal groups,
- <code>4_AST_MDRComparisonPerVisit_PatternsAsso</code>: comparing multi-drug resistance prevalence in the pre- vs post-withdrawal groups,
- <code>5_Fimbriae analysis</code>: comparing the fimbriae expressed in the pre- vs post-withdrawal groups,
- <code>6_Additional_picks_heatmaps</code>: displaying the phenotypic resistance profiles of the isolates using heatmaps (i.e. heatmaps for additional isolates not considered in the main analysis).


## References
[1] Mulchandani R, Wang Y, Gilbert M, Van Boeckel TP. Global trends in antimicrobial use in food-producing animals: 2020 to 2030. PLOS Glob Public Health. 2023 Feb 1;3(2):e0001305. doi: 10.1371/journal.pgph.0001305.

[2] Li X, Rensing C, Vestergaard G, Arumugam M, Nesme J, Gupta S, Brejnrod AD, Sørensen SJ. 2022. Metagenomic evidence for co-occurrence of antibiotic, biocide and metal resistance genes in pigs. Environ Int. doi: 10.1016/j.envint.2021.106899.

[3] “In Brief: Studying the Impact of the zinc oxide ban”. Veterinary Record 195 (1), p12, 6/13 July 2024; doi.org/10.1002/vetr.4446.

[4] EMA 2017, EMA-Zinc oxide-Annex II, EMEA: Amsterdam, The Netherlands.
