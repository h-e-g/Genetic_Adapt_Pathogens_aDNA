# Genetic_Adapt_Pathogens_aDNA

Code for reproducing main analyses of Kerner et al., Genetic adaptation to pathogens and increased risk of inflammatory disorders in post-Neolithic Europe (2022).

In this work, we use simulations and an approximate Bayesian computation approach based on time-series of allele frequency trajectories computed from publicly available ancient DNA (1240k capture array) to unravel candidate variants and loci under natural selection over the last 10,000 years in Europe.

For candidate regions, we estimate, at the variant level, selection intensity (s) assuming an additive model with a dominance coefficient of 0.5, and onset of selection (T), i.e., the age at which selection started. 

We discovered 89 loci, defined on the basis of LD, that are candidates for positive selection and 50 candidate missense variants in conserved positions of the genome that are candidates for negative selection. 

We perform enrichment analyses to assess the impact of the possibly positively-selected loci (see scripts) and functional analyses to assess the candidate pathogenic role of the negatively-selected variants (see paper). 

We describe here also a polygenic risk score-based approach to study the evolution of variants associated, in present Europeans by GWAS, to infectious and autoimmune traits. 

Estimations and summary statistics for all tested variants are reported in matrix format and are available for replication purposes. Other datasets used to assess enrichments, such as a curated list of immunity genes, are also provided alongside scripts.
