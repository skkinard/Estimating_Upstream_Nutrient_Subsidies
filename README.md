### Estimating Estuarine Assimilation (EA)

#### Stable Isotope Mixing Models (Bayesian, MCMC):
- Utilized Bayesian mixing models incorporating δ34S and δ13C to estimate estuarine assimilation.
- Excluded δ15N due to inability to differentiate source signatures between estuaries and streams.
- Validated model assumptions visually, ensuring <2% of mixtures exceeded source ranges.
- Employed two-source mixing model with regional averages and standard deviations.
- Excluded periphyton and coarse detritus from end-member determination due to confounding effects.
- Combined filamentous algae and aquatic macrophytes as 'Estuary Sources', while 'Stream Sources' included green leaves from riparian vegetation.
- Executed six mixing model frameworks across different site and taxonomic groupings.
- Utilized 'simmr_mcmc' function with 10,000 iterations to derive individual source proportions.

#### Transient Detection Models (Regression):
- Categorized fish and invertebrate species into euryhaline, diadromous, and stenohaline types.
- Employed Transient Presence Proportion (TPP) to assess distribution across study area.
- Investigated influence of annual rainfall and distance to estuaries on TPP using linear regression analysis.

#### Estuarine Assimilation Model (Multiple Regression):
- Developed linear regression models to predict estuarine assimilation.
- Considered variables including annual rainfall, estuary distance, elevation, and their interactions.
- Withheld data from Nueces River to assess effects of Calallen Dam.
- Compared models using Akaike Information Criterion (AIC).

#### Evaluation of Calallen Dam Impact (A/B Test):
- Compared model estimates with observations above and below Calallen Dam to evaluate nutrient movement.
- Utilized R software with 'stats' and 'MASS' packages for calculations.
