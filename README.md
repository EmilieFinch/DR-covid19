# DR-covid19

This repo contains code and data to model the COVID-19 pandemic in the Dominican Republic from 2020-2022. This analysis uses `covidm`, a model developed by Davies and colleagues to estimate the impact of NPIs on COVID-19 transmission in the UK. Model fitting requires installation of `Rcpp`, instructions for the installation of covidm can be found in [R/covidm_for_fitting](https://github.com/EmilieFinch/DR-covid19/tree/main/R/covidm_for_fitting)

Model set up and fitting is performed in the following scripts:
- [R/00_load-data.R](https://github.com/EmilieFinch/DR-covid19/blob/main/R/00_load-data.R)
- [R/01_covidm-set-up.R](https://github.com/EmilieFinch/DR-covid19/blob/main/R/01_covidm-set-up.R)
- [R/02_burden-processes.R](https://github.com/EmilieFinch/DR-covid19/blob/main/R/02_burden-processes.R)
- [R/03_fit-covidm.R](https://github.com/EmilieFinch/DR-covid19/blob/main/R/03_fit-covidm.R)

Then, scenario analysis is performed in scripts:
- [R/04_run-vaccine-counterfactuals.R](https://github.com/EmilieFinch/DR-covid19/blob/main/R/04_run-vaccine-counterfactuals.R)
- [R/05_run-counterfactual-trade-off.R](https://github.com/EmilieFinch/DR-covid19/blob/main/R/05_run-counterfactual-trade-off.R)

Figures included in the manuscript can be generated by knitting [reports/manuscript-figures.Rmd](https://github.com/EmilieFinch/DR-covid19/blob/main/reports/manuscript-figures.Rmd). These can also be generated from model fits and scenario results used in the manuscript, which are saved in the `output` folder. 

Data required for the analysis can be found in the `data` folder and includes:
- Serological data from [Nilles et al, 2021](https://www.sciencedirect.com/science/article/pii/S2667193X22002071)
- COVID-19 death data from the Dominican Republic's COVID-19 dashboard
- COVID-19 case data from [Johns Hopkins Coronavirus Resource Centre](https://github.com/CSSEGISandData/COVID-19)
- Hospital and ICU bed occupancy data, scraped from daily [COVID-19 bulletins](https://www.msp.gob.do/web/?page_id=6948) published by the Ministerio de Salud Pública y Asistencia Social
- Vaccinations from [Oxford's Our World in Data](https://github.com/owid/covid-19-data/tree/master/public/data/vaccinations)
- Shape files for the Caribbean and the Dominican Republic
- [GISAID sequence data](https://gisaid.org/) for the Dominican Republic
- Age-dependent susceptibility and symptomatic rates from [Davies et al 2020](https://pubmed.ncbi.nlm.nih.gov/32546824/)

Please note that data on the clusters shown in Figure 1 are not publicly available and are not included in this repo.
  


 
