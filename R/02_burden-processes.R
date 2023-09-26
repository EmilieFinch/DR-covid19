# -------------------------------------------------------------------------------------------------------
# Set up burden processes
# Author: Emilie Finch, adapted from code by Nick Davies in https://github.com/nicholasdavies/newcovid
# -------------------------------------------------------------------------------------------------------

# Infection-fatality rate -------------------------------------------------------------------------------
# From O'Driscoll et al:https://doi.org/10.1038/s41586-020-2918-0

IFR = fread(
  "Age, IFR
0-4,0.003
5-9,0.001
10-14,0.001
15-19,0.003
20-24,0.006
25-29,0.013
30-34,0.024
35-39,0.04
40-44,0.075
45-49,0.121
50-54,0.207
55-59,0.323
60-64,0.456
65-69,1.075
70-74,1.674
75-79,3.203
80+,8.292")

# Reformat so highest age category is 75+ weighting by age distribution used in O'Driscoll paper

IFR[16,2] <- IFR[16,2]*(133403/(133403+175611)) + IFR[17,2]*(175611/(133403+175611))
IFR[16,1] <- "75+"
IFR <- IFR[-17,]

ifr_odriscoll <- as.numeric(unlist(IFR[,2]/100))

# Infection-hospitalisation rate ---------------------------------------------------------------------------

# From Esposito et al: https://doi.org/10.1186/s12879-022-07262-0

isr_esposito <- exp(-7.28 + 0.076*0:85)/(1 + exp(-7.28 + 0.076*0:85))
ISR_ICR <- fread(
  "Age,isr,icr
0-4,	0.086,	0.0069
5-9,	0.13,	0.011
10-14,	0.18,	0.019
15-19,	0.27,	0.03
20-24,	0.39,	0.05
25-29,	0.56,	0.081
30-34,	0.82,	0.13
35-39,	1.2,	0.22
40-44,	1.7,	0.36
45-49,	2.5,	0.59
50-54,	3.7,	0.96
55-59,	5.2,	1.6
60-64,	7.5,	2.6
65-69,	10.5,	4.1
70-74,	14.6,	6.6
75-79,	20,	10.4
80-84,	26.6,	15.9
85+,	34.3,	23.4")

# Reformat so highest age category is 75+ 

ISR_ICR[16,2] <- sum(as.numeric(ISR_ICR[16,2]), as.numeric(ISR_ICR[17,2]), as.numeric(ISR_ICR[18,2]))/3
ISR_ICR[16,3] <- sum(as.numeric(ISR_ICR[16,3]), as.numeric(ISR_ICR[17,3]), as.numeric(ISR_ICR[18,3]))/3
ISR_ICR[16,1] <- "75+"
ISR_ICR <- ISR_ICR[-(17:nrow(ISR_ICR)),]

ISR_ICR[,2:3] <- ISR_ICR[,2:3] /100

# Amalgamate probabilities

probabilities = ISR_ICR
probabilities$ifr = ifr_odriscoll

# Create model burden processes ------------------------------------------------------------------------------

 P.hosp     = probabilities[, isr];
 P.critical = probabilities[, icr];
 # P.severe   = probabilities[, ihr * (1 - picu)];
 P.death    = probabilities[, ifr];

burden_processes = list(
  
    # Variant 1 burden processes
  
    list(source = "newI", type = "multinomial", names = c("to_hosp", "null"), report = c("", ""),
        prob = matrix(c(P.hosp, 1 - P.hosp), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(6.0 + 2.5, 0.71, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "newI", type = "multinomial", names = c("to_icu", "null"), report = c("", ""),
        prob = matrix(c(P.critical, 1 - P.critical), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(9.6 + 2.5, 1.91, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "to_hosp", type = "multinomial", names = "hosp", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(11.08, 1.202, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "to_icu", type = "multinomial", names = "icu", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(13.33, 1.25, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "newI", type = "multinomial", names = c("death", "null"), report = c("o", ""),
        prob = matrix(c(P.death, 1 - P.death), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_lnorm(15, 0.9, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    # Sero and PCR positivity
    
    list(source = "newII2", type = "multinomial", names = "to_lfia_positive", report = c(""),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(delay_normal(11.9 + 2.5, 5.3, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "to_lfia_positive", type = "multinomial", names = "lfia_positive", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(c(rep(0,799*4),1), nrow = 1, byrow = T)), 

    list(source = "to_hosp", type = "multinomial", names = "hosp_undetected", report = c("po"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(6, 1, 60, 0.25)$p, nrow = 1, byrow = T)),
    
    list(source = "newEE2", type = "multinomial", names = "to_pcr_positive", report = c(""),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(2.76, 4.79, 60, 0.25)$p), nrow = 1, byrow = T)),

    list(source = "to_pcr_positive", type = "multinomial", names = "pcr_positive", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(8.47, 1.96, 730, 0.25)$p, nrow = 1, byrow = T)),
    
    # Variant 2 burden processes
    
    list(source = "newI2", type = "multinomial", names = c("to_hosp2", "null"), report = c("", ""),
        prob = matrix(c(P.hosp, 1 - P.hosp), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(6.0 + 2.5, 0.71, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "newI2", type = "multinomial", names = c("to_icu2", "null"), report = c("", ""),
        prob = matrix(c(P.critical, 1 - P.critical), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(9.6 + 2.5, 1.91, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "to_hosp2", type = "multinomial", names = "hosp2", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(11.08, 1.202, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "to_icu2", type = "multinomial", names = "icu2", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(13.33, 1.25, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "newI2", type = "multinomial", names = c("death2", "null"), report = c("o", ""),
        prob = matrix(c(P.death, 1 - P.death), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_lnorm(15, 0.9, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),
    
    list(source = "to_hosp2", type = "multinomial", names = "hosp_undetected2", report = c("po"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(6, 1, 60, 0.25)$p, nrow = 1, byrow = T)),
    
    list(source = "newI", type = "multinomial", names = "test", report = c("o"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(2.5 + 3.9, 2, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "newI2", type = "multinomial", names = "test2", report = c("o"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(2.5 + 3.9, 2, 60, 0.25)$p, nrow = 1, byrow = T))
)

params$processes = burden_processes

rm(IFR, ISR_ICR, probabilities)
