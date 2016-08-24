# ====================================================
# Arquivo de configuracao global do Alerta Dengue
# ====================================================

#------- Pacotes para carregar 

library(foreign)
library("RPostgreSQL")
library(xtable)
library(knitr)
#library(AlertTools)
devtools::load_all()
source("../config/codigofiguras.R")
source("../config/configRelatorio.r")

# ====================================================
## Parametros globais do alerta
# ====================================================

## -------- Distribuicao do tempo de geracao da dengue:
gtdist="normal"; meangt=3; sdgt = 1.2   


## --------- Regras de mudança de nivel de alerta
# (criterio, duracao da condicao para turnon, turnoff)
criteria = list(
#crity = c("temp_min > tcrit | (temp_min < tcrit & inc > preseas)", 3, 2),
crity = c("temp_min > tcrit", 3, 1),
crito = c("p1 > 0.95 & inc > preseas & temp_min >= tcrit", 3, 1),
critr = c("inc > inccrit", 2, 2)
)
