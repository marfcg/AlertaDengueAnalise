# ====================================================
# Arquivo de configuracao global do Alerta Dengue
# ====================================================

#------- Pacotes para carregar 

pkgs <- c("foreign", "forecast", "RPostgreSQL", "xtable",
          "zoo","tidyverse","assertthat","AlertTools",
          "futile.logger")

lapply(pkgs, library, character.only = TRUE ,quietly = T)

# explicita o diretorio base
basedir = system("pwd", intern=TRUE)
print(basedir)

# codigos das figuras dos boletins
source("AlertaDengueAnalise/config/codigofiguras.R")
source("AlertaDengueAnalise/config/configRelatorio.r")

# INLA
INLA:::inla.dynload.workaround()
# ====================================================
## Parametros globais do alerta
# ====================================================

## -------- Distribuicao do tempo de geracao da dengue:
#gtdist="normal"; meangt=3; sdgt = 1.2   
#GT.Deng = list(gtdist="normal", meangt=3, sdgt = 1.2) 
## --------- Distri tempode geracao da chikungunya
#The mean serial interval ranged from 1.5 to 2.7 
# weeks for CHIKV and from 2.2 to 4.7 weeks for ZIKV according to the change in temperature over the periods
# Riou et al 2017
#GT.Chik = list(gtdist="normal", meangt=2, sdgt = 1)   

## --------- Regras de mudança de nivel de alerta
# (criterio, duracao da condicao para turnon, turnoff)
#criteria = list(
#crity = c("temp_min > tcrit | (temp_min < tcrit & inc > preseas)", 3, 2),
#crito = c("p1 > 0.9 & inc > preseas", 2, 2),
#critr = c("inc > inccrit", 1, 2)
#)

### DENGUE  ### 

#criteria = list(
  #crity = c("temp_min > tcrit | (temp_min < tcrit & inc > preseas)", 3, 2),
#  crity = c("temp_min > tcrit & casos > 0", 3, 1),
#  crito = c("p1 > 0.95 & inc > preseas", 3, 1),
#  critr = c("inc > inccrit & casos > 5", 2, 2)
#)

#criteriaCE = list(
  #crity = c("temp_min > tcrit | (temp_min < tcrit & inc > preseas)", 3, 2),
#  crity = c("umid_max > ucrit & casos > 0", 3, 1),
#  crito = c("p1 > 0.95", 3, 1),
#  critr = c("inc > inccrit & casos > 5", 2, 2)
#)

### CHIK ###

#criteriaCE.chik = list(
#  #crity = c("temp_min > tcrit | (temp_min < tcrit & inc > preseas)", 3, 2),
#  crity = c("umid_max > ucrit & casos > 0", 3, 1),
#  crito = c("p1 > 0.95", 3, 1),
#  critr = c("inc > inccrit & casos > 5", 2, 2)
#)

#criteriaChik = list(
#  crity = c("temp_min > tcrit & casos > 0", 3, 1),
#  crito = c("p1 > 0.95 & temp_min >= tcrit", 3, 1),
#  critr = c("inc > inccrit & casos > 5", 2, 2)
#)

## Limiares epidemicos
## Estao inscritos na tabela regional_saude
## Refazer anualmente


