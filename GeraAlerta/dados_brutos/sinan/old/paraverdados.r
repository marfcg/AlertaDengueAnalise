library(foreign)
d <- read.dbf("Dengue2013_BancoSINAN30_12_2013.dbf")
names(d)
sort(table(d$NM_BAIRRO))
d <- read.dbf("../DENGON2015_03_02_2015.dbf")
sort(table(d$NM_BAIRRO))
d <- read.dbf("Dengue2010_BancoSINAN16_04_2012_v09012014.dbf")
names(d)
plot(d$X,d$Y)
