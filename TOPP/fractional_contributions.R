library(tidyr)
library(dplyr)

tno.data = read.csv(file = "TNO_Fractional_Contributions.csv", header = TRUE)
ipcc.data = read.csv(file = "IPCC_Fractional_Contributions.csv", header = TRUE)
emep.data = read.csv(file = "EMEP_Fractional_Contributions.csv", header = TRUE)
de94.data = read.csv(file = "DE94_Fractional_Contributions.csv", header = TRUE)
gr95.data = read.csv(file = "GR95_Fractional_Contributions.csv", header = TRUE)
gr05.data = read.csv(file = "GR05_Fractional_Contributions.csv", header = TRUE)
uk98.data = read.csv(file = "UK98_Fractional_Contributions.csv", header = TRUE)
uk08.data = read.csv(file = "UK08_Fractional_Contributions.csv", header = TRUE)

extract.mozart.data = function (df, speciation) {
    df = tbl_df(df)
    df = df %>% select(MCM.Species, MOZART.Species, Fractional.Contribution.of.Species)
    df = df %>% rename(Contribution = Fractional.Contribution.of.Species, Mechanism.Species = MOZART.Species)
    df$Speciation = rep(speciation, length(df$MCM.Species))
    df$Mechanism = rep("MOZART", length(df$MCM.Species))
    return(df)
}

extract.radm2.data = function (df, speciation) {
    df = tbl_df(df)
    df = df %>% select(MCM.Species, RADM2.Species, Fractional.Contribution.of.Species.1)
    df = df %>% rename(Contribution = Fractional.Contribution.of.Species.1, Mechanism.Species = RADM2.Species)
    df$Speciation = rep(speciation, length(df$MCM.Species))
    df$Mechanism = rep("RADM2", length(df$MCM.Species))
    return(df)
}

tno.mozart = extract.mozart.data(tno.data, "TNO")
tno.radm2 = extract.radm2.data(tno.data, "TNO")
ipcc.mozart = extract.mozart.data(ipcc.data, "IPCC")
ipcc.radm2 = extract.radm2.data(ipcc.data, "IPCC")
emep.mozart = extract.mozart.data(emep.data, "EMEP")
emep.radm2 = extract.radm2.data(emep.data, "EMEP")
de94.mozart = extract.mozart.data(de94.data, "DE94")
de94.radm2 = extract.radm2.data(de94.data, "DE94")
gr95.mozart = extract.mozart.data(gr95.data, "GR95")
gr95.radm2 = extract.radm2.data(gr95.data, "GR95")
gr05.mozart = extract.mozart.data(gr05.data, "GR05")
gr05.radm2 = extract.radm2.data(gr05.data, "GR05")
uk98.mozart = extract.mozart.data(uk98.data, "UK98")
uk98.radm2 = extract.radm2.data(uk98.data, "UK98")
uk08.mozart = extract.mozart.data(uk08.data, "UK08")
uk08.radm2 = extract.radm2.data(uk08.data, "UK08")

data = rbind(tno.mozart, tno.radm2, ipcc.mozart, ipcc.radm2, emep.mozart, emep.radm2, de94.mozart, de94.radm2, gr95.mozart, gr95.radm2, gr05.mozart, gr05.radm2, uk98.mozart, uk98.radm2, uk08.mozart, uk08.radm2)
write.csv(data, file = "Fractional_Contributions_MOZ_RADM2_of_MCM_species.csv", row.names = FALSE)

data = data %>% select(Mechanism, Speciation, Mechanism.Species, MCM.Species)
write.csv(data, file = "Mapping_MOZ_RADM2_of_MCM_species.csv", row.names = FALSE)
