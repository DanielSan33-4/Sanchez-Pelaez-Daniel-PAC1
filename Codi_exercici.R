# Llibreries necessàries
library(GEOquery)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)

#Es carreguen les dades
load("C:/Users/dd28/OneDrive/Màster/08. Anàlisi de dades òmiques/Unitat 1. Les ómiques/PAC1/Sanchez-Pelaez-Daniel-PAC1/summarized_experiment.Rda")

# Comprovació de les dades
sum(is.na(assay(se)))
colSums(is.na(assay(se)))
rowSums(is.na(assay(se)))

" Es fa un gràfic resum on es comparen els dos grups factoritzats en les diferents variables del dataset.
S'extreu la matriu de dades i es trasposen per poder utiltizar ggplot2"
data_long <- as.data.frame(t(assay(se))) %>%
  # S'extreu la variable resposta i s'afegeix el ID de cada pacient com una columna
  mutate(Condition = colData(se)$Condition, PatientID = rownames(colData(se))) %>%
  pivot_longer(cols = -c(PatientID, Condition), names_to = "Molecule", values_to = "Concentration")
# Es crea un plot per cada variable 
resumdades_graf <- ggplot(data_long, aes(x = Condition, y = Concentration, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ Molecule, scales = "free") +
  theme_minimal() +
  labs(title = "Distribució de biomolècules segons la condició")

# Anàlisi de covariàncies de les dades centrades 
cachexicNum <- scale(t(assay(se)), center = TRUE, scale = FALSE)
apply(cachexicNum, 2, mean)

"Es calcula la matriu de variànces ajustant a dividir per $n$ en lloc de per $(n-1)$ per fer els resultats compatibles en cas que es considerès necessari 
emplear més d'un mètode, com càlcul mitjançant covariànces i/o mitnjançant funcions, com ara $princomp$ i $prcomp$."
n <- dim (se)[2]
S <- cov(cachexicNum)*(n-1)/n

# Es calcula la matriu de correlacions
R <- cor(cachexicNum)

# Càlcul mitjançant matriu de covariàncies
EIG <- eigen(S)

"A continuació es fa la transformació de les dades associada als components principals, 
obtinguda de la multiplicació de la matriu original per la matriu de vectors propis."
eigenVecs1 <- EIG$vectors
PCAS1 <- cachexicNum %*% eigenVecs1
vars1 <- EIG$values/sum(EIG$values)
round(vars1, 3)


"Com que en el següent gràfic ens interessa veure quines pertanyen a pacients amb cachexia i quins a pacients control, 
retallem lleugerament el gràfic per poder ampliar el gruix de les dades, encara que ens doni com a resultat deixar de veure 3 observacions del gràfic."
xlabel <- paste("PCA1 ", round(vars1[1]*100, 2),"%" )
ylabel <- paste("PCA2 ", round(vars1[2]*100,2),"%" )
bgSurv <- colSurv <- ifelse(se$Condition=="cachexic", "red", "blue")
pchSurv <- ifelse(se$Condition=="cachexic",1,2)
plot(PCAS1[,1], PCAS1[,2], main = "Primeres PCAs de biomolècules",
     xlab=xlabel, ylab=ylabel,
     pch=pchSurv, col=colSurv, bg=bgSurv)
legend( "topright"  , inset = c(0.01,0.01), cex =1,
        bty = "y", legend = c("Cachexic", "Control"),
        text.col = c("red", "blue"),
        col = c("red", "blue"),
        pt.bg = c("red","blue"),
        pch = c(1,2)
)

# Càlcul dels components principals mitjançant les funcions princomp i prcomp

# Mètode princomp
PCAS_met2 <- princomp (cachexicNum)
# Mètode prcomp
PCAS_met3 <- prcomp (cachexicNum)

# Amb aquest codi es genera el .txt amb les dades, el deixo comentat perquè no s'executi cada cop
# dades <- read.csv("~/metaboData/Datasets/2024-Cachexia/human_cachexia.csv")
# write.table(dades, file = "dades_cachexia.txt", sep = "\t", quote = FALSE, row.names = FALSE)

