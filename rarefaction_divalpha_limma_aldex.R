
library(tidyverse)
library(readxl)
library(edgeR)
library(Glimma)
library(ALDEx2)
library(cluster)

pracma::clear()


#### IMPORTAR Y LIMPIAR LOS DATOS ####

desktfg <- readRDS(file = '../Data_TFG/desk_deskdefault/fastqofchoice_counts_DESKTOP.RDS')


#### RAREFACCIÓN ####

rarelist <- vegan::rarecurve(desktfg %>% t() %>% as.data.frame(),
                             step = 100)

maxsam <- sapply(rarelist, length) %>% which.max()  # ¿qué número de muestra tiene más puntos?
maxptos <- rarelist[[maxsam]] %>% length()  # ¿cuántos tiene?
max_samsize <- attr(rarelist[[maxsam]], 'Subsample')  # valores de abscisa, tamaño de submuestreo de 100 en 100

# La idea es crear un dataframe con los tamaños de submuestra (eje X) en una columna
# y las especies halladas en cada submuestra (eje Y) en las demás.
# Adecuado para ggplot(), para representar cada grupo en un gráfico de una rejilla.
# Rellenar con NAs para igualar las diferencias de puntos.

raredf <- matrix(ncol = ncol(desktfg), nrow = maxptos) %>% as.data.frame()
colnames(raredf) <- colnames(desktfg)
for (sam in 1:ncol(desktfg)) {
  ptos <- rarelist[[sam]] %>% as.numeric()
  if (length(ptos) < maxptos) {
    raredf[, sam] <- c(
      ptos,
      rep(NA, maxptos-length(ptos))
    )
  } else {
    raredf[, sam] <- ptos
  }
}
raredf <- bind_cols(raredf, 
                    data.frame('sample_size' = max_samsize))

# Visión general
raredf %>% 
  pivot_longer(cols = -sample_size, 
               names_to = 'samID', 
               values_to = 'nspecies') %>%
  mutate('kk1' = samID) %>%
  separate(col = kk1,
           sep = '_', 
           into = c('group', 'kk2')) %>% 
  dplyr::select(-starts_with('kk')) %>% 
  ggplot(aes(x = sample_size, y = nspecies)) +
  geom_point() +
  facet_wrap(~group, scales = 'free')

# Grupo por grupo
# pBIC
raredf %>% 
  pivot_longer(cols = -sample_size, 
               names_to = 'samID', 
               values_to = 'nspecies') %>%
  mutate('kk1' = samID) %>%
  separate(col = kk1,
           sep = '_', 
           into = c('group', 'kk2')) %>% 
  dplyr::select(-starts_with('kk')) %>% 
  filter(group == 'pBIC') %>% 
  ggplot(aes(x = sample_size, y = nspecies, colour = samID)) +
  geom_point()
# pBIC10
raredf %>% 
  pivot_longer(cols = -sample_size, 
               names_to = 'samID', 
               values_to = 'nspecies') %>%
  mutate('kk1' = samID) %>%
  separate(col = kk1,
           sep = '_', 
           into = c('group', 'kk2')) %>% 
  dplyr::select(-starts_with('kk')) %>% 
  filter(group == 'pBIC10') %>% 
  ggplot(aes(x = sample_size, y = nspecies, colour = samID)) +
  geom_point()


#### CRIBAS ####

# Se establece que cada OTU tenga, como mínimo, tantas cuentas como muestras haya en el grupo más pequeño.

colnames(desktfg) %>% 
  as.data.frame() %>% 
  separate(col = '.',
           into = c('genotype', 'sample'),
           sep = '_') %>% 
  pull('genotype') %>% 
  table()

colnames(desktfg) %>% 
  as.data.frame() %>% 
  separate(col = '.',
           into = c('genotype', 'sample'),
           sep = '_') %>% 
  pull('genotype') %>% 
  table() %>% 
  as.numeric() -> tam_grupos

desktfg_criba <- desktfg[rowSums(desktfg) >= min(tam_grupos), ]


#### DIVERSIDAD ALFA ####

Simpson <- numeric(length = ncol(desktfg_criba))
InvSimpson <- numeric(length = ncol(desktfg_criba))
for(i in seq_along(colnames(desktfg_criba))){
  D <- 0
  for(j in 1:nrow(desktfg_criba)){
    D <- D + (desktfg_criba[j,i]/sum(desktfg_criba[,i]))^2
  }
  Simpson[i] <- 1-D
  InvSimpson[i] <- 1/D
}

Shannon <- numeric(length = ncol(desktfg_criba))
for(i in seq_along(colnames(desktfg_criba))){
  H <- 0
  for(j in 1:nrow(desktfg_criba)){
    p <- desktfg_criba[j,i]/sum(desktfg_criba[,i])
    if(p == 0){
      H <- H
    }
    else{
      H <- H + p*log(p)
    }
  }
  Shannon[i] <- H*(-1)
}

Observed <- numeric(length = ncol(desktfg_criba))
for(i in seq_along(colnames(desktfg_criba))){
  n_especies <- 0
  for(j in 1:nrow(desktfg_criba)){
    if(desktfg_criba[j,i] > 0){
      n_especies <- n_especies + 1
    }
  }
  Observed[i] <- n_especies
}

diversidad_alfa <- data.frame(Sample_ID = colnames(desktfg_criba),
                              Auxiliar = colnames(desktfg_criba)) %>% 
  separate(Auxiliar, into = c('Genotype', 'kk'), sep = "_") %>% 
  dplyr::select(Sample_ID, Genotype) %>% 
  mutate("Simpson" = Simpson, "InvSimpson" = InvSimpson, 
         "Shannon" = Shannon, "Observed" = Observed)

rm(D, H, InvSimpson, Observed, Shannon, Simpson)

diversidad_alfa %>% 
  pivot_longer(cols = c('Simpson', 'Shannon', 'InvSimpson', 'Observed'),
               names_to = 'Index_type',
               values_to = 'Alpha_diversity') %>%
  ggplot(aes(x = Genotype, y = Alpha_diversity, color = Genotype)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~Index_type, scales = 'free') +
  theme_minimal()



#### ANÁLISIS DE ABUNDANCIA DIFERENCIAL: limma ####

group <- as.factor(diversidad_alfa$Genotype)
group

# Normalize with TMMwsp:
dge <- DGEList(desktfg_criba, group = group)
dim(dge)
dge <- normLibSizes(dge, method = "TMMwsp")
dge$samples
# Limma
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
v <- voom(dge, design, plot = TRUE)
head(v$E)
# Fit linear model
fit<-lmFit(v)
names(fit)
# Create contrast matrix (for one comparison)
cont.matrix <- makeContrasts(pBIC_pBIC10 = pBIC10 - pBIC, 
                             levels=design)
cont.matrix
# Run statistics between contrasts and eBayes shrinkage:
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
# Check number of DA OTUs with different adjusted p-value thresholds
summary(decideTests(fit.cont, p.value = 0.1, lfc =1))
#
# Now get top DA OTUs for the contrast with desired values of adjusted p-value and logFC.
result.limma <- topTable(fit.cont, coef="pBIC_pBIC10", p.value = 0.1, lfc = 1, 
                         sort.by="logFC", resort.by = "logFC", number = 400)
dim(result.limma)
# Check heatmap
daos.limma <- rownames((result.limma))
heatgenes <- (log.counts.cpm[daos.limma,])
pheatmap::pheatmap(heatgenes, scale = "row", fontsize_row = 7, 
                   annotation_col = dplyr::select(diversidad_alfa, c(Sample_ID, Genotype)) %>% 
                     column_to_rownames('Sample_ID'),
                   legend = T)


#### ANÁLISIS DE ABUNDANCIA DIFERENCIAL: ALDEx2 ####

grupo1 <- diversidad_alfa %>% dplyr::select(c(Sample_ID, Genotype)) %>% filter(Genotype == 'pBIC10')
grupo2 <-  diversidad_alfa %>% dplyr::select(c(Sample_ID, Genotype)) %>% filter(Genotype == 'pBIC')
counts.aldex <- desktfg_criba[,c(grupo1$Sample_ID, grupo2$Sample_ID)]
colnames(counts.aldex)
conds <- c(grupo1$Genotype, grupo2$Genotype)
str(conds)
# Run ALDEx2:
x.all <- ALDEx2::aldex(counts.aldex, conds, mc.samples=256, test="t", effect=TRUE,
                       include.sample.summary=FALSE, denom="all", verbose=TRUE, paired.test=FALSE)
# Get results (DAOs below an FDR = 10%):
result.aldex <- x.all[x.all$we.eBH < 0.1,]
daos.aldex <- rownames(result.aldex)
heatgenes <- (log.counts.cpm[daos.aldex,])
pheatmap::pheatmap(heatgenes, scale = "row", fontsize_row = 7, 
                   annotation_col = dplyr::select(diversidad_alfa, c(Sample_ID, Genotype)) %>% 
                     column_to_rownames('Sample_ID'),
                   legend = T)

#### COMPARACIÓN ENTRE METODOLOGÍAS ####
daos.final <- intersect(daos.aldex, daos.limma)
heatgenes <- (log.counts.cpm[daos.final,])
pheatmap::pheatmap(heatgenes, scale = "row", fontsize_row = 10, 
                   annotation_col = dplyr::select(diversidad_alfa, c(Sample_ID, Genotype)) %>% 
                     column_to_rownames('Sample_ID'),
                   legend = T)
result.limma[daos.final,] %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.val))) +
  geom_point(size=3) +
  xlim(-max(result.limma[daos.final,]$logFC),max(result.limma[daos.final,]$logFC))
