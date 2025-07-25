library(readxl)
library(tidyverse)
library(FactoMineR)
library(ggthemes)
library(psych)
library(NbClust)
library(tidyverse)
library(GPArotation)
library(fpc)
library(cluster)
library(NbClust)
library(factoextra)

### Dataset
vote <- read_excel("Voting Motivations Data Set.xlsx") 
vote2 <- vote %>% select(Column1, Column2, Column3, Column4, Column5, Column6, Column7, Column8, Column9, Column10, Column11, Column12, Column13, Column14, Column15, Column16, Column17, Column18, Column19, Column20, Column21, Column22, Column23, Column24)
vote_motive <- as.data.frame(na.omit(vote2))

standardized <- scale(vote_motive, center=T, scale=T)
apply(standardized, 2, mean)    # Should be close to 0
apply(standardized, 2, sd) 

### Determining Sampling Adequacy
# Test for Sphericity
cortest.bartlett(R=cor(standardized), n=190)
KMO(r = standardized)

# The test for sphericity shows a very low p-value that is near 0, then we reject the null hypothesis that the variables are uncorrelated. Hence, factor analysis is appropriate for the data. Also, we will remove columns 16, 19, and 21 since they have MSA\<0.7.

### Filtered data according to KMO

vote3 <- vote %>% select(Column1, Column2, Column3, Column4, Column5, Column6, Column7, Column8, Column9, Column10, Column11, Column12, Column13, Column14, Column15, Column17, Column18, Column22, Column23, Column24)
real_vote_motive <- as.data.frame(na.omit(vote3))
standardized_motive <- scale(real_vote_motive, center=T, scale=T)

# Test for Sphericity
cortest.bartlett(R=cor(standardized_motive), n=190)
KMO(r = standardized_motive)
# We have obtained a KMO of 0.9 meaning that the results are interpreted as marvelous. Hence, it is appropriate to use FA. Moreover, there are no variables with MSA \< 0.7

### Determining number of factors
pca <- prcomp(x=standardized_motive)
pca_summary <- summary(pca)

squared_sd <- pca$sdev^2

result <- data.frame(
  PC = paste0("PC", seq_along(squared_sd)),
  Standard_Deviation = pca$sdev,
  variance = squared_sd,
  Proportion_of_Variance = pca_summary$importance[2, ],
  Cumulative_Proportion = pca_summary$importance[3, ])

print(result)
# If we set our cutoff at 0.7 for the cumulative proportion of variance, we are going to use 7 factors for our factor analysis. Using Kaiser's rule, we only need 7 factors.

## Scree Plots
tibble(eigenvalues=(pca$sdev)^2, PC=1:20) %>% 
  ggplot(aes(y=eigenvalues, x=PC)) +
  geom_point() +
  geom_bar(stat="identity", fill="skyblue", alpha=0.5) +
  geom_line() +
  scale_x_continuous(breaks=1:20) +
  ggthemes::theme_gdocs()

# However, using scree plot, the most noticeable elbow point is the one where we only have 4 factors.

### Factor Loading
as_tibble(round(pca$rotation, digits=4), rownames = "variable")

### Very Simple Structure
vss(standardized_motive)

### Doing the factor analysis
final_vote <- vote %>% select(Column1, Column2, Column3, Column4, Column5, Column6, Column7, Column8, Column9, Column12, Column13, Column14, Column15, Column17, Column18)
final_vote_motive <- as.data.frame(na.omit(final_vote))
final_standardized_motive <- scale(final_vote_motive)

fa <- psych::fa(r = final_standardized_motive,
                nfactors = 4,
                rotate = "Promax",
                scores = "regression",
                SMC = T,
                fm = "pa")
fa.sort(fa)

cortest.bartlett(R=cor(final_standardized_motive), n=190)

KMO(r = final_standardized_motive)

factors <- as_tibble(fa$scores)

euc_dist <- cluster::daisy(x = cluster_final_standardized_motive, metric = "euclidean")
str(euc_dist)

### HIERARCHICHAL ###
pc <- prcomp(fa$scores)
summary(pc)
pc$rotation

subset <- cluster_final_standardized_motive %>%
  bind_cols(pc$x) %>%
  select(PC1, PC2)

dist_from_center <- mahalanobis(x = subset, center = c(0, 0), cov = var(subset))

cluster_final_standardized_motive %>%
  bind_cols(pc$x) %>%
  mutate(id = row_number(),
         dist = dist_from_center,
         tag = if_else(dist > qchisq(0.999, 2), "might be outlier", "ok")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_text(aes(label = id, col = tag)) +
  ggthemes::theme_gdocs()

# Creating a Function for the Mahalanobis Distance
mahal <- function(x, cx = NULL) {
  x <- as.data.frame(x)
  if(is.null(cx)) cx <- cov(x)
  out <- lapply(1:nrow(x), function(i) {
    mahalanobis(x = x, 
                center = do.call("c", x[i, ]),
                cov = cx,
                tol=1e-20)
  })
  return(as.dist(do.call("rbind", out)))
}

mahal_dist <- mahal(cluster_final_standardized_motive)
agnes_mahal_comp<- cluster::agnes(x = mahal_dist, diss = T, method = "ward")
plot(agnes_mahal_comp, which.plots = 2)

survey_mahal <- mahal(fa$scores)

### Ward's and Mahalanobis
best_ch <- NbClust(data = factors,
                   distance = NULL,
                   method = "ward.D2",
                   index = "ch")

best_silhouette <- NbClust(data = factors,
                           diss = survey_mahal,
                           distance = NULL,
                           method = "ward.D2",
                           index = "silhouette")

best_db <- NbClust(data = factors,
                   diss = survey_mahal,
                   distance = NULL,
                   method = "ward.D2",
                   index = "db")

best_ratkowsky <- NbClust(data = factors,
                          diss = survey_mahal,
                          distance = NULL,
                          method = "ward.D2",
                          index = "ch")

best_kl <- NbClust(data = factors,
                   diss = survey_mahal,
                   distance = NULL,
                   method = "ward.D2",
                   index = "kl")

best_ptbiserial <- NbClust(data = factors,
                           diss = survey_mahal,
                           distance = NULL,
                           method = "ward.D2",
                           index = "ptbiserial")

best_gap <- NbClust(data = factors,
                    diss = survey_mahal,
                    distance = NULL,
                    method = "ward.D2",
                    index = "gap")


best_frey <- NbClust(data = factors,
                     diss = survey_mahal,
                     distance = NULL,
                     method = "ward.D2",
                     index = "frey")

best_mcclain <- NbClust(data = factors,
                        diss = survey_mahal,
                        distance = NULL,
                        method = "ward.D2",
                        index = "mcclain")

best_gamma <- NbClust(data = factors,
                      diss = survey_mahal,
                      distance = NULL,
                      method = "ward.D2",
                      index = "gamma")

best_gplus <- NbClust(data = factors,
                      diss = survey_mahal,
                      distance = NULL,
                      method = "ward.D2",
                      index = "gplus")

best_dunn <- NbClust(data = factors,
                     diss = survey_mahal,
                     distance = NULL,
                     method = "ward.D2",
                     index = "dunn")

best_sdindex <- NbClust(data = factors,
                        diss = survey_mahal,
                        distance = NULL,
                        method = "ward.D2",
                        index = "sdindex")

best_sdbw <- NbClust(data = factors,
                     diss = survey_mahal,
                     distance = NULL,
                     method = "ward.D2",
                     index = "sdbw")

best_cindex <- NbClust(data = factors,
                       diss = survey_mahal,
                       distance = NULL,
                       method = "ward.D2",
                       index = "cindex")

best_hartigan <- NbClust(data = factors,
                         diss = survey_mahal,
                         distance = NULL,
                         method = "ward.D2", 
                         index = "hartigan")

best_scott <- NbClust(data = factors,
                      diss = survey_mahal,
                      distance = NULL,
                      method = "ward.D2", 
                      index = "ball")

best_tau <- NbClust(data = factors,
                      diss = survey_mahal,
                      distance = NULL,
                      method = "ward.D2", 
                      index = "tau")
best_ch$Best.nc
best_dunn$Best.nc
best_frey$Best.nc
best_gamma$Best.nc
best_gap$Best.nc
best_kl$Best.nc
best_ptbiserial$Best.nc
best_ratkowsky$Best.nc
best_silhouette$Best.nc
best_tau$Best.nc
best_db$Best.nc
best_gplus$Best.nc
best_mcclain$Best.nc
best_sdbw$Best.nc
best_sdindex$Best.nc
best_sdindex$Best.nc
best_hartigan$Best.nc
best_scott$Best.nc

### Ward's and Euclidean
best_ch <- NbClust(data = factors,
                   distance = "euclidean",
                   method = "ward.D2",
                   index = "ch")

best_silhouette <- NbClust(data = factors,
                           distance = "euclidean",
                           method = "ward.D2",
                           index = "silhouette")

best_db <- NbClust(data = factors,
                   distance = "euclidean",
                   method = "ward.D2",
                   index = "db")

best_ratkowsky <- NbClust(data = factors,
                          distance = "euclidean",
                          method = "ward.D2",
                          index = "ch")

best_kl <- NbClust(data = factors,
                   distance = "euclidean",
                   method = "ward.D2",
                   index = "kl")

best_ptbiserial <- NbClust(data = factors,
                           distance = "euclidean",
                           method = "ward.D2",
                           index = "ptbiserial")

best_gap <- NbClust(data = factors,
                    distance = "euclidean",
                    method = "ward.D2",
                    index = "gap")


best_frey <- NbClust(data = factors,
                     distance = "euclidean",
                     method = "ward.D2",
                     index = "frey")

best_mcclain <- NbClust(data = factors,
                        distance = "euclidean",
                        method = "ward.D2",
                        index = "mcclain")

best_gamma <- NbClust(data = factors,
                      distance = "euclidean",
                      method = "ward.D2",
                      index = "gamma")

best_gplus <- NbClust(data = factors,
                      distance = "euclidean",
                      method = "ward.D2",
                      index = "gplus")

best_dunn <- NbClust(data = factors,
                     distance = "euclidean",
                     method = "ward.D2",
                     index = "dunn")

best_sdindex <- NbClust(data = factors,
                        distance = "euclidean",
                        method = "ward.D2",
                        index = "sdindex")

best_sdbw <- NbClust(data = factors,
                     distance = "euclidean",
                     method = "ward.D2",
                     index = "sdbw")

best_cindex <- NbClust(data = factors,
                       distance = "euclidean",
                       method = "ward.D2",
                       index = "cindex")

best_hartigan <- NbClust(data = factors,
                         distance = "euclidean",
                         method = "ward.D2", 
                         index = "hartigan")

best_scott <- NbClust(data = factors,
                      distance = "euclidean",
                      method = "ward.D2", 
                      index = "ball")

best_tau <- NbClust(data = factors,
                      distance = "euclidean",
                      method = "ward.D2", 
                      index = "tau")
best_ch$Best.nc
best_dunn$Best.nc
best_frey$Best.nc
best_gamma$Best.nc
best_gap$Best.nc
best_kl$Best.nc
best_ptbiserial$Best.nc
best_ratkowsky$Best.nc
best_silhouette$Best.nc
best_tau$Best.nc
best_db$Best.nc
best_gplus$Best.nc
best_mcclain$Best.nc
best_sdbw$Best.nc
best_sdindex$Best.nc
best_sdindex$Best.nc
best_hartigan$Best.nc
best_scott$Best.nc

ward_cluster6 <- cutree(tree = agnes_mahal_comp, k = 5)
table(ward_cluster6)

ward_cluster8 <- cutree(tree = agnes_ward, k = 5)
table(ward_cluster8)

## Euclidean
agnes_ward <- cluster::agnes(x = cluster_final_standardized_motive, metric ="euclidean", method = "ward")
plot(agnes_ward, which.plots = 2)

ward_clust5 <- cutree(tree = agnes_ward, k = 6)
table(ward_clust5)

## Manhattan
man_agnes_ward <- cluster::agnes(x = cluster_final_standardized_motive, metric ="manhattan", method = "ward")
plot(man_agnes_ward, which.plots = 2)

ward_clust5 <- cutree(tree = man_agnes_ward, k = 6)
table(ward_clust5)

as_tibble(factors) %>%
  mutate(cluster = ward_cluster6) %>%
  group_by(cluster) %>%
  summarise_all(.funs = mean) %>%
  rename(Aksesibilidad = PA4,
         Tindig = PA2,
         Sistema = PA3,
         Pangsarili = PA1)

fviz_gap_stat(clusGap(factors,
                       FUN = clara,
                       K.max = 6,
                       B = 100))

fviz_nbclust(factors, clara, method = "wss")
fviz_nbclust(factors, clara, method = "silhouette")
fviz_nbclust(factors, clara, method = "gap_stat")
