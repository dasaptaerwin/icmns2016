# Kode R untuk menganalisis data petrografi batuan
# Judul skripsi: ....
# Penulis: ....
# Kode ini dikembangkan dari berbagai contoh dari situs: www.....

#reading data

df <- read.csv("20161005_data_volcanic_Jejen.csv") 

row.names(df) <- df$Kode # setting row names

df2 <- df[4:8] # exclude Location name
install.packages("tidyverse")
library("tidyverse")

install.packages("ggcorrplot")
library('ggcorrplot')

install.packages('cluster')
library(cluster) # for cluster analysis

install.packages('FactoMineR')
library('FactoMineR')

install.packages('factoextra')
library('factoextra')

install.packages('vegan')
library(vegan)

cor.tab <- cor(df2)
cor.tab

ggcorrplot(cor.tab)              # making heatmap

## locating and imputing missing data
install.packages('mice')
library('mice')
md.pattern(df2) # from package 'mice'
df2 <- mice(df2,m=5,maxit=50,meth='pmm',seed=500)
df2 <- as.data.frame(df2)
  
# Run PCA 
res.pca <- PCA(df2, graph = FALSE)
res.pca <- prcomp(df2, center = TRUE, scale = TRUE, na.action = na.omit)
res.pca
fviz_pca(res.pca, choix = "var", col.var="contrib", jitter=list(what = "label", width =NULL, height = NULL))

fviz_pca(res.pca, choix = "var", col.var="contrib", repel=TRUE)

################# old lines ################

summary(jejen)
str(jejen)

is.numeric(jejen$Plagioklas)
is.numeric(jejen$opak)
is.numeric(jejen$Piroksen)
is.numeric(jejen$olivin)
is.numeric(jejen$gelas)
is.numeric(jejen$Xenolith)
is.numeric(jejen$Nephelin)
is.numeric(jejen$Kuarsa)
is.numeric(jejen$Tuff)
is.numeric(jejen$Fragmen.batuan)
is.numeric(jejen$Kalsit)
is.numeric(jejen$Hornblenda)
is.numeric(jejen$Porositas)
is.numeric(jejen$Mineral.Lempung)
is.numeric(jejen$Oksida.besi)
is.numeric(jejen$Klorit)
is.numeric(jejen$Pseudomorf)
is.numeric(jejen$K.Feldspar)
is.numeric(jejen$Foraminifera)
is.numeric(jejen$Ganggang)


jejen2 = jejen[4:8]
cor(jejen2, use="complete.obs", method="kendall")
jejen2 = jejen[9:15]
cor(jejen2, use="pairwise.complete.obs", method="kendall")
jejen2 = jejen[16:17]
cor(jejen2, use="pairwise.complete.obs", method="kendall")
jejen2 = jejen[18:23]
cor(jejen2, use="pairwise.complete.obs", method="kendall")


# Prepare Data
jejen2 <- na.omit(jejen2) # listwise deletion of missing
jejen2 <- scale(jejen2) # standardize variables 

#cov(mydata, use="complete.obs") 
pairs(jejen2)
pairs(jejen)

#PLOT
plot(jejen2$Plagioklas~jejen2$opak)
plot(jejen2$Plagioklas)
hist(jejen2$Plagioklas)
hist(jejen2$Plagioklas, col = "red",xlab = "nourut",ylab = "persentase", main = "plagioklas")

# Prepare Data
jejen2 <- na.omit(jejen2) # listwise deletion of missing
jejen2 <- scale(jejen2) # standardize variables 

# K-Means Cluster Analysis
fit <- kmeans(jejen2, 5) # 5 cluster solution
# get cluster means
aggregate(jejen2,by=list(fit$cluster),FUN=mean)
# append cluster assignment
jejen2 <- data.frame(jejen2, fit$cluster) 


#PCA
pca <- prcomp(jejen2)
summary(pca)
pca$sdev
screeplot(pca, type="lines")


install.packages("devtools")
library(devtools)
install_github("kassambara/factoextra")
library("factoextra")
install.packages("FactoMineR")
install_github("kassambara/FactoMineR")
library("FactoMineR")
pca.res <- PCA(jejen2, graph=FALSE)
print (pca.res)
eigenvalues <- pca.res$eig
head(eigenvalues)
fviz_screeplot(pca.res, ncp=10)
fviz_pca_var(pca.res)
plot(pca.res, choice="ind")
fviz_pca_contrib(pca.res, choice="var", axes=1:2)

library(fpc)
plotcluster(jejen2, fit$cluster)
# vary parameters for most readable graph
library(cluster)
clusplot(jejen2, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
abline(h=0, v=0, col="red")
table(jejen2[,1],fit$cluster)

#hirearki clustering
d <- dist(jejen2, method = "euclidean") # Euclidean distance matrix.
H.fit <- hclust(d, method="ward.D2")

plot(H.fit) # display dendogram
groups <- cutree(H.fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(H.fit, k=5, border="red") 
table(jejen2[,1],groups)

#studi
age<-read.table('VUL.KUARTER.csv', sep = ';', header=T)
head(age)
compositions<-read.table('VUL.KUARTER.CSV', sep = ';', header=T)
head(compositions)

library(cluster)
D=daisy(jejen2, metric='gower')
H.fit <- hclust(D, method="ward.D2")

plot(H.fit) # display dendrogram
groups <- cutree(H.fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters
rect.hclust(H.fit, k=4, border="red") 

clusplot(jejen2, groups, color=TRUE, shade=TRUE,labels=2, lines=0, main= 'MINERAL')
