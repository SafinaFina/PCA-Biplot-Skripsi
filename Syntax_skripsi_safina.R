##CALL DATA##
data=read.csv(file.choose(),sep=";",dec=",",header=T)
head(data)
str(data)

##PREPARE DATA##
Y <- as.matrix(data[,2:9])
X <- scale(X, center = TRUE, scale = FALSE)
head(X)
column_variance <- apply(X, 2, var)

##PACKAGES##
library(ggplot2)
library(reshape2)
library(rrcov)
library(rgl)
library(writexl)

##SINGULAR VALUE DECOMPOSITION##
Xt <- t(X)
XtX <- Xt%*%X
pca <- eigen(cov)
eigenvalues <- pca$values
eigenvectors <- pca$vectors

L <- diag(sqrt(eigenvalues)) 
A <- eigenvectors
U <- X %*% A %*% solve(L)

# Variance Explained
total_variance <- sum(eigenvalues)
variance_explained <- eigenvalues / total_variance
cumulative_variance <- cumsum(variance_explained)

##H and G matrixes##
alpha <- 0.5
G  <- U%*%L^alpha
Ht <- L^(1-alpha)%*%(A)   
H <- t(Ht)

##Korelasi PC dan variabel utama##
cov_matrix <- CH$cov 
var_matrix <- diag(cov)
sqrt_nilaieigen=sqrt(eigenvalues)
rho_matrix = t(t(eigenvectors) * sqrt_nilaieigen) / (sqrt_var)##BIPLOT 2D & 3D##
# VISUALISASI 3D
install.packages("rgl")
library(rgl)
# Analisis Biplot
G12 <- G[,c(1,2)]
H12 <- H[,c(1,2)]
G13 <- G[,c(1,3)]
H13 <- H[,c(1,3)]
G23 <- G[,c(2,3)]
H23 <- H[,c(2,3)] 

# PC 1 VS PC 2
biplot(G12,H12,cex=0.8,main="PC1 vs PC2",xlab="PC 1
(40,7%)",ylab="PC2 (22,5%)")
abline(h=0)
abline(v=0)
# PC 1 VS PC 3
biplot(G13,H13,cex=0.8,main="PC1 vs PC3",xlab="PC 1
(40,7%)",ylab="PC3 (14,4%)")
abline(h=0)
abline(v=0)
# PC 2 VS PC 3
biplot(G23,H23,cex=0.8,main="PC2 vs PC3",ylab="PC 3
(14,4%)",xlab="PC2 (22,5%)")
abline(h=0)
abline(v=0)
# BIPLOT 3D
G3 <- G[,1:3]
H3 <- H[,1:3]
scaling_factor <- 2
H3_scaled <- H3 * scaling_factor

# VISUALISASI 3D
# KELOMPOK KECAMATAN
G3_C <- as.data.frame(G3)
G3_C[,4] <- NA
G3_C[,4] <- ifelse((G3_C[,1]>0 & G3_C[,2]>0 & G3_C[,3]>0),1,
 ifelse((G3_C[,1]>0 & G3_C[,2]<0&G3_C[,3]>0),2,
 ifelse((G3_C[,1]>0& G3_C[,2]>0&G3_C[,3]<0),3,
 ifelse((G3_C[,1]>0 & G3_C[,2]<0 & G3_C[,3]<0),4,
 ifelse((G3_C[,1]<0 & G3_C[,2]>0 & G3_C[,3]>0),5,
 ifelse((G3_C[,1]<0 & G3_C[,2]>0 & G3_C[,3]<0),6,
 ifelse((G3_C[,1]<0 & G3_C[,2]<0 & G3_C[,3]>0),7,
 ifelse((G3_C[,1]<0 & G3_C[,2]<0 & G3_C[,3]<0),8,
 0))))))))
rownames(G3_C) <- data[,1]
G3_C[order(G3_C$V4),]
objek <- paste0("O", 1:75)
variabel <- paste0("X", 1:9)
plot3d(H3)
points3d(G3,color=G3_C[,4],size=5)
text3d(G3+0.09, texts=objek,col=G3_C[,4],size=5)
points3d(H3_scaled,color='blue', size=5)
text3d(H3_scaled+0.09, texts=variabel,col="blue")
coords <- NULL
for (i in 1:nrow(H3_scaled)) {
 coords <- rbind(coords,rbind(c(0,0,0),H3_scaled[i,1:3]))

}
lines3d(coords,col="blue",lwd=2)
grid3d(c("x", "y", "z"))
abclines3d(0,0,0,a=diag(3),col="black",lwd=1) 

round(rho_matrix,3)

##Keragaman Variabel
panjangvektor <- function(x) {
 p <- nrow(x)
 j <- matrix(,nrow=p,ncol=1)
 for (i in 1:p) {
 j[i,] <- sqrt(t(x[i,])%*%x[i,])
 }
 print(j)
}
panjangvektor(H3)

##Korelasi Antar Variabel
cor_var <- function(x) {
 p <- nrow(x)
 y <- matrix(,nrow=p,ncol=p)
 for (i in 1:p) {
 for (j in 1:p) {
 y[i,j] <- sum(x[i,]*x[j,]) / (sqrt(sum(x[i,]^2)) *
 sqrt(sum(x[j,]^2)))
 }
 }
 print(y)
}

## Nilai Variabel Pada Suatu Objek
korelasi_var_ob <- function(X, Y) {
 p <- nrow(X)
 q <- nrow(Y)
 r <- matrix(, nrow = p, ncol = q)
 for (i in 1:p) {
 for (j in 1:q) {
 r[i, j] <- sum(X[i, ]* Y[j, ]) / (sqrt(sum(X[i, ]^2)) *
sqrt(sum(Y[j, ]^2)))
 }
 }
 print(r)
}

# INTERPRETASI MATRIKS KOSINUS
int_sudut <- function(x) {
  p <- nrow(x)
  k <- ncol(x)
  y <- matrix(NA, nrow = p, ncol = k)

  for (i in 1:p) {
    for (j in 1:k) {
      y[i, j] <- ifelse((x[i, j] >= 0.90 & x[i, j] <= 1), "Sangat erat positif",
                 ifelse((x[i, j] >= 0.70 & x[i, j] < 0.90), "Erat positif",
                 ifelse((x[i, j] >= 0.40 & x[i, j] < 0.70), "Cukup erat positif",
                 ifelse((x[i, j] >= 0.20 & x[i, j] < 0.40), "Hubungan lemah positif",
                 ifelse((x[i, j] >= 0.00 & x[i, j] < 0.20), "Tidak berhubungan",
                 ifelse((x[i, j] > -0.20 & x[i, j] < 0.00), "Tidak berhubungan negatif",
                 ifelse((x[i, j] >= -0.40 & x[i, j] <= -0.20), "Hubungan lemah negatif",
                 ifelse((x[i, j] >= -0.70 & x[i, j] < -0.40), "Cukup erat negatif",
                 ifelse((x[i, j] >= -0.90 & x[i, j] < -0.70), "Erat negatif",
                 "Sangat erat negatif")))))))))
    }
  }

  return(y)
}

# PEMBOBOTAN VARIABEL BERDASARKAN PCA
# Bobot formatif
q1 <- eigenvalues[1] / sum(eigenvalues)
q2 <- eigenvalues[2] / sum(eigenvalues)
q3 <- eigenvalues[3] / sum(eigenvalues)
q4 <- eigenvalues[4] / sum(eigenvalues)
q5 <- eigenvalues[5] / sum(eigenvalues)
q6 <- eigenvalues[6] / sum(eigenvalues)
q7 <- eigenvalues[7] / sum(eigenvalues)
q8 <- eigenvalues[8] / sum(eigenvalues)
q <- data.frame(PC = c("PC1","PC2",
                       "PC3","PC4","PC5","PC6", "PC7", "PC8"),q =
                  c(q1,q2,q3,q4,q5,q6,q7,q8))

PC1 <- q1
PC2 <- PC1+q2
PC3 <- PC2+q3
PC4 <- PC3+q4
PC5 <- PC4+q5
PC6 <- PC5+q6
PC7 <- PC6+q7
PC8 <- PC7+q8
PC <- data.frame(Cummulative =
                   c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8))                 
cbind(q,PC)

sqrtvektor = sqrt(eigenvectors^2)
round(sqrtvektor,4)
rata_vektor <- colMeans(sqrtvektor)
bobot_formatif <- sweep(sqrtvektor, 2, rata_vektor, FUN = "/") * q1
round(bobot_formatif,4)

# Bobot terstandar
jumlah_bobot <- colSums(bobot_formatif)
round(jumlah_bobot,4)
bobot_terstandar <- sweep(bobot_formatif, 2, jumlah_bobot, FUN = "/")
round(bobot_terstandar,4)

# Bobot akhir
bobot_akhir <- sweep(bobot_terstandar, 2, q$q, FUN = "*")
print(bobot_akhir)
round(bobot_akhir,4)
round(colSums(bobot_akhir),4) #proporsi varians
bobot_variabel = round(rowSums(bobot_akhir),4) ; bobot_variabel






