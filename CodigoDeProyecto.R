##### VICTOR ANDRES AMAYA CARVAJAL
##### Proyecto de TDA.

#Librarias que vamos utilizar:
require(stats)
library("TDAmapper")
require(fastcluster) 
library(igraph)
library(dplyr)
library(distances)


#Función con la cual vamos a reescalar cada una de las
#columnas de los datos
reescaling <- function(MAT){
  mat = MAT
  for( i in 1:ncol(mat) ) {
    valor <- (mat[,i] / (sd(mat[,i])) )
    mat[,i] <- valor
  }
  return(mat)
}


#leyendo los datos.
#dir <- "/home/vikan/Dropbox/Victor-Victor/Proyecto de verano/codigo"
dir  <-  "C://Users//victo//Dropbox//Victor-Victor//Proyecto de verano//codigo"
setwd(dir)
data <- read.table("outfile.txt", header = F)
variables <- c("ID","dist(2)","ndist(3)","bdist(4)","fflux(5)","fminor(6)","fmajor(7)",
            "fposang(8)","nflux(9)","nminor(10)","nmajor(11)","nposang(12)","nu(13)",
            "ng(14)","nr(15)","ni(16)","nz(17)","bu(18)","bg(19)","br(20)","bi(21)",
            "bz(22)","z(23)")

#Le ponemos nombre a las columnas.
colnames(data) <- variables

#Toamos nada mas las variables en las cuales estamos interesados.
subdata <- data[,c("ID","fflux(5)","nflux(9)","nu(13)","ng(14)","nr(15)",
                   "nz(17)","z(23)")]


#ahora vamos a eliminar los datos incompletos. (los que tienen -99.00)
subdata <- subdata[subdata[,"fflux(5)"] != -99.00,]
subdata <- subdata[subdata[,"nflux(9)"] != -99.00,]
subdata <- subdata[subdata[,"nu(13)"] != -99.00,]
subdata <- subdata[subdata[,"ng(14)"] != -99.00,]
subdata <- subdata[subdata[,"nr(15)"] != -99.00,]
#subdata <- subdata[subdata[,"ni(16)"] != -99.00,]
subdata <- subdata[subdata[,"nz(17)"] != -99.00,]

#cambiado la escala de las magnitudes de los filtros u,g,r,i,z.
subdata[,c("nu(13)","ng(14)","nr(15)","nz(17)")] <- 
    10^( (-0.4)*(subdata[,c("nu(13)","ng(14)","nr(15)","nz(17)")]) )

#tamano de los datos.
N = 1000

#Aquí hacemos un copia de los datos, por si hacemos algo más con la base de datos
#que vamso a trabajar.
subdata2 <- subdata[sample(1:nrow(subdata),N,replace=FALSE),]
subdata3 = subdata2

subdata3[,2:ncol(subdata2)] <- reescaling(subdata2[,2:ncol(subdata2)])
rownames(subdata3) <- as.character(subdata3[,"ID"])
subdata3 <- subdata3[,-1]


###################################################################################
# flujos como funcion de filtro.#

#vamos como se distribuyen los valores de brillo en radio.
hist(log(subdata3[,"fflux(5)"]) )
hist(log(subdata3[,"nflux(9)"]) )

#aquí vamos a convertir los datos a otro formato de modo que no perdamos la información
#del nombre de cada una de las galaxias.
P <- log(subdata3[,"fflux(5)"])
P <- matrix(P, ncol = 1)
rownames(P) <- rownames(subdata3)
P <- sapply(seq_len(nrow(P)), function(i) P[i,])

# aplicamos Mapper.
m1 <- mapper1D(
  distance_matrix = distances(subdata3),
  filter_values = P,
  num_intervals = 10,
  percent_overlap = 40,
  num_bins_when_clustering = 15)

# Graficamos con igraph
g1 <- graph.adjacency(m1$adjacency, mode="undirected")
plot( g1, layout = layout.auto(g1) ,main=paste("Datos Reescalados; filtro Fflux; N=", N))


#Ahora, procederemos a guardar el nombre de las galaxias que hay en cada vertice
#en archivos .txt para poder acceder a ellos mas facilmente.
vertices <- list()
N <- m1$num_vertices

for(i in 1:N){
  vertices[[i]]<-(as.numeric(names(m1$points_in_vertex[[i]]))) 
  write.table(vertices[i], paste(i, "vertice.txt"))
}



###################################################################################
############################    MAPPER   2D

P1 <- log(subdata3[,"nflux(9)"])
P1 <- matrix(P1, ncol = 1)
rownames(P1) <- rownames(subdata3)
P1 <- sapply(seq_len(nrow(P1)), function(i) P1[i,])


m2 <- mapper2D(
  distance_matrix = distances(subdata3),
  filter_values = list(P , P1 ),
  num_intervals = c(5,5),
  percent_overlap = 10,
  num_bins_when_clustering = 10)

g2 <- graph.adjacency(m2$adjacency, mode="undirected")
plot(g2, layout = layout.auto(g2) ,main=paste("Datos Reescalados; filtro=(Fflux,Nflux); N=", N))


#Vamos a guardar
vertices1 <- list()
N <- m2$num_vertices

for(i in 1:N){
  vertices1[[i]]<-(as.numeric(names(m2$points_in_vertex[[i]]))) 
  write.table(vertices1[i], paste(i, "vertice.txt"))
}



####################################################
library("TDA")

DiagAlphaComplex <- alphaComplexDiag(X = subdata3, printProgress = TRUE)
plot(DiagAlphaComplex[["diagram"]])


#probar con complejos testigo.
Diag <- ripsDiag(X = subdata3, 2, 1,
                 library = "GUDHI", printProgress = TRUE)
#par(mfrow = c(1, 2), mai=c(0.8, 0.8, 0.3, 0.3))
#plot(X, pch = 16, xlab = "",ylab = "")
plot(Diag[["diagram"]])


DiagLim <- 5
maxdimension <- 1
## rips diagram
Diag <- ripsDiag(X, maxdimension, DiagLim, printProgress = TRUE)
#plot
par(mfrow = c(1, 3))

plot(Diag[["diagram"]])

plot(Diag[["diagram"]], rotated = TRUE)

plot(Diag[["diagram"]], barcode = TRUE)





####################################################################3

#######################################       MAPPER 1D   PCA
principal <- prcomp(subdata3,scale. = F, center = F)
#[,2:ncol(subdata3)]
P <- principal$x[,1]

plot(P)
max(P)
min(P)

P[P == max(P)] = 0
P[P == min(P)] = 0


# Apply mapper
m1 <- mapper1D(
  distance_matrix = distances(subdata3),
  filter_values = P,
  num_intervals = 30,
  percent_overlap = 40,
  num_bins_when_clustering = 15)

# create and plot mapper graph
g1 <- graph.adjacency(m1$adjacency, mode="undirected")
plot( g1, layout = layout.auto(g1) ,main=paste("Datos Reescalados; filtro PCA1; N=", N))

#, vertex.size=c(rep(10,10))


#hola <- unique(c(as.numeric(names(m1$points_in_vertex[[15]])),
#                 as.numeric(names(m1$points_in_vertex[[19]]))))

#aqui vamos a sacar los nombres de las galaxias en cada cluster (de forma manual)
cluster1 <- c(as.numeric(names(m1$points_in_vertex[[1]])))

cluster2 <- c(as.numeric(names(m1$points_in_vertex[[2]])),
              as.numeric(names(m1$points_in_vertex[[3]])))

cluster3 <- c(as.numeric(names(m1$points_in_vertex[[11]])),
              as.numeric(names(m1$points_in_vertex[[13]])))

cluster4 <- c(as.numeric(names(m1$points_in_vertex[[4]])),
              as.numeric(names(m1$points_in_vertex[[5]])),
              as.numeric(names(m1$points_in_vertex[[6]])))

cluster5 <- c(as.numeric(names(m1$points_in_vertex[[9]])),
              as.numeric(names(m1$points_in_vertex[[7]])))

cluster6 <- c(as.numeric(names(m1$points_in_vertex[[10]])),
              as.numeric(names(m1$points_in_vertex[[8]])),
              as.numeric(names(m1$points_in_vertex[[12]])))

cluster1 <- unique(cluster1)
cluster2 <- unique(cluster2)
cluster3 <- unique(cluster3)
cluster4 <- unique(cluster4)
cluster5 <- unique(cluster5)
cluster6 <- unique(cluster6)




####################################################
library("TDA")

DiagAlphaComplex <- alphaComplexDiag(X = subdata3, printProgress = TRUE)
plot(DiagAlphaComplex[["diagram"]])

#probar con complejos testigo.

Diag <- ripsDiag(X = subdata3, 2, 1,
                 library = "GUDHI", printProgress = TRUE)
#par(mfrow = c(1, 2), mai=c(0.8, 0.8, 0.3, 0.3))
#plot(X, pch = 16, xlab = "",ylab = "")
plot(Diag[["diagram"]])


#XX1 <- circleUnif(30)
#XX2 <- circleUnif(30, r = 2) + 3
#XX <- rbind(XX1, XX2)
DiagLim <- 5
maxdimension <- 1
## rips diagram
Diag <- ripsDiag(X, maxdimension, DiagLim, printProgress = TRUE)
#plot
par(mfrow = c(1, 3))

jpeg("Diagram3000Datos")
plot(Diag[["diagram"]])
dev.off()

plot(Diag[["diagram"]], rotated = TRUE)

jpeg("Barcode3000Dato")
plot(Diag[["diagram"]], barcode = TRUE)
dev.off()


#x <- matrix(rnorm(100*3), ncol = 3)
#stopifnot(mahalanobis(x, 0, diag(ncol(x))) == rowSums(x*x))
##- Here, D^2 = usual squared Euclidean distances

#Sx <- cov(x)
#D2 <- mahalanobis(x, colMeans(x), Sx)
#plot(density(D2, bw = 0.5),
#     main="Squared Mahalanobis distances, n=100, p=3") ; rug(D2)
#qqplot(qchisq(ppoints(100), df = 3), D2,
#       main = expression("Q-Q plot of Mahalanobis" * ~D^2 *
#                           " vs. quantiles of" * ~ chi[3]^2))
#abline(0, 1, col = 'gray')
