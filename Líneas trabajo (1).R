# Primero emigración y después inmigración
# 
# Elementos
# d: vector de capacidades de dispersión de las especies. Largo S
# M: Matriz de migración, con comunidades en la columnas y especies en las filas
# E: Matriz de emigración
# Hay una poisson
# Mmetacom : matriz d ela metacomunidad (especie por sitio)
# rpois(n = 1, lambda = ver)
# 
# M: matriz metacomunitaria (S x C)
# d: vector de capacidades de dispersion (S). Emtre 0 y 1.
# Md: matriz de dispersión (C x C). Depende del paisaje (acá meter variaciones)
# Mf: matriz de filtro local (S x C). Representa la capacidad de establecimiento. Valores entre 0 y 1.
# Epot: matriz de emigración potencial
# 
# en la dinámica, M va variando: Mt, Mt+1, etc
#sample(x = 1:200, 50)
M <- matrix(50, 5, 10) # todas con 50 individuos, 5 ssp x 10 comunidades
M  
d <- rep(0.1, 5)
d <- seq(from = 0.1, to = 0.9, length.out = 5)
d
Md <- matrix(sample(c(0, 1), 100, replace = TRUE), ncol = ncol(M), nrow = ncol(M)) 
Md ## es paisaje
## si pongo un pool regional le sumo el mismo valor a cada elemanto del vector
#### 
Mmigra <- Md/(apply(Md, 1, sum)) ### Md estandarizada
Mmigra
rowSums(Mmigra)
Epot <- M*d 
Epot
#Epot<-apply(Epot, 2, function(x)x/sum(x)) # NO ESTANDARIZAR!!

## o puede ser sin un for
# set.seed(3)
#  vect <- NULL
#  for (i in Epot)  {
#   vect <- c(vect, rpois(n = 1, lambda = i))
# }
#  vect
#  

#set.seed(3)
vect <-  rpois(n = length(Epot), lambda = Epot)
vect
Efect <- matrix(vect, ncol = ncol(Epot), nrow = nrow(Epot), byrow = FALSE)
Efect

#Efect<-apply(Efect, 2, function(x)x/sum(x))

## mutriplicamos por la dispersión
# Mfpot <- Efect %*% Md
# Mfpot

Mfpot <- (Efect %*% Mmigra)
Mfpot ### Mato llamó Dis
sum(Efect) == sum(Mfpot) ### tienen que dar lo mismo

Filtro <- matrix(1, nrow = nrow(M), ncol = ncol(M)) #### valores 0 filtra mucho, valores 1 pasa todo
Filtro <- matrix(sample(seq(from = 0.01, to = 1, length.out = length(M)), replace = T),
                 nrow = nrow(M), ncol = ncol(M))
# Filtro <- matrix(abs(rnorm(n = length(M), mean = 0, sd = 0.5)), 
#                  nrow = nrow(M), ncol = ncol(M))



Ipot <- Mfpot * Filtro
Ipot == Mfpot ## cuando Filtro son 1

vect1 <-  rpois(n = length(Ipot), lambda = Ipot)
vect1
Inmigra <- matrix(vect1, ncol = ncol(Ipot), nrow = nrow(Ipot), byrow = FALSE)
Inmigra

Emigra <- Efect
Emigra
Inmigra

###### PARTE 2 ####
# r es el vector de tasas de crecimiento per cápita, largo S, mismo que d
r <- rep(0.01, nrow(M))
r
M
M1 = M - Emigra + Inmigra + exp(M*r)
M; sum(M)
M1; sum(M1)

### nicho y filtro
### matrix de nicho (2 columnas: media y desvío, filas son las sopp
desvio <- 0.3
medias.spp <- sample(seq(from = desvio, to = 1- desvio, length.out = nrow(M)), replace = TRUE)
medias.spp 
desvio.spp <- rep(desvio, nrow(M))
desvio.spp

medias.com <- sample(seq(from = desvio, to = 1- desvio, length.out = ncol(M)), replace = TRUE) ### o con runif
medias.com
desvio.com <- rep(desvio, ncol(M))
desvio.com

medias <- c(medias.spp, medias.com)
desvios <- c(desvio.spp, desvio.com)
medias;desvios

#### funcion overlap (acá usamos ov.matrix_vm, que es una modificacióin de Mari y Vero a la función de Mato)
Ov <- ov.matrix_vm(medias, desvios, M)
# out  
# Ov.spp <- out[1:nrow(M), 1:nrow(M)]
# Ov.comm.spp <- out[1:nrow(M), (nrow(M) + 1): ncol(out)]
# OUT <- list(Ov = out, Ov.spp = Ov.spp, Ov.comm.spp = Ov.comm.spp)
# OUT
Ov.spp <- Ov$Ov.spp # 5 x 5 especies
Ov.comm.spp <- Ov$Ov.comm.spp ### 5 especies por 10 comunidades

#### generamos Kmax
Kmax <- runif(ncol(M), min = 500, max = 1000) ### Kmax es el K de cada comunidad, largo comm
K <- t(apply(Ov.comm.spp, 1, function(x) x * Kmax)) ### capacidades de carga por especie
### ver si no hay que dividir entre el totol de cada comunidad, porque sumando los ind de las spp
## de cada sitio supero la capacidad de carga
# Ov.comm.spp[1,1] * Kmax[1]
# Ov.comm.spp[2,1] * Kmax[1]
# Ov.comm.spp[3,1] * Kmax[1]

#### parte 3 modelo logístico  ### 
### los betas
### 
# los alfas sj están en Ov.spp, los Ks están en K
b <- function(comunidad) {
  out <- matrix(NA, nrow = length(comunidad), ncol = length(comunidad))
  for (i in 1: length(comunidad))
    {a <- comunidad/comunidad[i]
    out[i, ] <- a
  }
  diag(out) <- 0
  out
}

### prueba
comunidad <- c(5, 6, 7, 8)
b(comunidad)
comunidad[1] / comunidad[1]
comunidad[1] / comunidad[2]
comunidad[1] / comunidad[3]
comunidad[2] / comunidad[1]

b(K[,1]) ## ok
apply(K, 2, b) ## no funka, junta todo en una matriz
listaKdivs <- lapply(as.data.frame(K), b) ### esto da la lista de matrices Kij/ks para cada parche (comm)
## multiplico por alfas
diag(Ov.spp) <- 0
listaKdivs[[1]] * Ov.spp ## acá ojo que esta multiplicación está bien solo porque la matriz de overlaps es simétricp
listabetas <- lapply(listaKdivs, function(x) x * Ov.spp)
listabetas[[1]] ## ok

listabetas[[2]] == listaKdivs[[2]] * Ov.spp
## para el parche 1 M1[,1]
competenciainter <- sum(listabetas[[1]][1, ] * M1[,1], na.rm = TRUE) ## esto es lo que pierde por competencia cada especie dentro de la comunidad 1
## esto habría que repetirlo para cada especies

## esto en realidad es una multiplicación matricial, como sigue:
M1[,1] %*% t(listabetas[[1]] )
Meta.compet <- matrix(0, ncol = ncol(M), nrow = nrow(M))
for(i in 1:ncol(M)) {
  Meta.compet[,i] <-  M1[,i] %*% t(listabetas[[i]] )
}
Meta.compet

### ahora hacemos la matriz de competencia 
competencias <- matrix(0, ncol = ncol(M), nrow = nrow(M))

## logístico







