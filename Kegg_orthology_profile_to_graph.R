##Author:Niraj kumar khemka ##########
##Contact:nirajkhe@gmail.com #########


a <-"Kegg.profile"  #enter the profile name
n <- readline(prompt="Enter an integer 1 : for Jaccard distance , 2: for Dollo distance   ")
  n <- as.integer(n)

library(igraph)  
  
## profile read  
profile_main <- read.delim(a) # read data
b <- nrow(profile_main)           # Enter the number of orthologs
rownames(profile_main) <- profile_main[,1]                  # rownames
profile_main <- profile_main[,-1]                           # removing rownames from matrix

m <- profile_main
m <- t(m)                                                 # transpose for make orthologs vs ortholgs matrix not required for organism matrix
c1<-combn(seq_len(ncol(m)),2)
mat1<- matrix(0,ncol=b,nrow=b,dimnames=list(colnames(m),colnames(m)))

# For 1 1  pair presence
vec1<-unlist(lapply(seq_len(ncol(c1)),function(i) {m1<-m[,c1[,i]];  length(which(m1[,1]!=0 & m1[,2]!=0)) }))
mat1[lower.tri(mat1)] <- vec1
mat1 <- mat1 + t(mat1)     
pair_1_1 <- mat1

# For  1 0  pair presence
mat1<- matrix(0,ncol=b,nrow=b,dimnames=list(colnames(m),colnames(m)))
vec1<-unlist(lapply(seq_len(ncol(c1)),function(i) {m1<-m[,c1[,i]];  length(which(m1[,1]!=0 & m1[,2]==0)) }))
mat1[lower.tri(mat1)]<-vec1
mat1 <- mat1 + t(mat1)
pair_1_0 <- mat1

# For  0  1  pair presence
mat1<- matrix(0,ncol=b,nrow=b,dimnames=list(colnames(m),colnames(m)))
vec1<-unlist(lapply(seq_len(ncol(c1)),function(i) {m1<-m[,c1[,i]]; length(which(m1[,1]==0 & m1[,2]!=0)) }))
mat1[lower.tri(mat1)]<-vec1
mat1 <- mat1 + t(mat1)
pair_0_1 <- mat1


# For  0  0  pair presence
mat1<- matrix(0,ncol=b,nrow=b,dimnames=list(colnames(m),colnames(m)))
vec1<-unlist(lapply(seq_len(ncol(c1)),function(i) {m1<-m[,c1[,i]]; length(which(m1[,1]==0 & m1[,2]==0)) }))
mat1[lower.tri(mat1)]<-vec1
mat1 <- mat1 + t(mat1)
pair_0_0 <- mat1


if(n==2){

dollo_dist <- log (((pair_1_1 + pair_0_1) * ( pair_1_1 + pair_1_0))/ ((pair_1_1)^2) )
diag(dollo_dist) <- 0
dollo_dist <- round(dollo_dist,3)
dollo <- dollo_dist/max(dollo_dist)  # normalize the data
dollo_sim <- 1 - dollo
write.table(dollo_sim,file="dollo_sim.txt",sep = "\t")  # write matrix into data file
g <- graph.adjacency(as.matrix(dollo_sim),weighted= TRUE,mode="undirected") 
V(g)$label <- V(g)$name                                 # label node of graph
E(g)$label <- E(g)$weight  
write.graph(g,file="dollo_sim.gml",format="gml")
}
else(n==1){
# for frequency calculation of pair_1_1
jaccard <- pair_1_1 / (pair_1_1 + pair_1_0  + pair_0_1)
diag(jaccard) <- 0
jaccard <- round(jaccard,3)
write.table(jaccard,file="jaccard.txt",sep = "\t")  # write matrix into data file

g <- graph.adjacency(as.matrix(jaccard),weighted= TRUE,mode="undirected") 
V(g)$label <- V(g)$name                                 # label node of graph
E(g)$label <- E(g)$weight  
write.graph(g,file="jaccard.gml"),format="gml")
}
