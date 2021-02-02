
# Point 1 -----------------------------------------------------------------

load("../Connect-your-Brain/data/hw3_data.RData") # import data
library(ggraph)
library(igraph)

asd_sel # autistc brain list
td_sel # tipiccaly developed brain list
zfisher <- function(x) 1/2*log((1+x)/(1-x)) # Z-Fisher Transform function



# Point 0 -----------------------------------------------------------------

# let's check t on every patient:
t_asd <- rep(NA,12)
t_td <- rep(NA,12)
for (n in 1:12){
  
  asd_cor <- cor(asd_sel[[n]])
  td_cor <- cor(td_sel[[n]])
  t_asd[n] <- quantile(abs(asd_cor), probs = c(0.8)) 
  t_td[n] <- quantile(abs(td_cor), probs = c(0.8)) 
  
}
#hist(t_asd,breaks = 10)


# Point 2 -----------------------------------------------------------------

# APPROACH 1:

asd_bind = do.call(rbind, asd_sel)
td_bind = do.call(rbind, td_sel)

asd_bind_cor = cor(asd_bind)
td_bind_cor = cor(td_bind)

asd_td_bind_matrix = rbind(asd_bind_cor, td_bind_cor) 

t_bind <- quantile(abs(asd_td_bind_matrix), probs = c(0.8)) 
t_bind_zfisher = zfisher(t_bind)
t_zfisher = zfisher(t)

CI_function <- function(N,rho,alpha,t){
  se = 1/sqrt(N-3)
  lower = rho - qnorm(1-(alpha/2)/choose(116,2))*se
  upper = rho + qnorm(1-(alpha/2)/choose(116,2))*se
  if (-t>= upper | t<= lower) return(1)
  else return(0)
}




asd_bind_matrix_zfisher = zfisher(asd_bind_cor)

find_edges <- function(matrix,N,t){
  for (i in 1:nrow(matrix)){
    for (j in 1:ncol(matrix)){
      if (i != j) matrix[i,j] = CI_function(N, matrix[i,j], 0.05,t)
      else matrix[i,j] = 0
    }
  }
  return(matrix)
}
asd_bind_matrix_zfisher = zfisher(asd_bind_cor)
asd_bind_matrix_zfisher <- find_edges(asd_bind_matrix_zfisher, 145*12, t_bind_zfisher )
graph_asd <- graph_from_adjacency_matrix(asd_bind_matrix_zfisher, mode = c("undirected") )
plot(graph_asd, vertex.size = 1, edge.size = 0.001, vertex.label=NA,
     vertex.frame.color=NA, vertex.color ='red', edge.color='black')

td_bind_matrix_zfisher = zfisher(td_bind_cor)
td_bind_matrix_zfisher <- find_edges(td_bind_matrix_zfisher, 145*12, t_bind_zfisher )
graph_td <- graph_from_adjacency_matrix(td_bind_matrix_zfisher, mode = c("undirected") )
plot(graph_td, vertex.size = 1, edge.size = 0.001, vertex.label=NA,
     vertex.frame.color=NA, vertex.color ='red', edge.color='black')

c(length(V(graph_asd)),length(V(graph_td)))
c(length(E(graph_asd)),length(E(graph_td)))
plot(degree_distribution(graph_asd, cumulative = FALSE))
plot(degree_distribution(graph_td, cumulative = FALSE))

asd_eigen_centrality = eigen_centrality(graph_asd)
plot(sort(asd_eigen_centrality$vector))
asd_betweenness <- betweenness(graph_asd)
plot(sort(asd_betweenness), ylim = c(0,50))

td_eigen_centrality = eigen_centrality(graph_td)
plot(sort(td_eigen_centrality$vector))
td_betweenness <- betweenness(graph_td)
plot(sort(td_betweenness), ylim = c(0,50))

#write_graph(graph_asd, "graph_asd.txt", format = c("edgelist"))
#write_graph(graph_td, "graph_td.txt", format = c("edgelist"))

######################### APPRACH 2 ###########################

asd_cor = lapply(asd_sel, cor) # list of correlation matrices
td_cor = lapply(td_sel, cor) # list of correlation matrices

asd_matrix = do.call(rbind, asd_cor) # combine correlation matrices into 1 matrix
td_matrix = do.call(rbind, td_cor)  # combine correaltion matrices into 1 matrix

# combine preiovusly obtained correlation matrices into 1 big corr-matrix
asd_td_matrix = rbind(asd_matrix, td_matrix) 

# get 80% percentile from the previously obtained matrix
t <- quantile(abs(asd_td_matrix), probs = c(0.8)) 

asd_fisher_matrix <- lapply(asd_cor, zfisher) # apply Z-Fisher Transform to all cor
td_fisher_matrix <- lapply(td_cor, zfisher) # pply Z-Fisher Transoform to all cor







