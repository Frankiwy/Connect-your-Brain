
# Libraries, Imports & Functions -------------------------------------------
set.seed(42)
library(igraph)
load("../Connect-your-Brain/data/hw3_data.RData") # import data
dir.create("images_asd", showWarnings = TRUE)
dir.create("images_td", showWarnings = TRUE)


zfisher <- function(x) 1/2*log((1+x)/(1-x)) # Z-Fisher Transform function

CI_function <- function(N,rho,alpha,t,bonferroni=TRUE){
  se = 1/sqrt(N-3)
  
  if (bonferroni){
    lower = rho - qnorm(1-(alpha/2)/choose(116,2))*se
    upper = rho + qnorm(1-(alpha/2)/choose(116,2))*se
  }  
  else{
    lower = rho - qnorm(1-(alpha/2))*se
    upper = rho + qnorm(1-(alpha/2))*se
  }
  #print(c(rho - qnorm(1-(alpha/2))*se,rho - qnorm(1-(alpha/2)/choose(116,2))*se))
  if (-t>= upper | t<= lower) return(1)
  else return(0)
}

find_edges <- function(matrix,N,t,bonferroni=TRUE){
  for (i in 1:nrow(matrix)){
    for (j in 1:ncol(matrix)){
      if (i != j) matrix[i,j] = CI_function(N, matrix[i,j], 0.05,t,bonferroni)
      else matrix[i,j] = 0
    }
  }
  return(matrix)
}

f_plot<- function(g,name=NA,width_edges=NA){
  if (is.na(width_edges)) width_edges=(0.1+E(g)$weight)^1.5*8
  
  if (!is.na(name)){
    png(filename=name, width = 10000, height = 4720)
    plot(g, vertex.size = 1, edge.size = 0.001,
         edge.width= width_edges,
         vertex.frame.color=NA, vertex.color ='red',
         edge.color='black', layout=layouts)
    dev.off()
  }
  
  else {
    plot(g, vertex.size = 1, edge.size = 0.001,
         edge.width= width_edges,
         vertex.frame.color=NA, vertex.color ='red',
         edge.color='black', layout=layouts)
  }
}

get_graphs <- function(list_df,t,bonferroni=TRUE){
  graph_list <- list()
  n = 1
  for (df in list_df){
    df_zfisher <- zfisher(df)
    edges <- find_edges(df_zfisher,nrow(df_zfisher),t,bonferroni)
    graph_get <- graph_from_adjacency_matrix(edges, mode = c("undirected") )
    graph_list[[n]] <- graph_get
    n = n+1
  }
  return(graph_list)
}

get_aggregate_graph <- function(list_df,t,bonferroni=TRUE){
  graph_matrix_list <- list()
  n = 1
  for (df in list_df){
    df_zfisher <- zfisher(df)
    edges <- find_edges(df_zfisher,nrow(df_zfisher),t,bonferroni)
    graph_matrix_list[[n]] <- edges
    n = n+1
    #print(sum(edges))
  }
  #print(graph_matrix_list)
  graph_matrix <- Reduce('+', graph_matrix_list)/length(graph_matrix_list)
  return(graph_matrix)
}

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

# Point 2 -----------------------------------------------------------------

# APPROACH 1:

asd_bind = do.call(rbind, asd_sel) # aggregate asd data
td_bind = do.call(rbind, td_sel) # aggregate td data

asd_bind_cor = cor(asd_bind) # compute asd cor
td_bind_cor = cor(td_bind) # compute td cor

asd_td_bind_matrix = rbind(asd_bind_cor, td_bind_cor) # aggregate asd & td cor

t_bind <- quantile(abs(asd_td_bind_matrix), probs = c(0.8)) # get threshold
t_bind_zfisher = zfisher(t_bind) # apply zfisher on threshold


asd_bind_matrix_zfisher_cor = zfisher(asd_bind_cor) # apply zfisher on asd
asd_bind_matrix_zfisher <- find_edges(asd_bind_matrix_zfisher_cor, 145*12, t_bind_zfisher )
graph_asd <- graph_from_adjacency_matrix(asd_bind_matrix_zfisher, mode = c("undirected") )

layouts <- layout.davidson.harel(graph_asd) # get general layout

f_plot(graph_asd,name='images_asd/asd_1.png',width_edges = 1)

td_bind_matrix_zfisher_cor = zfisher(td_bind_cor) # apply zfisher on td
td_bind_matrix_zfisher <- find_edges(td_bind_matrix_zfisher_cor, 145*12, t_bind_zfisher)
graph_td <- graph_from_adjacency_matrix(td_bind_matrix_zfisher, mode = c("undirected") )

f_plot(graph_td,name='images_td/td_1.png',width_edges = 1)


c(length(V(graph_asd)),length(V(graph_td))) # count verteces
c(length(E(graph_asd)),length(E(graph_td))) # counte edges
plot(degree_distribution(graph_asd, cumulative = TRUE), main = 'asd', ylab='density asd')
plot(degree_distribution(graph_td, cumulative = TRUE), main='td', ylab='density td')

asd_eigen_centrality = eigen_centrality(graph_asd)
plot(sort(asd_eigen_centrality$vector), main='asd', ylab = 'asd centrality')
#asd_betweenness <- betweenness(graph_asd)
#plot(sort(asd_betweenness), main='asd', ylab='asd betweenness')
td_eigen_centrality = eigen_centrality(graph_td)
plot(sort(td_eigen_centrality$vector), main='td', ylab='td centrality')
#td_betweenness <- betweenness(graph_td)
#plot(sort(td_betweenness), ylim = c(0,50))

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
t_zfisher = zfisher(t)

asd_fisher_matrix <- lapply(asd_cor, zfisher) # apply Z-Fisher Transform to all cor
td_fisher_matrix <- lapply(td_cor, zfisher) # pply Z-Fisher Transoform to all cor


edges_asd_2 <- get_aggregate_graph(asd_cor, t_bind_zfisher)
graph_asd_2 <- graph_from_adjacency_matrix(edges_asd_2,
                                           mode = c("undirected"),
                                           weighted = TRUE)

f_plot(graph_asd_2,name='images_asd/asd_2.png')

edges_td_2 <- get_aggregate_graph(td_cor, t_bind_zfisher)
graph_td_2 <- graph_from_adjacency_matrix(edges_td_2,
                                           mode = c("undirected"),
                                           weighted = TRUE)

f_plot(graph_td_2,name='images_td/td_2.png')


graph_asd_list <- get_graphs(asd_cor, t_bind_zfisher) #graph per each asd person
graph_td_list <- get_graphs(td_cor, t_bind_zfisher) #graph per each td person


m=1
for (g in graph_asd_list){
  f_plot(g,name=paste("images_asd/asd_2_person_",m,".png",sep=""),width_edges = 1)
  m=m+1 
}

m=1
for (g in graph_td_list){
  f_plot(g,name=paste("images_td/td_2_person_",m,".png",sep=""),width_edges = 1)
  m=m+1 
}



# part 1 & 2 without Bonferroni Correction --------------------------------

asd_bind_matrix_zfisher_no_bonferroni <- find_edges(asd_bind_matrix_zfisher_cor, 145*12, t_bind_zfisher,bonferroni=FALSE )
graph_asd_no_bonferroni <- graph_from_adjacency_matrix(asd_bind_matrix_zfisher_no_bonferroni, mode = c("undirected") )
f_plot(graph_asd_no_bonferroni,name='images_asd/asd_1_no_bonferroni.png',width_edges = 1)


td_bind_matrix_zfisher_no_bonferroni <- find_edges(td_bind_matrix_zfisher_cor, 145*12, t_bind_zfisher,bonferroni = FALSE)
graph_td_no_bonferroni <- graph_from_adjacency_matrix(td_bind_matrix_zfisher_no_bonferroni, mode = c("undirected") )
f_plot(graph_td_no_bonferroni,name='images_td/td_1_no_bonferroni.png',width_edges = 1)

sum(asd_bind_matrix_zfisher - asd_bind_matrix_zfisher_no_bonferroni)
sum(td_bind_matrix_zfisher - td_bind_matrix_zfisher_no_bonferroni)

######### APPROACH 2 #############


edges_asd_2_no_bonferroni <- get_aggregate_graph(asd_cor, t_bind_zfisher, bonferroni = FALSE)
graph_asd_2_no_bonferroni <- graph_from_adjacency_matrix(edges_asd_2_no_bonferroni,
                                           mode = c("undirected"),
                                           weighted = TRUE)

f_plot(graph_asd_2_no_bonferroni,name='images_asd/asd_2_new_no_bonferroni.png')

edges_asd_2 <- get_aggregate_graph(asd_cor, t_bind_zfisher)
graph_asd_2 <- graph_from_adjacency_matrix(edges_asd_2,
                                           mode = c("undirected"),
                                           weighted = TRUE)

f_plot(graph_asd_2,name='images_asd/asd_new_2.png')






edges_td_2_no_bonferroni <- get_aggregate_graph(td_cor, t_bind_zfisher, bonferroni=FALSE)
graph_td_2_no_bonferroni <- graph_from_adjacency_matrix(edges_td_2_no_bonferroni,
                                          mode = c("undirected"),
                                          weighted = TRUE)

f_plot(graph_td_2_no_bonferroni,name='images_td/td_2_no_bonferroni.png')


graph_asd_list_no_bonferroni <- get_graphs(asd_cor, t_bind_zfisher, bonferroni = FALSE) #graph per each asd person
graph_td_list_no_bonferroni <- get_graphs(td_cor, t_bind_zfisher, bonferroni = FALSE) #graph per each td person


m=1
for (g in graph_asd_list_no_bonferroni){
  f_plot(g,name=paste("images_asd/asd_2_person_",m,"_no_bonferroni.png",sep=""),width_edges = 1)
  m=m+1 
}

m=1
for (g in graph_td_list_no_bonferroni){
  f_plot(g,name=paste("images_td/td_2_person_",m,"_no_bonferroni.png",sep=""),width_edges = 1)
  m=m+1 
}




