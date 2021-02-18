
# Libraries, Imports & Functions -------------------------------------------
set.seed(42)
library(igraph)
library(tensorr)
load("../Connect-your-Brain/data/hw3_data.RData") # import data
dir.create("images_asd", showWarnings = TRUE)
dir.create("images_asd/no_Bonferroni", showWarnings = TRUE)
dir.create("images_asd/Bonferroni", showWarnings = TRUE)
dir.create("images_td", showWarnings = TRUE)
dir.create("images_td/no_Bonferroni", showWarnings = TRUE)
dir.create("images_td/Bonferroni", showWarnings = TRUE)


zfisher <- function(x) 1/2*log((1+x)/(1-x)) # Z-Fisher Transform function

CI_function <- function(N,rho,alpha,t,bonferroni=TRUE){
  se = 1/sqrt(N-3) # standard error
  
  if (bonferroni){ # compute CI using Bonferroni correction
    lower = rho - qnorm(1-(alpha/2)/choose(116,2))*se 
    upper = rho + qnorm(1-(alpha/2)/choose(116,2))*se
  }  
  else{ # compute CI wihout Bonferroni correction
    lower = rho - qnorm(1-(alpha/2))*se
    upper = rho + qnorm(1-(alpha/2))*se
  } # check if an edge has to be added or not
  if (-t>= upper | t<= lower) return(1) 
  else return(0)
}

find_edges <- function(matrix,N,t,bonferroni=TRUE){
  for (i in 1:nrow(matrix)){
    for (j in 1:ncol(matrix)){
      # compute CI and check for edging iff we are not dealing with the 
      # diagonal matrix
      if (i != j) matrix[i,j] = CI_function(N, matrix[i,j], 0.05,t,bonferroni)
      else matrix[i,j] = 0 # add zero if diagonal matrix
    }
  }
  return(matrix)
}

f_plot<- function(g,name=NA,width_edges=NA){
  # the function is used to plot the graphs and, if a name is given, saves it
  if (is.na(width_edges)) width_edges=(0.1+E(g)$weight)^1.5*8
  
  if (!is.na(name)){
    png(filename=name, width = 10000, height = 4720) # open image
    plot(g, vertex.size = 1, edge.size = 0.001,
         edge.width= width_edges,
         vertex.frame.color=NA, vertex.color ='red',
         edge.color='black', layout=layouts)
    dev.off() # close and save image
  }
  
  else {
    plot(g, vertex.size = 1, edge.size = 0.001,
         edge.width= width_edges,
         vertex.frame.color=NA, vertex.color ='red',
         edge.color='black', layout=layouts)
  }
}

get_graphs <- function(list_df,t,bonferroni=TRUE){
  # the function takes in input a list of correlation matrices and returns
  # a list of graphs
  graph_list <- list() # list where store the graphs
  n = 1
  for (df in list_df){
    df_zfisher <- zfisher(df) # compute z fisher on all correlations
    edges <- find_edges(df_zfisher,nrow(df_zfisher),t,bonferroni) # get adjacency matrix
    # use the adjacency matrix to get the graph
    graph_get <- graph_from_adjacency_matrix(edges, mode = c("undirected") ) 
    graph_list[[n]] <- graph_get
    n = n+1
  }
  return(graph_list)
}

get_aggregate_graph <- function(list_df,t,bonferroni=TRUE){
  # the function takes in input a list of matrices and returns an adjacency matrix
  graph_matrix_list <- list() # list where store all 12 adjacency matrix
  n = 1
  for (df in list_df){
    df_zfisher <- zfisher(df) # apply z-fisher on each corr matrix
    edges <- find_edges(df_zfisher,nrow(df_zfisher),t,bonferroni)
    graph_matrix_list[[n]] <- edges # store adjacency matrix inside the list 
    n = n+1
  }
  # return a normalized weighted graph 
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

f_plot(graph_asd,name='images_asd/Bonferroni/asd_1.png',width_edges = 1)

td_bind_matrix_zfisher_cor = zfisher(td_bind_cor) # apply zfisher on td
td_bind_matrix_zfisher <- find_edges(td_bind_matrix_zfisher_cor, 145*12, t_bind_zfisher)
graph_td <- graph_from_adjacency_matrix(td_bind_matrix_zfisher, mode = c("undirected") )

f_plot(graph_td,name='images_td/Bonferroni/td_1.png',width_edges = 1)


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


edges_asd_2 <- get_aggregate_graph(asd_cor, t_zfisher)
graph_asd_2 <- graph_from_adjacency_matrix(edges_asd_2,
                                           mode = c("undirected"),
                                           weighted = TRUE)

f_plot(graph_asd_2,name='images_asd/Bonferroni/asd_2.png')

edges_td_2 <- get_aggregate_graph(td_cor, t_zfisher)
graph_td_2 <- graph_from_adjacency_matrix(edges_td_2,
                                           mode = c("undirected"),
                                           weighted = TRUE)

f_plot(graph_td_2,name='images_td/Bonferroni/td_2.png')


graph_asd_list <- get_graphs(asd_cor, t_zfisher) #graph per each asd person
graph_td_list <- get_graphs(td_cor, t_zfisher) #graph per each td person


m=1
for (g in graph_asd_list){
  f_plot(g,name=paste("images_asd/Bonferroni/asd_2_person_",m,".png",sep=""),width_edges = 1)
  m=m+1 
}

m=1
for (g in graph_td_list){
  f_plot(g,name=paste("images_td/Bonferroni/td_2_person_",m,".png",sep=""),width_edges = 1)
  m=m+1 
}



# part 1 & 2 without Bonferroni Correction --------------------------------

asd_bind_matrix_zfisher_no_bonferroni <- find_edges(asd_bind_matrix_zfisher_cor, 145*12, t_bind_zfisher,bonferroni=FALSE )
graph_asd_no_bonferroni <- graph_from_adjacency_matrix(asd_bind_matrix_zfisher_no_bonferroni, mode = c("undirected") )
f_plot(graph_asd_no_bonferroni,name='images_asd/no_Bonferroni/asd_1_no_bonferroni.png',width_edges = 1)


td_bind_matrix_zfisher_no_bonferroni <- find_edges(td_bind_matrix_zfisher_cor, 145*12, t_bind_zfisher,bonferroni = FALSE)
graph_td_no_bonferroni <- graph_from_adjacency_matrix(td_bind_matrix_zfisher_no_bonferroni, mode = c("undirected") )
f_plot(graph_td_no_bonferroni,name='images_td/no_Bonferroni/td_1_no_bonferroni.png',width_edges = 1)

sum(asd_bind_matrix_zfisher - asd_bind_matrix_zfisher_no_bonferroni)
sum(td_bind_matrix_zfisher - td_bind_matrix_zfisher_no_bonferroni)

######### APPROACH 2 #############


edges_asd_2_no_bonferroni <- get_aggregate_graph(asd_cor, t_zfisher, bonferroni = FALSE)
graph_asd_2_no_bonferroni <- graph_from_adjacency_matrix(edges_asd_2_no_bonferroni,
                                           mode = c("undirected"),
                                           weighted = TRUE)

f_plot(graph_asd_2_no_bonferroni,name='images_asd/no_Bonferroni/asd_2_no_bonferroni.png')



edges_td_2_no_bonferroni <- get_aggregate_graph(td_cor, t_zfisher, bonferroni=FALSE)
graph_td_2_no_bonferroni <- graph_from_adjacency_matrix(edges_td_2_no_bonferroni,
                                          mode = c("undirected"),
                                          weighted = TRUE)

f_plot(graph_td_2_no_bonferroni,name='images_td/no_Bonferroni/td_2_no_bonferroni.png')


graph_asd_list_no_bonferroni <- get_graphs(asd_cor, t_zfisher, bonferroni = FALSE) #graph per each asd person
graph_td_list_no_bonferroni <- get_graphs(td_cor, t_zfisher, bonferroni = FALSE) #graph per each td person


m=1
for (g in graph_asd_list_no_bonferroni){
  f_plot(g,name=paste("images_asd/no_Bonferroni/asd_2_person_",m,"_no_bonferroni.png",sep=""),width_edges = 1)
  m=m+1 
}

m=1
for (g in graph_td_list_no_bonferroni){
  f_plot(g,name=paste("images_td/no_Bonferroni/td_2_person_",m,"_no_bonferroni.png",sep=""),width_edges = 1)
  m=m+1 
}

############### BONUS ########################

# we have decided to stick with 12 people in each bootstrap, in order 
# to maintain lower the variance (ADD THIS COMMENT)


S <- 1000

boots_delta_res <- list()
dims <- c(116,116,S)
vals <- array(rep(0,116*116*S),dims)
boots_tensor_res <- dtensor(vals)

for (i in 1:S){
  
  pick_sample_asd <- sample(1:12, 12, replace=T) #pick a sample of 12 individuals from the asd population
  pick_sample_td <- sample(1:12, 12, replace=T) #pick a sample of 12 individuals from the td population
  
  #We'll now follow the first approach and put the all the data in two 1740x116 matrices.
  boots_asd <- do.call(rbind,asd_sel[pick_sample_asd])
  boots_td <- do.call(rbind,td_sel[pick_sample_td])
  
  #We calculate the correlation matrix for both the asd and the 
  boots_asd_cor <- cor(boots_asd)
  boots_td_cor <- cor(boots_td)
  boots_delta <- boots_asd_cor - boots_td_cor
  boots_tensor_res[,,i] <- boots_delta #We store the computed matrix with all the difference as one slice of a tensor.
  if (i%%100 == 0) print(i)
}

#This function calculates the CI interval
myfun_CI <- function(datax, alpha){
  return(quantile(datax, c(alpha/2, 1-alpha/2)))
}

find_p_value <-function(datax){
  data=sort(datax)
  #The error due to this is negligible for S large enugh
  data=append(data,-Inf,after=0)
  data=append(data,Inf,after=length(data))

  if (length(data)%%2==1) {
    print ("not implemented for odd lenght")
    return (0)
  }
  else {
    counter=2
    index_1=length(data)/2
    index_2=index_1+1
    while (!(data[index_1]<0 && data[index_2]>0)){
      index_1=index_1-1
      index_2=index_2+1
      }
    return (1-(index_2-index_1+1)/(length(data)-2))
    }
  }
  
  
  
find_p_value_2 <-function(datax){
  
  lower=0
  upper=1
  
  for (i in 0:32){
    CI=myfun_CI(datax, (lower+upper)/2)
    if (CI[1]<0 && CI[2]>0) lower=(upper+lower)/2
    else upper=(upper+lower)/2
  }
  
  return ((lower+upper)/2)
}

results=rep(0,6670)
TEST_MATRIX=asd_cor[[1]]


i=1

for (a in 1:116) {
  for (b in a:116){
    if (a!=b){
      results[i]=find_p_value_2(array(boots_tensor_res[a,b,],dim=S))
      TEST_MATRIX[a,b]=results[i]
      TEST_MATRIX[b,a]=results[i]
      i=i+1}
    else TEST_MATRIX[a,b]=0
    
  }
}

sorted_results=sort(results)

find_j <-function(data,alpha=0.1){
  t_bon=alpha/length(data)
  k_max=-1
  for (k in 1:length(data)){
    if (data[k]< k*t_bon) k_max=k
    
  }
  return (data[k_max])
}

t_bh=find_j(sorted_results,alpha=0.1)

for (a in 1:116){
  for (b in 1:116){
    if (a!=b){
      if (TEST_MATRIX[a,b]<=t_bh) TEST_MATRIX[a,b]=1
      
      else TEST_MATRIX[a,b]=0
    }
  }
}

difference_graph<- graph_from_adjacency_matrix(TEST_MATRIX,
                                           mode = c("undirected")
                                  )

f_plot(difference_graph,name='difference_graph.png')


##############################

gsize(graph_asd)
gsize(graph_td)
asd_sum <- 0
td_sum <- 0
for (i in 1:12) {
  asd_sum <- asd_sum + gsize(graph_asd_list[[i]])
  td_sum <- td_sum + gsize(graph_td_list[[i]])
}

graph_asd_test <- graph_asd %s% graph_td

f_plot(graph_asd_test,name='images_td/intersection.png')

as_edgelist(graph_asd)

asd_sum
td_sum





edge_attr(graph_asd_2, index =c(8102) )

