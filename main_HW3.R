
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
dir.create("images_bootstrap", showWarnings = TRUE)

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

layouts <- layout.circle(graph_asd) # get general layout

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


S <- 500

difference_matrix <- asd_bind_cor - td_bind_cor
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

png(filename='images_bootstrap/Difference_distribution.png', width = 1000, height = 600) # open image
par(mfrow = c(2, 3))
for (i in 100:105){
  if (i!=105){
  y <- array(boots_tensor_res[1,i+1,])
  hist(y, freq=FALSE, col = 'CornflowerBlue', ylim = c(0,max(density(y)$y)+1),
       main = bquote(paste('Distribution of ROI'[1],' and ROI',''[.(i+1)] )),
       font.main=4,
       xlab = 'Correlation Values',
       cex.lab=1.7,
       cex.axis=1.5,
       cex.main=2)
  curve(dnorm(x,mean=mean(y),sd=sd(y)), add=TRUE,col="red", lwd=3)
  lines(density(y),col='green', lwd=3, lty=3)
  }
  else{
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("center", legend=c("Normal Distribution", "Actual Distribution"),
           col=c("red", "green"), lty=2:3, pt.cex=5, cex=2.5, bty='o')
  }
}
dev.off() # close and save image

#This function calculates the Percentile CI interval
percentile_CI <- function(datax, alpha,point_estimate){
  # here the point_estimate has been inserted in order to avoid an if statement
  # in the diff_matrix function
  return(quantile(datax, c(alpha/2, 1-alpha/2)))
}


#This function calculates the Normal CI interval
normal_CI <-function(datax,alpha,point_estimate){
  z <- qnorm(1-alpha/2)
  rho <- point_estimate
  se <- sqrt(var(datax))
  intervals <- c(rho - z*se, rho + z*se)
  return(intervals)
}


TEST_MATRIX_percentile=asd_cor[[1]] # here we call asd_cor just to copy the correct names on each row and col.

diff_matrix <- function(interval.function){
  TEST_MATRIX=asd_cor[[1]] # here we call asd_cor just to copy the correct names on each row and col.
  for (a in 1:116) {
    for (b in a:116){
      if (a!=b){
        intervals=interval.function(array(boots_tensor_res[a,b,],dim=S),
                                    .05/choose(116,2),
                                    point_estimate=difference_matrix[a,b]) # compute intervals
        if (0>= intervals[2] | 0<= intervals[1]){
          #if zero is not in the interval put 1 in the matrix
          TEST_MATRIX[a,b]=1
          TEST_MATRIX[b,a]=1 
        }
        else {#if zero is in the interval put 0 in the matrix
          TEST_MATRIX[a,b]=0
          TEST_MATRIX[b,a]=0
        }
      }
      else TEST_MATRIX[a,b]=0 # put zero on the diagonal
      
    }
  }
  return(TEST_MATRIX) #return edge matrix
}


# compute difference graph using Normal CI
difference_graph_percentile<- graph_from_adjacency_matrix(diff_matrix(percentile_CI),
                                                          mode = c("undirected"))
# compute difference graph using Percentile CI
difference_graph_normal<- graph_from_adjacency_matrix(diff_matrix(normal_CI),
                                                      mode = c("undirected"))
# save Normal graph
f_plot(difference_graph_percentile,name='images_bootstrap/difference_graph_percentile.png')
# save Percentile graph
f_plot(difference_graph_normal,name='images_bootstrap/difference_graph_normal.png')




sum(diff_matrix(normal_CI))/2
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

