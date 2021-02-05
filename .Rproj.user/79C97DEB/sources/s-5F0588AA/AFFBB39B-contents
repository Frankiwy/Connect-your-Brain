library(igraph)
library(sdpt3r)
library(glue)
library(mc2d)
library(ggplot2)
library(knitr)

set.seed(161821)


ER_graph = erdos.renyi.game(10,p=1/3, type = 'gnp', directed=FALSE) # Erdős–Rényi
VS_graph = graph_from_literal(1--2, 1--4, 2--3, 4--3, 3--5, 5--6, 5--7) # the very small one

include_nodes <- c(2,3,5)
plot(VS_graph, vertex.size=25, vertex.frame.color="black",
     vertex.label.color="black",
     vertex.color=c("purple",'green')[1+(V(VS_graph)%in%include_nodes)],
     vertex.label.cex=1.5, vertex.label.dist=0, edge.curved=0, edge.color=('black'),
     xlim=c(-1,1), ylim=c(-1,1))
title("'Very' small one", adj = 0, line = -5)

list_nodes <- c(3,4,5,10)
plot(ER_graph, vertex.color=c("purple",'green')[1+(V(VS_graph)%in%list_nodes)],
     vertex.size=30, vertex.frame.color="black", vertex.label.color="black",
     vertex.label.cex=1, vertex.label.dist=0,
     edge.curved=0, edge.color='black',
     xlim=c(-1,1), ylim=c(-1,1))
title("Erdős–Rényi", adj = 0, line = -5)


ER_matrix = as.matrix(as_adj(ER_graph))
VS_matrix = as.matrix(as_adj(VS_graph))


print(glue('ER maxcut: {round(maxcut(ER_matrix)$pobj,2)} & VS maxcut: {round(maxcut(VS_matrix)$pobj,2)}'))


for (m in 1:2){
  
  g_matrices <-list(ER_matrix, VS_matrix)
  name_graphs <- c('ER_graph', 'VS_graph')
  OPT_names <- c('OPT_sdpt3r', 'OPT_E')
  g_dimension = dim(g_matrices[[m]])[1]
  
  sum = 0
  for (n in 1:100000){
    mask = rbern(n = g_dimension, prob=.5) 
    vector = c(1:g_dimension)
    indeces_U = vector[!!mask] # this is a little workaround to have a mask with TRUE and FALSE instead of 1 0
    indeces_not_U = vector[!mask] 
    num_edges = sum(g_matrices[[m]][indeces_U, indeces_not_U]) 
    sum = sum+num_edges
  }
  exp = sum/100000
  print(glue('{name_graphs[m]}: EXP= {exp} & {OPT_names[m]}/2= {round(-maxcut(g_matrices[[m]])$pobj/2,2)}'))
}


Paul_Erdos_function = function(gen_graph){
  
  starting = proc.time()
  g_matrix = as.matrix(as_adj(gen_graph))
  g_dimension = dim(g_matrix)
  
  sum = 0
  for (n in 1:10000){
    mask = rbern(n = g_dimension, prob=.5) 
    vector = c(1:g_dimension)
    indeces_U = vector[!!mask] # this is a little workaround to have a mask with TRUE and FALSE instead of 1 0
    indeces_not_U = vector[!mask] 
    num_edges = sum(g_matrix[indeces_U, indeces_not_U]) 
    sum = sum+num_edges
  }
  exp = sum/10000
  end_time = proc.time() - starting 
  return(end_time[[3]])
}

maxcut_function = function(gen_graph){
  
  starting = proc.time()
  
  g_matrix = as.matrix(as_adj(gen_graph))
  maxcut_value = maxcut(g_matrix)$pobj
  
  end_time = proc.time() - starting 
  return(end_time[[3]])
}


time_elapsed_PE <- c()

for (n in seq(from = 10, to = 200, by = 10)) {
  intermediated_time = 0
  for (l in 1:5){
    intermediated_time = intermediated_time +  
      Paul_Erdos_function(erdos.renyi.game(n,p=1/3, type = 'gnp', directed=FALSE))
  }
  time_elapsed_PE <- c(time_elapsed_PE, round(intermediated_time/5,3))
}

time_elapsed_maxcut <- c()
for (n in seq(from = 10, to = 200, by = 10)) {
  intermediated_time = 0
  for (l in 1:5){
    intermediated_time = intermediated_time +  
      maxcut_function(erdos.renyi.game(n,p=1/3, type = 'gnp', directed=FALSE))
  }
  time_elapsed_maxcut <- c(time_elapsed_maxcut, round(intermediated_time/5,3))
}

PE_performances <- data.frame('n_of_nodes' = seq(from = 10, to = 200, by = 10),
                              'elapsed_time_PE' = time_elapsed_PE,
                              'elapsed_time_maxcut' = time_elapsed_maxcut)

write.csv(PE_performances,"performances.csv", row.names = FALSE)

PE_performances <- read.csv("performances.csv")

kable(PE_performances, caption= "Time colplexity Table")


ggplot(PE_performances, aes(x=n_of_nodes))+
  
  geom_line(aes(y=time_elapsed_PE, color="#09611D"), size=1) +
  geom_point(aes(y=time_elapsed_PE), shape=21, color="#09611D", fill="#F34835", size=3) + 
  
  geom_line(aes(y=time_elapsed_maxcut, color="#F3C035"), size=1) + 
  geom_point(aes(y=time_elapsed_maxcut), color="#F3C035", shape=21, fill="purple", size=3) +
  
  labs(title = "Time Complexity Plot" , x = 'number of nodes', y= 'elapsed time [s]') +
  scale_color_identity(name = "Algorithms",
                       breaks = c("#09611D", "#F3C035"),
                       labels = c("Paul Erdos", "Maxcut"),
                       guide = "legend") +
  coord_cartesian(xlim = c(0,201), ylim = c(0,20))







