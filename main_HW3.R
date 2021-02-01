
# Point 1 -----------------------------------------------------------------

load("../Connect-your-Brain/data/hw3_data.RData") # import data

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

asd_td_bin_matrix = rbind(asd_bind_cor, td_bind_cor) 



# APPRACH 2:
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







