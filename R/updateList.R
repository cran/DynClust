updateList <-
function
### 'updateList' eliminates the clustered pixels x time from 'pixtoclust'
### a vector containing the indexes of the pixels left for clustering   
(fp.pixtoclust,
### a vector containing the indexes of the pixels left for clustering   
fp.neighinclust
### a vector containing the indexes of the neighboors of the current pixel to clust
##keyword<< internal 
){
  #returns the indexes of clusted pixels x time that are still in the list of pixels x time to clust
  fp.whoisin    <- which(fp.pixtoclust %in% fp.neighinclust)
  #if there are pixels x time to eliminate from list do
  if(length(fp.whoisin)>0){
    #returns the updated list
    return(fp.pixtoclust[-fp.whoisin])}
  #else returns the list
  fp.pixtoclust
  ### a vector containing the indexes of the pixels left for clustering
}
