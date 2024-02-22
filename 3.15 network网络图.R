library(Hmisc)
library(igraph)
library(ggClusterNet)

setwd('C:/Users/mahoon/Desktop/rawdata/')
bacteria <- read.csv(file = 'table_16s.csv')

taxon_bac <- read.csv(file = 'taxon.csv')
taxon_bac <- taxon_bac[, c('asv', 'Domain', 'Class')]
taxon_bac <- tibble::column_to_rownames(taxon_bac, var = 'asv')

voyage <- levels(factor(bacteria$voyage))

i = 1
for (i in 1:length(voyage)) {
  bac_asv <- bacteria[which(bacteria$voyage == voyage[i]), 19:dim(bacteria)[2]]
  bac_asv <- bac_asv[, which(colSums(bac_asv) > 0)]
  rownames(bac_asv) <- bacteria[which(bacteria$voyage == voyage[i]), ]$id
  
  bac_asv <- prop.table(as.matrix(bac_asv), margin = 1)
  bac_asv <- as.data.frame(t(bac_asv))
  
  taxa1 <- bac_asv
  taxa1[taxa1 > 0] <- 1
  bac_asv <- bac_asv[which(rowSums(taxa1) >= dim(taxa1)[2]/4), ]
  remove(taxa1)
  
  asv_table <- bac_asv
  
  taxa_coor <- rcorr(t(asv_table), type = 'spearman')
  r <- taxa_coor$r
  r[abs(r) < 0.6] <- 0
  p <- taxa_coor$P
  p <- p.adjust(p, method = 'BH')
  p[p >= 0.05] <- -1
  p[p < 0.05 & p >= 0] <- 1
  p[p == -1] <- 0
  z <- r * p
  diag(z) <- 0   
  
  g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
  g <- simplify(g)
  g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
  
  E(g)$correlation <- E(g)$weight
  E(g)$weight <- abs(E(g)$weight)
  E(g)$edge_color <- ifelse(E(g)$correlation > 0, 'blue', ifelse(E(g)$correlation < 0, 'red', 'gray'))
  
  
  taxon <- taxon_bac[as.character(V(g)$name), ]
  V(g)$Domain <- taxon$Domain
  V(g)$Class <- taxon$Class
  V(g)$degree <- degree(g)
  
  com = cluster_fast_greedy(g)
  V(g)$membership <- com$membership
  
#  res <- ZiPiPlot(igraph = g, method = 'cluster_fast_greedy')
  
#  label <- res[[2]][["label"]]
#  zi <- res[[2]][["z"]]
#  pi <- res[[2]][["p"]]
#  zi_pi <- data.frame(label, zi, pi)
  
#  P_zi_pi <- ggplot(zi_pi, aes(x = pi, y = zi)) +
#    geom_point(size = 3, alpha =0.3) +
#    geom_label_repel(label = label) +
#    geom_hline(aes(yintercept = 2.5), linetype = 'dashed') +
#    geom_vline(aes(xintercept = 0.62), linetype = 'dashed') +
#    theme_test()
#  ggsave(filename =  paste0('C:/Users/mahoon/Desktop/', voyage[i], '_keystone.pdf'), dpi = 600)
#  write.csv(zi_pi, file = paste0('C:/Users/mahoon/Desktop/', voyage[i], '_zi_pi_bac.csv'))
  write.graph(g, file = paste0('C:/Users/mahoon/Desktop/', voyage[i], '-bac-0.6-25%.graphml'), format = 'graphml')
  
  net_prop <- as.data.frame(net_properties.2(g))
  write.csv(net_prop, file = paste0('C:/Users/mahoon/Desktop/', voyage[i], '_net_prop.csv'))
}

###1000次随机取样

library(Hmisc)
library(igraph)
library(ggClusterNet)

net_properties.3 <- function (igraph, n.hub = FALSE) 
{
  num.edges <- length(E(igraph))
  num.edges
  num.vertices <- length(V(igraph))
  num.vertices
  connectance <- edge_density(igraph, loops = FALSE)
  average.degree <- mean(igraph::degree(igraph))
  average.degree
  if (!is.null(E(igraph)$weight)) {
    igraph.weight <- E(igraph)$weight
    E(igraph)$weight = abs(E(igraph)$weight)
  }
  average.path.length <- average.path.length(igraph)
  average.path.length
  diameter <- diameter(igraph, directed = FALSE, unconnected = TRUE, 
                       weights = NULL)
  diameter
  if (!is.null(E(igraph)$weight)) {
    E(igraph)$weight = igraph.weight
  }
  edge.connectivity <- edge_connectivity(igraph)
  edge.connectivity
  clustering.coefficient <- transitivity(igraph, type = "average")
  clustering.coefficient
  no.clusters <- no.clusters(igraph)
  no.clusters
  centralization.degree <- centralization.degree(igraph)$centralization
  centralization.degree
  centralization.betweenness <- centralization.betweenness(igraph)$centralization
  centralization.betweenness
  centralization.closeness <- centralization.closeness(igraph)$centralization
  centralization.closeness
  if (!is.null(E(igraph)$weight)) {
    num.pos.edges <- sum(E(igraph)$correlation > 0)
    num.neg.edges <- sum(E(igraph)$correlation < 0)
  }
  else {
    num.pos.edges <- 0
    num.neg.edges <- 0
  }
  modularity_igraph = function(net, method = "cluster_walktrap") {
    if (method == "cluster_walktrap") {
      fc <- cluster_walktrap(net)
    }
    if (method == "cluster_edge_betweenness") {
      fc <- cluster_edge_betweenness(net)
    }
    if (method == "cluster_fast_greedy") {
      fc <- cluster_fast_greedy(net)
    }
    if (method == "cluster_spinglass") {
      fc <- cluster_spinglass(net)
    }
    modularity <- modularity(net, membership(fc))
    return(modularity)
  }
  mod1 = modularity_igraph(igraph, method = "cluster_walktrap")
  rand.g <- erdos.renyi.game(length(V(igraph)), length(E(igraph)), 
                             type = "gnm")
  mod2 = modularity_igraph(rand.g, method = "cluster_walktrap")
  RM = (mod1 - mod2)/mod2
  if (n.hub) {
    res = ZiPiPlot(igraph = igraph, method = "cluster_walktrap")
    data = res[[2]]
    head(data)
    n.hub = data$roles[data$roles != "Peripherals"] %>% 
      length()
  }
  else {
    n.hub = "Not.calculated"
  }
  igraph.network.pro <- rbind(num.edges, num.pos.edges, num.neg.edges, 
                              num.vertices, connectance, average.degree, average.path.length, 
                              diameter, edge.connectivity, clustering.coefficient, 
                              no.clusters, centralization.degree, centralization.betweenness, 
                              centralization.closeness, RM, n.hub)
  rownames(igraph.network.pro) <- c("num.edges(L)", "num.pos.edges", 
                                    "num.neg.edges", "num.vertices(n)", "Connectance(edge_density)", 
                                    "average.degree(Average K)", "average.path.length", 
                                    "diameter", "edge.connectivity", "mean.clustering.coefficient(Average.CC)", 
                                    "no.clusters", "centralization.degree", "centralization.betweenness", 
                                    "centralization.closeness", "RM(relative.modularity)", 
                                    "the.number.of.keystone.nodes")
  colnames(igraph.network.pro) <- "value"
  return(igraph.network.pro)
}

######

bacteria <- read.csv(file = 'table_16s.csv')
regional <- levels(factor(bacteria$voyage))

##############################

net_propertice <- as.data.frame(matrix(nrow = 16, ncol = 0))

i = 4
23*22*21*20/(4*3*2*1)

aim_regional <- bacteria[which(bacteria$voyage == regional[i]), ]
bac_asv_start <- aim_regional[, 19:dim(aim_regional)[2]]
bac_asv_start <- bac_asv_start[, which(colSums(bac_asv_start) > 0)]
rownames(bac_asv_start) <- aim_regional$id

for (j in 1:1000) {
  bac_asv <- bac_asv_start[sample(1:dim(bac_asv_start)[1], 19, replace = FALSE), ]
  bac_asv <- prop.table(as.matrix(bac_asv), margin = 1)
  bac_asv <- as.data.frame(t(bac_asv))
  taxa1 <- bac_asv
  taxa1[taxa1 > 0] <- 1
  bac_asv <- bac_asv[which(rowSums(taxa1) >= dim(taxa1)[2]/4), ]
  taxa_coor <- rcorr(t(bac_asv), type = 'spearman')
  r <- taxa_coor$r
  r[abs(r) < 0.6] <- 0
  p <- taxa_coor$P
  p <- p.adjust(p, method = 'BH')
  p[p >= 0.05] <- -1
  p[p < 0.05 & p >= 0] <- 1
  p[p == -1] <- 0
  z <- r * p
  diag(z) <- 0   
  g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
  g <- simplify(g)
  g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
  E(g)$correlation <- E(g)$weight
  E(g)$weight <- abs(E(g)$weight)
  net_prop <- as.data.frame(net_properties.3(g))
  colnames(net_prop) <- paste0(regional[i], j)
  net_propertice <- cbind(net_propertice, net_prop)
}

write.csv(net_propertice, file = paste0('C:/Users/mahoon/Desktop/', regional[i], '.csv'))
net_propertice <- as.data.frame(matrix(nrow = 16, ncol = 0))



