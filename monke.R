setwd("C:/Users/mlove/Desktop/MSc-Paleobio/semester-2/Phylogenetics/Project/dataset")
library(phytools)
library(qpcR)
library(geiger)
library(adiv)

library(DiagrammeR) #packages for visual construction
library(flextable)
library(officer)
library(reshape2)
library(ggplot2)

#set seed
set.seed(2)

#load in data
characters<-read.csv("Data-Data.csv")
row.names(characters) <- characters$Species
characters <- lapply(characters[,4:15], setNames, nm = rownames(characters))
characters[[12]]

#load in nexus files
#trees <- read.nexus("TreeBlock_10kTrees_Primates_Version3.nex")
tree <- read.nexus("consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
tree

#crossreference species
chk<-name.check(tree,characters[[1]])
summary(chk)

################################################################################
#plot first trait
group_size<-setNames(characters[[1]],
                  names(characters[[1]]))
plot.phylo(tree,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Group Size Trait Distribution")
X<-strsplit(setNames(as.character(group_size),names(group_size)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(tree),4,
             dimnames=list(tree$tip.label,
                           c("0","1","2","3")))
for(i in 1:Ntip(tree)){ 
  pies[tree$tip.label[i],
       X[[tree$tip.label[i]]]]<-
  rep(1/length(X[[tree$tip.label[i]]]),
      length(X[[tree$tip.label[i]]]))
}
cols<-c("red", "green", "blue", "khaki")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Solitary","2-15","16-30",">30"),
       pt.cex=1, cex=0.7,pch=16,col=cols)


#fit to model (ordered T/F)
# Fit ER model
group.ER<-fitpolyMk(tree,group_size,
                    model="ER",ordered=FALSE)
# Fit SYM model
group.SYM <-fitpolyMk(tree,group_size,
                      model="SYM",ordered=FALSE)

# Fit ARD model
group.ARD<-fitpolyMk(tree,group_size,
                     model="ARD",ordered=FALSE)

# Compare all AICs together
aic_table <- AIC(group.ARD,group.SYM,group.ER)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
group.model<-fitpolyMk(tree,group_size,
                             model="ER",ordered=FALSE)
plot(group.model)
title(main=
        "Fitted polymorphic trait evolution model\nfor primate group size",
      font.main=3)

#stochastic
Q<-as.Qmatrix(group.model)
print(Q,digits=3)
X<-group.model$data
head(X)
group.maps<-make.simmap(tree,x=X,Q=Q,
                            nsim=100, ncores=4)
group.maps

obj<-summary(group.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(group.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))
# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main="Reconstructed Phylogenetic Tree of Primate Group Size")
leg_group_map <- c("0" = "Solitary", "1" = "2-15", "2" = "16-30", "3" = ">30")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(legend_groups) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = legend_groups[present_states],
       col = map.cols[present_states],
       pt.cex=1, cex=0.7,pch=16, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)

#calculate signal
logL_obs1 <- logLik(group.model)

# Permutation test
n_sim <- 1000
logL_rand1 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  group_rand <- setNames(sample(group_size), names(group_size))
  
  # Fit the same model
  fit_rand1 <- fitpolyMk(tree, group_rand, model="ER", ordered=FALSE)
  logL_rand1[i] <- logLik(fit_rand1)
}
# Compare observed vs null distribution
p_value <- (sum(logL_rand1 >= logL_obs1)+1)/(length(logL_rand1)+1)

cat("Observed logLik:", logL_obs1, "\n")
cat("Mean randomized logLik:", mean(logL_rand1), "\n")
cat("p-value for phylogenetic signal:", p_value, "\n")


##############################################################################################
#plot second trait
seasonal_breeding<-setNames(characters[[2]],
                     names(characters[[2]]))

plot.phylo(tree,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Seasonal Breeding Trait Distribution")
X<-strsplit(setNames(as.character(seasonal_breeding),names(seasonal_breeding)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(tree),2,
             dimnames=list(tree$tip.label,
                           c("0","1")))
for(i in 1:Ntip(tree)){ 
  pies[tree$tip.label[i],
       X[[tree$tip.label[i]]]]<-
    rep(1/length(X[[tree$tip.label[i]]]),
        length(X[[tree$tip.label[i]]]))
}
cols<-c("red", "blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Non-seasonal", "Seasonal"),
       pt.cex=1, cex=0.7,pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
seasonal.ER<-fitpolyMk(tree,seasonal_breeding,
                    model="ER",ordered=TRUE)
# Fit SYM model
seasonal.SYM <-fitpolyMk(tree,seasonal_breeding,
                      model="SYM",ordered=TRUE)

# Fit ARD model
seasonal.ARD<-fitpolyMk(tree,seasonal_breeding,
                     model="ARD",ordered=TRUE)

# Compare all together
aic_table <- AIC(seasonal.ARD,seasonal.SYM,seasonal.ER)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
plot(seasonal.ARD)
title(main=
        "Fitted polymorphic trait evolution model\nfor primate group size",
      font.main=3)

#stochastic
Q<-as.Qmatrix(seasonal.ARD)
print(Q,digits=3)
X<-seasonal.ARD$data
head(X)
seasonal.maps<-make.simmap(tree,x=X,Q=Q,
                        nsim=100, ncores=4)
seasonal.maps

obj<-summary(seasonal.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(seasonal.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main= "Phylogenetic Tree of Primate Seasonal Breeding Behavior")
leg_group_map <- c("0"="Non-Seasonal","1"="Seasonal", "1+2"="Variable")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(leg_group_map) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = leg_group_map[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)


#calculate signal
logL_obs2 <- logLik(seasonal.ARD)

# Permutation test
n_sim <- 1000
logL_rand2 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  season_rand <- setNames(sample(seasonal_breeding), names(seasonal_breeding))
  
  # Fit the same model
  fit_rand2 <- fitpolyMk(tree, season_rand, model="ARD", ordered=TRUE)
  logL_rand2[i] <- logLik(fit_rand2)
}
# Compare observed vs null distribution
p_value2 <- (sum(logL_rand2 >= logL_obs2)+1)/(length(logL_rand2)+1)

cat("Observed logLik:", logL_obs2, "\n")
cat("Mean randomized logLik:", mean(logL_rand2), "\n")
cat("p-value for phylogenetic signal:", p_value2, "\n")


#################################################################################
#plot third trait
mating_system<-setNames(characters[[3]],
                     names(characters[[3]]))

plot.phylo(tree,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Mating System Trait Distribution")
X<-strsplit(setNames(as.character(mating_system),names(mating_system)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(tree),4,
             dimnames=list(tree$tip.label,
                           c("0","1","2","3")))
for(i in 1:Ntip(tree)){ 
  pies[tree$tip.label[i],
       X[[tree$tip.label[i]]]]<-
    rep(1/length(X[[tree$tip.label[i]]]),
        length(X[[tree$tip.label[i]]]))
}
cols<-c("red", "green", "blue", "khaki")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Monogomous","Polygynous","Polyandrous","Promiscuous"),
       pt.cex=1,cex=0.6, pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
mating.ER<-fitpolyMk(tree,mating_system,
                       model="ER",ordered=FALSE)
# Fit SYM model
mating.SYM <-fitpolyMk(tree,mating_system,
                         model="SYM",ordered=FALSE)

# Fit ARD model
mating.ARD<-fitpolyMk(tree,mating_system,
                        model="ARD",ordered=FALSE)

# Compare all together
aic_table <- AIC(mating.ARD, mating.SYM, mating.ER)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
plot(mating.ER)
title(main=
        "Fitted polymorphic trait evolution model\nfor primate mating system",
      font.main=3)

#stochastic
Q<-as.Qmatrix(mating.ER)
print(Q,digits=3)
X<-mating.ER$data
head(X)
mating.maps<-make.simmap(tree,x=X,Q=Q,
                           nsim=100, ncores=4)
mating.maps

obj<-summary(mating.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(mating.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main= "Phylogenetic Tree of Primate Mating System")
leg_group_map <- c("0"="Monogomous","1"="Polygynous","2"="Polyandrous","3"="Promiscuous")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(legend_groups) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = legend_groups[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)


#calculate signal
logL_obs3 <- logLik(mating.ER)

# Permutation test
n_sim <- 1000
logL_rand3 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  mating_rand <- setNames(sample(mating_system), names(mating_system))
  
  # Fit the same model
  fit_rand3 <- fitpolyMk(tree, mating_rand, model="ER", ordered=FALSE)
  logL_rand3[i] <- logLik(fit_rand3)
}
# Compare observed vs null distribution
p_value3 <- (sum(logL_rand3 >= logL_obs3)+1)/(length(logL_rand3)+1)

cat("Observed logLik:", logL_obs3, "\n")
cat("Mean randomized logLik:", mean(logL_rand3), "\n")
cat("p-value for phylogenetic signal:", p_value3, "\n")


################################################################################
#plot fourth trait
repro_age<-setNames(characters[[4]],
                        names(characters[[4]]))

par(mar=c(1,1,1,1))
plot.phylo(tree,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Reproductive Age Trait Distribution")
X <- setNames(as.character(characters[[4]]), names(characters[[4]]))
pies<-matrix(0,Ntip(tree),5,
             dimnames=list(tree$tip.label,
                           c("0","1","2","3","4")))
for(i in 1:Ntip(tree)){ 
  tip <- tree$tip.label[i]   # skip if missing
  state <- X[[tip]]                # single state (e.g. "0","1",…)
  pies[tip, state] <- 1
}
cols<-c("red", "green", "blue", "khaki", "orange")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("<3 y.o.","3.1-4.0 y.o.","4.1-6.0 y.o.","6.1-10.0 y.o.",">10.1 y.o."),
       pt.cex=0.8, cex=0.6,pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
repro.ER<-fitMk(tree,repro_age,
                     model="ER",ordered=TRUE)
# Fit SYM model
repro.SYM <-fitMk(tree,repro_age,
                       model="SYM",ordered=TRUE)

# Fit ARD model
repro.ARD<-fitMk(tree,repro_age,
                      model="ARD",ordered=TRUE)

# Compare all together
aic_table <- AIC(repro.ARD, repro.SYM, repro.ER)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
plot(repro.SYM)
title(main=
        "Fitted polymorphic trait evolution model\nfor primate mating system",
      font.main=3)

#stochastic
Q<-as.Qmatrix(repro.SYM)
print(Q,digits=3)
X<-repro.SYM$data
head(X)
repro.maps<-make.simmap(tree,x=X,Q=Q,
                         nsim=100, ncores=4)
repro.maps

obj<-summary(repro.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(repro.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main="Phylogenetic Tree of Primate Reproductive Age")
leg_group_map <- c("0"="<3 y.o.","1"="3.1-4.0 y.o.","2"="4.1-6.0 y.o.","3"="6.1-10.0 y.o.", "4"=">10.1 y.o.")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(legend_groups) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = legend_groups[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)


#calculate signal
logL_obs4 <- logLik(repro.SYM)

# Permutation test
n_sim <- 298
logL_rand4 <- numeric(n_sim)
for (i in 1:n_sim) {
  # Randomize states across the tree tips
  repro_rand <- setNames(sample(repro_age), names(repro_age))
  
  # Fit the same model
  fit_rand4 <- fitMk(tree, repro_rand, model="SYM", ordered=TRUE)
  logL_rand4[i] <- logLik(fit_rand4)
  
  
}
# Compare observed vs null distribution
p_value4 <- (sum(logL_rand4 >= logL_obs4)+1)/(length(logL_rand4)+1)

cat("Observed logLik:", logL_obs4, "\n")
cat("Mean randomized logLik:", mean(logL_rand4), "\n")
cat("p-value for phylogenetic signal:", p_value4, "\n")



################################################################################
#plot fifth trait
rankbias<-setNames(characters[[5]],
                        names(characters[[5]]))

plot.phylo(tree,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Male Rank Bias Trait Distribution")
X<-strsplit(setNames(as.character(rankbias),names(rankbias)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(tree),2,
             dimnames=list(tree$tip.label,
                           c("0","1")))
for(i in 1:Ntip(tree)){ 
  pies[tree$tip.label[i],
       X[[tree$tip.label[i]]]]<-
    rep(1/length(X[[tree$tip.label[i]]]),
        length(X[[tree$tip.label[i]]]))
}
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Male Rank Bias Absent", "Male Rank Bias Present"),
       pt.cex=1,cex=0.6, pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
rank.ER<-fitpolyMk(tree,rankbias,
                     model="ER",ordered=TRUE)
# Fit SYM model
rank.SYM <-fitpolyMk(tree,rankbias,
                       model="SYM",ordered=TRUE)

# Fit ARD model
rank.ARD<-fitpolyMk(tree,rankbias,
                      model="ARD",ordered=TRUE)

# Compare all together
aic_table <- AIC(rank.ER,rank.SYM,rank.ARD)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
plot(rank.ARD)
title(main=
        "Fitted polymorphic trait evolution model\nfor primate male rank bias",
      font.main=3)

#stochastic
Q<-as.Qmatrix(rank.ARD)
print(Q,digits=3)
X<-rank.ARD$data
head(X)
rank.maps<-make.simmap(tree,x=X,Q=Q,
                         nsim=100, ncores=4)
rank.maps

obj<-summary(rank.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(rank.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main="Ancestral Reconstruction of Male Rank Bias")
leg_group_map <- c("0"="Bias Absent","1"="Bias Present", "1+2"="Variable")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(leg_group_map) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = leg_group_map[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)


#calculate signal
logL_obs5 <- logLik(rank.ARD)

# Permutation test
n_sim <- 1000
logL_rand5 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  rank_rand <- setNames(sample(rankbias), names(rankbias))
  
  # Fit the same model
  fit_rand5 <- fitpolyMk(tree, rank_rand, model="ER", ordered=TRUE)
  logL_rand5[i] <- logLik(fit_rand5)
}
# Compare observed vs null distribution
p_value5 <- (sum(logL_rand5 >= logL_obs5)+1)/(length(logL_rand5)+1)

cat("Observed logLik:", logL_obs5, "\n")
cat("Mean randomized logLik:", mean(logL_rand5), "\n")
cat("p-value for phylogenetic signal:", p_value5, "\n")


################################################################################
#plot sixth trait
paternal <- setNames(characters[[6]], names(characters[[6]]))

pies <- matrix(0, Ntip(tree_clean), 3,
               dimnames = list(tree_clean$tip.label, c("0","1","2")))

#omit species where data is missing (?)
valid_tips <- names(paternal)[paternal != "?"]

#prune the tree by cross-referencing the tree tips and the valid species
tree_clean <- drop.tip(tree, setdiff(tree$tip.label, valid_tips))

#redefine the character matrix to only include the valid species
paternal_clean <- paternal[valid_tips]

for (i in 1:Ntip(tree_clean)) {
  tip <- tree_clean$tip.label[i]
  
  if (!(tip %in% names(paternal))) next   # skip if no data
  if (paternal[[tip]] == "?") next        # skip missing
  
  state <- as.character(paternal[[tip]])
  pies[tip, state] <- 1
}

plot.phylo(tree_clean,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Male Care of Young Trait Distribution")
cols<-c("red", "green", "blue", "khaki", "orange")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Absent","Infrequent","Frequent"),
       pt.cex=0.8, cex=0.6,pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
paternal.ER<-fitMk(tree_clean,paternal_clean,
                model="ER",ordered=TRUE)
# Fit SYM model
paternal.SYM <-fitMk(tree_clean,paternal_clean,
                  model="SYM",ordered=TRUE)

# Fit ARD model
paternal.ARD<-fitMk(tree_clean,paternal_clean,
                 model="ARD",ordered=TRUE)

# Compare all together
aic_table <- AIC(paternal.ER,paternal.SYM,paternal.ARD)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
plot(paternal.ER)
title(main=
        "Fitted polymorphic trait evolution model\nfor primate paternal care",
      font.main=3)

#stochastic
Q<-as.Qmatrix(paternal.ER)
print(Q,digits=3)
X<-paternal.ER$data
head(X)
paternal.maps<-make.simmap(tree_clean,x=X,Q=Q,
                        nsim=100, ncores=4)
paternal.maps

obj<-summary(paternal.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(paternal.maps)

#Get the child node for each edge
child_nodes <- tree_clean$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree_clean),
                        tree_clean$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree_clean, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main="Ancestral Reconstruction of Male Care of Young")
leg_group_map <- c("0"="Absent","1"="Infrequent","2"="Frequent")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(legend_groups) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = legend_groups[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)



#calculate signal
logL_obs6 <- logLik(paternal.ER)

# Permutation test
n_sim <- 1000
logL_rand6 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  rank_rand <- setNames(sample(paternal_clean), names(paternal_clean))
  
  # Fit the same model
  fit_rand6 <- fitMk(tree_clean, rank_rand, model="ER", ordered=TRUE)
  logL_rand6[i] <- logLik(fit_rand6)
}
# Compare observed vs null distribution
p_value6 <- (sum(logL_rand6 >= logL_obs6)+1)/(length(logL_rand6)+1)

cat("Observed logLik:", logL_obs6, "\n")
cat("Mean randomized logLik:", mean(logL_rand6), "\n")
cat("p-value for phylogenetic signal:", p_value6, "\n")

 
################################################################################
#plot seventh trait
allo <- setNames(characters[[7]], names(characters[[7]]))
plot.phylo(tree_clean,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Alloparenting Trait Distribution")
X<-strsplit(setNames(as.character(allo),names(allo)),"+",
            fixed=TRUE)
pies <- matrix(0, Ntip(tree), 3,
               dimnames = list(tree$tip.label, c("0","1","2")))

for(i in 1:Ntip(tree)){ 
  pies[tree$tip.label[i],
       X[[tree$tip.label[i]]]]<-
    rep(1/length(X[[tree$tip.label[i]]]),
        length(X[[tree$tip.label[i]]]))
}

cols<-c("red", "green", "blue", "khaki", "orange")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Absent","Infrequent","Frequent"),
       pt.cex=0.8, cex=0.6,pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
allo.ER<-fitpolyMk(tree,allo,
                   model="ER",ordered=TRUE)
# Fit SYM model
allo.SYM <-fitpolyMk(tree,allo,
                     model="SYM",ordered=TRUE)

# Fit ARD model
allo.ARD<-fitpolyMk(tree,allo,
                    model="ARD",ordered=TRUE)

# Compare all together
aic_table <- AIC(allo.ER,allo.SYM,allo.ARD)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
plot(allo.ER)
title(main=
        "Fitted polymorphic trait evolution model\nfor primate allopatric care",
      font.main=3)

#stochastic
Q<-as.Qmatrix(allo.ER)
print(Q,digits=3)
X<-allo.ER$data
head(X)
allo.maps<-make.simmap(tree,x=X,Q=Q,
                           nsim=100, ncores=4)
allo.maps

obj<-summary(allo.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(allo.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree_clean, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main="Ancestral Reconstruction of Alloparenting Behavior")
leg_group_map <- c("0"="Absent","1"="Infrequent","2"="Frequent")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(legend_groups) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = legend_groups[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)


#calculate signal
logL_obs7 <- logLik(allo.ER)

# Permutation test
n_sim <- 1000
logL_rand7 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  allo_rand <- setNames(sample(allo), names(allo))
  
  # Fit the same model
  fit_rand7 <- fitpolyMk(tree, allo_rand, model="ER", ordered=TRUE)
  logL_rand7[i] <- logLik(fit_rand7)
}
# Compare observed vs null distribution
p_value7 <- (sum(logL_rand7 >= logL_obs7)+1)/(length(logL_rand7)+1)

cat("Observed logLik:", logL_obs7, "\n")
cat("Mean randomized logLik:", mean(logL_rand7), "\n")
cat("p-value for phylogenetic signal:", p_value7, "\n")



################################################################################
#plot eigth trait
male <- setNames(characters[[8]], names(characters[[8]]))
pies <- matrix(0, Ntip(tree), 2,
               dimnames = list(tree$tip.label, c("0","1")))

for (i in 1:Ntip(tree)) {
  tip <- tree$tip.label[i]
  state <- as.character(male[[tip]])
  pies[tip, state] <- 1
}

plot.phylo(tree_clean,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Male Same-Sex Sexual Behavior Trait Distribution")
cols<-c("red", "green", "blue", "khaki", "orange")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Not Observed","Observed"),
       pt.cex=0.8, cex=0.6,pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
male.ER<-fitMk(tree,male,
                   model="ER",ordered=FALSE)
# Fit SYM model
male.SYM <-fitMk(tree,male,
                     model="SYM",ordered=FALSE)

# Fit ARD model
male.ARD<-fitMk(tree,male,
                    model="ARD",ordered=FALSE)

# Compare all together
aic_table <- AIC(male.ER,male.SYM,male.ARD)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use simplest small AIC
plot(male.ARD)
title(main=
        "Fitted polymorphic trait evolution model\nfor male primate same-sex sexual behavior",
      font.main=3)

#stochastic
Q<-as.Qmatrix(male.ARD)
print(Q,digits=3)
X<-male.ARD$data
head(X)
male.maps<-make.simmap(tree,x=X,Q=Q,
                       nsim=100, ncores=4)
male.maps

obj<-summary(male.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(male.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree_clean, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main="Ancestral Reconstruction of Male Same-Sex Behavior")
leg_group_map <- c("0"="Not Observed","1"="Observed")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(legend_groups) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = legend_groups[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)


#calculate signal
logL_obs8 <- logLik(male.ARD)

# Permutation test
n_sim <- 1000
logL_rand8 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  male_rand <- setNames(sample(male), names(male))
  
  # Fit the same model
  fit_rand8 <- fitMk(tree, male_rand, model="ARD", ordered=FALSE)
  logL_rand8[i] <- logLik(fit_rand8)
}
# Compare observed vs null distribution
p_value8 <- (sum(logL_rand8 >= logL_obs8)+1)/(length(logL_rand8)+1)

cat("Observed logLik:", logL_obs8, "\n")
cat("Mean randomized logLik:", mean(logL_rand8), "\n")
cat("p-value for phylogenetic signal:", p_value8, "\n")



################################################################################
#plot ninth trait
female <- setNames(characters[[9]], names(characters[[9]]))
pies <- matrix(0, Ntip(tree), 2,
               dimnames = list(tree$tip.label, c("0","1")))

for (i in 1:Ntip(tree)) {
  tip <- tree$tip.label[i]
  state <- as.character(female[[tip]])
  pies[tip, state] <- 1
}

plot.phylo(tree_clean,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Female Same-Sex Behavior Trait Distribution")
cols<-c("red", "green", "blue", "khaki", "orange")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Not Observed","Observed"),
       pt.cex=0.8, cex=0.6,pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
female.ER<-fitMk(tree,female,
               model="ER",ordered=FALSE)
# Fit SYM model
female.SYM <-fitMk(tree,female,
                 model="SYM",ordered=FALSE)

# Fit ARD model
female.ARD<-fitMk(tree,female,
                model="ARD",ordered=FALSE)

# Compare all together
aic_table <- AIC(female.ER,female.SYM,female.ARD)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
plot(female.ER)
title(main=
        "Fitted polymorphic trait evolution model\nfor male primate same-sex sexual behavior",
      font.main=3)

#stochastic
Q<-as.Qmatrix(female.ER)
print(Q,digits=3)
X<-female.ER$data
head(X)
female.maps<-make.simmap(tree,x=X,Q=Q,
                       nsim=100, ncores=4)
female.maps

obj<-summary(female.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(female.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree_clean, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main="Ancestral Reconstruction of Female Same-Sex Sexual Behaviors")
leg_group_map <- c("0"="Not Observed","1"="Observed")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(legend_groups) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = legend_groups[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)


#calculate signal
logL_obs9 <- logLik(female.ER)

# Permutation test
n_sim <- 1000
logL_rand9 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  female_rand <- setNames(sample(female), names(female))
  
  # Fit the same model
  fit_rand9 <- fitMk(tree, female_rand, model="ER", ordered=FALSE)
  logL_rand9[i] <- logLik(fit_rand9)
}
# Compare observed vs null distribution
p_value9 <- (sum(logL_rand9 >= logL_obs9)+1)/(length(logL_rand9)+1)

cat("Observed logLik:", logL_obs9, "\n")
cat("Mean randomized logLik:", mean(logL_rand9), "\n")
cat("p-value for phylogenetic signal:", p_value9, "\n")




################################################################################
#plot tenth trait
trans <- setNames(characters[[10]], names(characters[[10]]))
pies <- matrix(0, Ntip(tree), 2,
               dimnames = list(tree$tip.label, c("0","1")))

for (i in 1:Ntip(tree)) {
  tip <- tree$tip.label[i]
  state <- as.character(trans[[tip]])
  pies[tip, state] <- 1
}

plot.phylo(tree_clean,type="phylogram", cex=0.7, show.tip.label=TRUE, 
           label.offset = 0.007, main="Trangender Behavior Trait Distribution")
cols<-c("red", "green", "blue", "khaki", "orange")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.3)
par(fg="black")
legend(x="bottomleft",legend=c("Not Observed","Observed"),
       pt.cex=0.8, cex=0.6,pch=16,col=cols)

#fit to model (ordered T/F)
# Fit ER model
trans.ER<-fitMk(tree,trans,
               model="ER",ordered=FALSE)
# Fit SYM model
trans.SYM <-fitMk(tree,trans,
                 model="SYM",ordered=FALSE)

# Fit ARD model
trans.ARD<-fitMk(tree,trans,
                model="ARD",ordered=FALSE)

# Compare all together
aic_table <- AIC(trans.ER,trans.SYM,trans.ARD)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table$weights <- akaike.weights(aic_table$AIC)$weights

print(aic_table)

#use smallest AIC
plot(trans.ARD)
title(main=
        "Fitted polymorphic trait evolution model\nfor primate trangender behavior",
      font.main=3)

#stochastic
Q<-as.Qmatrix(trans.ARD)
print(Q,digits=3)
X<-trans.ARD$data
head(X)
trans.maps<-make.simmap(tree,x=X,Q=Q,
                       nsim=100, ncores=4)
trans.maps

obj<-summary(trans.maps)
groups <- colnames(obj$ace)
map.cols<-setNames(rainbow(length(groups)), groups)

# Summarize the simmaps to get average time per state on each edge
summary_obj <- summary(trans.maps)

#Get the child node for each edge
child_nodes <- tree$edge[, 2]
edge_children <- ifelse(child_nodes <= Ntip(tree),
                        tree$tip.label[child_nodes],  # Tips: get label
                        as.character(child_nodes))

#Reorder summary_obj$ace to match tree edge order
edge_state_times <- summary_obj$ace[as.character(edge_children), ]

modal_states <- apply(edge_state_times, 1, function(x) names(which.max(x)))

# Assign colors
edge_colors <- map.cols[modal_states]

# Plot the consensus tree with modal-state edge colors
plot.phylo(tree, edge.color=edge_colors, show.tip.label=TRUE, cex=0.6,
           main="Ancestral Reconstruction of Transgender Behavior")
leg_group_map <- c("0"="Not Observed","1"="Observed")
legend_groups <- sapply(strsplit(groups, "\\+"), function(parts) {
  paste(leg_group_map[parts], collapse = " + ")
})
names(legend_groups) <- groups 
present_states <- intersect(groups, unique(modal_states))
legend("bottomleft",
       legend = legend_groups[present_states],
       pch = 15,
       col = map.cols[present_states],
       pt.cex = 0.7, cex = 0.6, bty = "n", ncol = 2, text.col = "black")

# Node pies still from the same summary
nodelabels(pie=summary_obj$ace[1:tree$Nnode,], piecol=map.cols, cex=0.3)


#calculate signal
logL_obs10 <- logLik(trans.ARD)

# Permutation test
n_sim <- 1000
logL_rand10 <- numeric(n_sim)

for (i in 1:n_sim) {
  # Randomize states across the tree tips
  trans_rand <- setNames(sample(trans), names(trans))
  
  # Fit the same model
  fit_rand10 <- fitMk(tree, trans_rand, model="ARD", ordered=FALSE)
  logL_rand10[i] <- logLik(fit_rand10)
}
# Compare observed vs null distribution
p_value10 <- (sum(logL_rand10 >= logL_obs10)+1)/(length(logL_rand10)+1)

cat("Observed logLik:", logL_obs10, "\n")
cat("Mean randomized logLik:", mean(logL_rand10), "\n")
cat("p-value for phylogenetic signal:", p_value10, "\n")


#####################################

#adjust for multiple comparisons
p_values_vec <- c(p_value,p_value2,p_value3,p_value4,p_value5,p_value6,p_value7,p_value8,p_value9,p_value10)  # vector of 10 permutation p-values
p_holm <- p.adjust(p_values_vec, method = "holm")
alpha <- 0.05
sig_holm <- p_holm < alpha

data.frame(p_values_vec, p_holm, sig_holm)

###############################
grViz("
digraph literature_research_flow {
  graph [rankdir = TB, labelloc = t, fontsize = 18, fontname = Arial,
         label=< <B>Literature Research Workflow</B> >]

  node [shape = box, style = 'rounded,filled', fillcolor = lavender, fontname = Arial]

  # Start and End nodes
  start [label=<<B>Start</B>> fillcolor=lightblue]
  end   [label=<<B>End</B>> fillcolor=lightblue]

  # Workflow nodes
  define [label=<
    <B>Define Traits Clearly</B><BR/>
    • Scientific definitions<BR/>
    • Quantify/define discrete character states<BR/>
    • Define binary, multi-state, polymorph criteria
  >]

  search [label=<
    <B>Systematic Literature Search</B><BR/>
    • Databases: Scopus, Google Scholar<BR/>
    • Keywords = species + trait
  >]

  screen [label=<
    <B>Screen &amp; Select Studies</B><BR/>
    • Exclude captive-only behavior<BR/>
    • Exclude behavior only seen in adolescents<BR/>
    • Exclude anecdotal/unpublished reports
  >]

  assess [label=<
    <B>Assess Reliability</B><BR/>
    • Ensure behaviors are observed repeatedly within a population<BR/>
    • Ensure multiple reports across populations<BR/>
    • Study length of one full season or more
  >]

  organize [label=<
    <B>Organize in Character Matrix</B><BR/>
    • Species &#215; Traits &#215; Source<BR/>
    • Mark uncertainty (unknown vs. absent)
  >]

  # Workflow edges
  start -> define
  define -> search
  search -> screen
  screen -> assess
  assess -> organize
  organize -> end
}
")
###############################################
trait_info <- data.frame(
  Trait = c("Group Size", "Seasonal Breeding", "Mating System", "Reproductive Age",
            "Male Rank Bias", "Male Care", "Alloparenting", 
            "Male SSB", "Female SSB", "Transgender"),
  Data_Type = c("Categorical (ordinal)", "Binary", "Categorical (nominal)", "Continuous (binned)",
                "Binary", "Ordinal", "Ordinal", "Binary", "Binary", "Binary"),
  Polymorphic = c("Yes","Yes","Yes","No","Yes","No","Yes","Yes","Yes","Yes"),
  stringsAsFactors = FALSE
)

# build flextable
ft <- flextable(trait_info)

caption_par <- as_paragraph(
  as_chunk("Table 1. Summary of the behavioral traits used in the analysis, including the trait name, data type, and polymorphic status.", 
           props = fp_text(italic = TRUE, font.size = 11)))

ft <- set_caption(ft, caption = caption_par)

# some cosmetic tweaks
ft <- bold(ft, part = "header")            # bold header row
ft <- autofit(ft)                          # adjust column widths
ft <- theme_vanilla(ft)                    # nice default theme

# show in viewer (RStudio will display it)
ft

#############################################
model_info <- data.frame(
  Trait = c("Group Size", "Seasonal Breeding", "Mating System", "Reproductive Age",
            "Male Rank Bias", "Male Care", "Alloparenting", 
            "Male SSB", "Female SSB", "Transgender"),
  ER_Model_AIC = c("164.737","81.191", "161.132","99.253","78.952","52.873","92.571","43.180","30.476","35.723"),
  SYM_Model_AIC = c("194.079","82.532", "180.676","94.570","80.864","54.564","93.831","43.180","30.476","35.723"),
  ARD_Model_AIC = c("237.576","78.093", "221.714","111.245","77.744","53.914","94.655","43.078","37.037","34.521"),
  Ordered_Model = c("No","Yes", "No","Yes","Yes","Yes","Yes","No","No","No"),
  stringsAsFactors = FALSE
)

# build flextable
ft <- flextable(model_info)

caption_par <- as_paragraph(
  as_chunk("Table 2. Comparison of AIC values for each trait across ER, SYM, and ARD models, with indication of ordered traits.", 
           props = fp_text(italic = TRUE, font.size = 11)))

ft <- set_caption(ft, caption = caption_par)

# some cosmetic tweaks
ft <- bold(ft, part = "header")            # bold header row
ft <- autofit(ft)                          # adjust column widths
ft <- theme_vanilla(ft)                    # nice default theme

# show in viewer (RStudio will display it)
ft

########################################
logLik_info <- data.frame(
  Trait = c("Group Size", "Seasonal Breeding", "Mating System", "Reproductive Age",
            "Male Rank Bias", "Male Care", "Alloparenting", 
            "Male SSB", "Female SSB", "Transgender"),
  Observed_logLik = c("-81.368","-35.046", "-79.566","-37.285","-34.872","-25.436","-45.285","-19.539","-14.238","-15.261"),
  Randomized_logLik = c("-103.673","-45.702", "-105.225","-59.236","-47.813","-38.647","-83.647","-27.237","-28.014","-21.233"),
  p_value = c("9.999*10","8.991*10", "9.990*10","3.334*10","9.990*10","9.990*10","9.990*10","9.990*10","9.990*10","2.997*10"),
  Sup = c("-4","-3", "-4","-3","-4","-4","-4","-4","-4","-3"),
  Holm_adj = c("9.999*10","9.999*10", "9.999*10","9.999*10","9.999*10","9.999*10","9.999*10","9.999*10","9.999*10","9.999*10"),
  Sup2 = c("-3","-3", "-3","-3","-3","-3","-3","-3","-3","-3"),
  Is_significant = c("Yes","Yes", "Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes"),
  stringsAsFactors = FALSE
)

# build flextable
ft <- flextable(logLik_info)

ft <- compose(
  ft, j = "p_value", part = "body",
  value = as_paragraph(
    as_chunk(logLik_info$p_value),
    as_chunk(logLik_info$Sup, props = fp_text(vertical.align = "superscript"))
  )
)

ft <- compose(
  ft, j = "Holm_adj", part = "body",
  value = as_paragraph(
    as_chunk(logLik_info$Holm_adj),
    as_chunk(logLik_info$Sup2, props = fp_text(vertical.align = "superscript"))
  )
)

caption_par <- as_paragraph(
  as_chunk("Table 2. Observed log-likelihoods versus randomized log-likelihoods for each trait, including raw and Holm–Bonferroni corrected p-values. Significant traits are indicated in the table.", 
           props = fp_text(italic = TRUE, font.size = 11)))

ft <- set_caption(ft, caption = caption_par)

# some cosmetic tweaks
ft <- bold(ft, part = "header")            # bold header row
ft <- autofit(ft)                          # adjust column widths
ft <- theme_vanilla(ft)                    # nice default theme

# show in viewer (RStudio will display it)
ft <- delete_columns(ft, "Sup")
ft <- delete_columns(ft, "Sup2")

ft

##############################
logLik_info$Observed_logLik <- as.numeric(logLik_info$Observed_logLik)
logLik_info$Randomized_logLik <- as.numeric(logLik_info$Randomized_logLik)

df_long <- melt(logLik_info[,c("Trait","Observed_logLik","Randomized_logLik")],
                id.vars = "Trait", variable.name = "Type", value.name = "LogLik")

# Plot
plot <- ggplot(df_long, aes(x = Trait, y = LogLik, color = Type, group = Trait)) +
  geom_point(size = 3) +
  geom_line(color = "gray60", linetype = "dashed") +
  geom_text(data = logLik_info, aes(x = Trait, y = pmin(Observed_logLik) +5,
                                  label = ifelse(Is_significant=="Yes","*", "")),
            inherit.aes = FALSE, color = "goldenrod3", size = 6) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_markdown(size=10, hjust = 0.5)
  ) +
  labs(title = "Observed vs Randomized log-likelihood per trait",
       subtitle = "Stars indicate Holm-Bonferroni significant results",
       caption = "**Figure 3.** Observed versus randomized log-likelihoods for 10 behavioral traits. Stars indicate significant<br>phylogenetic signal after Holm-Bonferroni correction.",
       y = "Log-Likelihood", x = "", ) +
  scale_color_discrete(name = "", labels = c("Observed log-likelihood", "Mean Randomized log-likelihood"))
plot

