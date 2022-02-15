#' ARG Pipeline
#' Filename: inference functions.R
#' Written by: Lauren Castro
#' Created:   Feb  7, 2019
#' Edited on:  Oct 23, 2020
#' Description: Functions that help run the entire inference pipeline. Steps correspond to where they show up in 
#' "patient_inference_full_modified.R"

library(compiler)

## Step 4 Functions: Running the ARG simulation ####
run_trial_sensitivity <- function(activate.rate, recombo.rate, 
                                  rep.times, selection.strength, 
                                  selection.frequency, final.time,
                                  experiment.path, sample.times, sample.sizes, 
                                  starting.k,
                                  latent.pool) {  
  
  
  
  num.collapses = length(seq(from = 0, to = final.time, by = selection.frequency))
  
  sim_params <- default_params(final.time = final.time, 
                               collapse.length = 5,
                               num.collapses = num.collapses,
                               sample.sizes = sample.sizes, 
                               #sample.sizes = c(11,11,7,8,10,11,10,10,9), 
                               sample.times = rev(sample.times*30), 
                               #sample.times = rev(c(3,9,18,24,30,36,42,48,54)*30), 
                               re.rate = recombo.rate,
                               latent.size = latent.pool,
                               activate.rate = activate.rate,
                               pop.size = selection.strength,
                               starting.k = starting.k)

  
  ## Data frame indicating when bottlenecks will occur 
  population.boundaries = generate_pophistory_linear_current(num.collapses = sim_params$num.collapses,
                                                             pop.size = round(sim_params$pop.size), 
                                                             collapse.length = sim_params$collapse.length, 
                                                             final.time = sim_params$final.time,
                                                             selection.frequency = selection.frequency)
  # rep.times: defunct way of running 
  trials=seq(1:rep.times)
  trial <- sim_arg_full(sim_params, population.boundaries)
  return(trial)
}

## Step 5 Functions:  Get distance matrix of the binary tree of each residue  ####
return_arg <- function(trial, arg.path) {
  trial$arg$arg.list[[trial$arg$num.list]] <- trial$arg$current.arg
  arg = trial$arg$arg.list %>% bind_rows()%>%filter(time > 0)
  return(arg)
}


# Step 6 Functions: Calculate average distance matrix ####

# Below functions are no longer needed 
clean_binary_trees = function(binary.tree) {
  binary.tree = find_hanging_edges(binary.tree)
  interior.nodes = binary.tree[, .(count = .N), by = to.node][count == 1][order(to.node)][,to.node]
  max.node = max(as.numeric(as.character(binary.tree$to.node)))
  
  short.tree = max.node %in% interior.nodes 
  
  if(short.tree == TRUE) {
    interior.nodes = interior.nodes[-which(interior.nodes == max.node)]
    rows.to.delete = binary.tree[to.node == max.node, which = TRUE]
    set_column_NA(DT = binary.tree, del.indxs = rows.to.delete)
  }
  for(first.node in interior.nodes) {
    binary.tree = remove_interior_nodes(first.node = first.node, edge.list = binary.tree)
  }
  binary.tree = na.omit(binary.tree, cols = "from.node")
  return(binary.tree)
}
find_hanging_edges = function(edge.list) { 
  only.out = edge.list[which(!(from.node %in% to.node))][which(from.node.type == "I")][, from.node] 
  setkey(edge.list, from.node)
  edge.list <- edge.list[!only.out]
  return(edge.list)
}
remove_interior_nodes = function(first.node, edge.list) { 
  
  selected.rows = edge.list[from.node == first.node | to.node == first.node][order(desc(time))]
  new.row = data.table(from.node = selected.rows$from.node[1],
                       to.node = selected.rows$to.node[2],
                       edge.length = sum(selected.rows$edge.length),
                       from.node.type = selected.rows$from.node.type[1],
                       edge.type = selected.rows$edge.type[1],
                       recombo = selected.rows$recombo[1],
                       time = selected.rows$time[2],
                       edge.length.w=sum(selected.rows$edge.length.w))
  
  ## want to make these be NA and then update
  del.indx1 = edge.list[first.node, on ="from.node", which = TRUE]
  del.indx2 = edge.list[first.node, on = "to.node", which = TRUE]
  
  set_column_NA(DT = edge.list, del.indxs = c(del.indx1, del.indx2))
  edge.list[del.indx2, names(edge.list):=as.list(new.row)]
  return(edge.list[order(desc(time))])
}
full_clean_process = function(bt.list, bt.path, num.cores) {
  site.names = names(bt.list)
  bt.clean = purrr::map(bt.list, clean_binary_trees)
  names(bt.clean)=site.names
  return(bt.clean)
}

clean_binary_trees <- cmpfun(clean_binary_trees)
find_hanging_edges <- cmpfun(find_hanging_edges)
remove_interior_nodes <- cmpfun(remove_interior_nodes)
full_clean_process <- cmpfun(full_clean_process)

full_dm_trial <- function(bt.clean, num.samples, sample.sizes, bt.path, dm.path, num.cores, starting.k) {
  pwd = bt.clean %>% purrr::map(calculate_dm_short, num.samples = num.samples, sample.sizes = sample.sizes, starting.k=starting.k)
  names(pwd)=names(bt.clean)
  return(pwd)
  
}
calculate_dm_short <- function(edge.list, num.samples, sample.sizes, starting.k) { 

  #create node id 
  root = max(as.numeric(as.character(edge.list$to.node)))
  tip.nodes = edge.list[from.node.type == "T"][,.(from.node)]
  
  sample.times = seq(1:(num.samples+1))
  samp.sizes = c(starting.k, sample.sizes)
  time.samples = rep(x=sample.times, samp.sizes)
  
  tip.nodes[,from.node := as.numeric(as.character(from.node))] %>% 
    .[order(from.node)] %>% 
    .[, time.sample:=time.samples] -> tip.nodes
 
  colnames(edge.list)[8] = "weight"
  
  #### Need to assign Time samples 
  g2=graph.data.frame(edge.list, directed = FALSE)
  
  dm = shortest.paths(g2, weights = edge.list$weight)
  dm %<>% as.data.table(keep.rownames = TRUE) 
  pwd = melt(dm, variable.name = "to", value.name = "distance", id.vars = "rn")
  setnames(pwd, old = "rn", "rowname")
  pwd.short  = return_pwd_tips(pwd, tip.nodes = tip.nodes)
  return(list(pwd = pwd.short, tip.nodes = tip.nodes))
}


# Step 7 Functions: Calculate Average Tree  ####
combine_pwd_tips = function(pwd.list.entry){
  pwd = pwd.list.entry$pwd
  tip.nodes = pwd.list.entry$tip.nodes
  
  pwd %>%
    .[,rowname:=as.numeric(as.character(rowname))] %>%
    .[,to:=as.numeric(as.character(to))]
  
  pwd = merge(pwd, tip.nodes, by.x="rowname", by.y="from.node", all.x=T)
  setnames(pwd, old = "time.sample", new = "from.time.sample")
  pwd = merge(pwd, tip.nodes, by.x = "to", by.y  = "from.node", all.x = T)
  setnames(pwd, old = "time.sample", new = "to.time.sample")
  return(pwd)
}
load_pwd = function(dm, dm.path) {
  
  full.r.list = purrr::map(.x = dm, combine_pwd_tips)
  ## combines all sites and all simulations for parameter set
  full.r = rbindlist(l=full.r.list, idcol = "r.site") 
  return(full.r)
}
calculate_avg_dm = function(dm.entry) { 
  dm.entry[,.(avg.pwd = mean(distance)), by = .(rowname, to)] %>% 
    dcast(., rowname~to, value.var = "avg.pwd") -> avg.matrix
  tip.names=avg.matrix$rowname
  #avg.matrix[,tip.names:= rowname]
  avg.matrix[,rowname:=NULL]
  avg.matrix = as.matrix(avg.matrix)
  dimnames(avg.matrix) = list(tip.names, tip.names)
  return(avg.matrix)
}
format_dm_long <- function(dm.entry) { 
  dm.entry %>% 
    as.data.frame %>%
    rownames_to_column %>% 
    gather(key = tip2, value = distance, -rowname) 
}
get_avg_tree <- function(dm, dm.path, activate.rate, sample.sizes) {

  ## Step 1: Load
  pwd.full = load_pwd(dm = dm, dm.path) ### at this point -- sample time one refers to the latest grabbed 
  ## Step 2: Calculate the average dm
  pwd.avg = calculate_avg_dm(pwd.full) 
  
  ## Step 3: Plot and save tree  
  samples = colnames(pwd.avg)
  tips.length = length(samples)
  sample.sizes = sample.sizes
  time.samples = samples %>% as.data.frame %>% 
    mutate(time = rep(x = seq(1:length(sample.sizes)), sample.sizes),
           sample.time = rep(x = rev(1:length(sample.sizes)), sample.sizes))
  
  ## Calculate tree 
  tr <- fastme.bal(pwd.avg)
  ## Rename tips so that they correspond to the time samples where 1 is the first sampled time in forward time 
  tip.order = tr$tip.label
  new.tips = tip.order %>% as.data.frame %>% left_join(time.samples, by = '.') %>% pull(sample.time)
  tr$tip.label = new.tips
  
  anc <- ape::getMRCA(phy=tr, tip = which(tr$tip.label == 1))
  if(anc == (length(tr$tip.label)+1)) {
    rooted.tree = root(phy = tr, node = anc)
  } else{
    rooted.tree = ape::root(phy = tr, node =  anc, resolve.root = TRUE)
  }
  return(rooted.tree)
}
get_avg_dm <- function(dm, dm.path, sample.sizes, sample.times) {
  pwd.full = load_pwd(dm=dm, dm.path = dm.path)
  
  ## Step 2: Calculate the average dm
  pwd.avg = calculate_avg_dm(pwd.full)
  rm(pwd.full)
  
  ## Step 3: Format into long version with the time samples as columns
  samples = colnames(pwd.avg)
  tips.length = length(samples)
  
  tips.length = length(samples)
  sample.sizes = rev(sample.sizes)
  time.samples = samples %>% as.data.frame %>% 
    mutate(#time = rep(x = seq(1:length(sample.sizes)), sample.sizes),
      sample.time = rep(x = rev(1:length(sample.sizes)), sample.sizes))
  
  colnames(time.samples)[1] = "tip"
  
  pwd.l = format_dm_long(pwd.avg) 
  #pwd.l = as.data.table(pwd.l)
  
  sample.months = sample.times
  sample.months = as.data.frame(cbind(sample.days = sample.months*30, time.point = seq(1:length(sample.times)), sample.sizes = rev(sample.sizes))) 
  
  pwd.l$rowname = factor(pwd.l$rowname, levels =  rev(sort(unique(as.numeric(as.character(pwd.l$rowname))))))
  pwd.l$tip2 = factor(pwd.l$tip2, levels =  rev(sort(unique(as.numeric(as.character(pwd.l$tip2))))))
  
  tip.matrix = pwd.l %>% select(rowname, tip2, distance) %>% spread(key = tip2, value = distance)
  rownames(tip.matrix) = c()  
  tip.matrix = tip.matrix %>% column_to_rownames(var = "rowname") 
  
  day.order = rep(sample.months$sample.days, times = sample.months$sample.sizes)
  
  column.names = colnames(tip.matrix)
  tip.names = as.data.frame(cbind(column.names, day.order)) %>% mutate(name = paste0(column.names, "_", day.order)) %>% pull(name)
  
  colnames(tip.matrix) = tip.names
  rownames(tip.matrix) = tip.names
  return(tip.matrix)
}
get_avg_tree_modified <- function(avg.dm, sample.sizes) {

  ## Step 1: Plot and save tree 
  samples = seq(1:ncol(avg.dm))
  tips.length = length(samples)
  #sample.sizes = c(11,11,11,7,8,10,11,10,10,9) # backwards
  time.samples = samples %>% as.data.frame %>% 
    mutate(time = rep(x = seq(1:length(sample.sizes)), sample.sizes),
           sample.time = rep(x = rev(1:length(sample.sizes)), sample.sizes))
  
  ## Calculate tree 
  tr <- fastme.bal(avg.dm)
  ## Rename tips so that they correspond to the time samples where 1 is the first sampled time in forward time 
  tip.order = as.numeric(as.character(tr$tip.label))
  new.tips = tip.order %>% as.data.frame %>% left_join(time.samples, by = '.') %>% pull(sample.time)
  tr$tip.label = new.tips
  
  anc <- ape::getMRCA(phy=tr, tip = which(tr$tip.label == 1))
  if(anc == (length(tr$tip.label)+1)) {
    rooted.tree = root(phy = tr, node = anc)
  } else{
    rooted.tree = ape::root(phy = tr, node =  anc, resolve.root = TRUE)
  }
  rooted.tree = ladderize(rooted.tree)
  return(rooted.tree)
}

# Step 8: Get Metrics and Save Score ####
get_tree_metrics <- function(tree, sample.times) { 
  
  ### sackin.index 
  sackin.index = sackin.index(tree, norm = T)
  
  ### Ratio of internal to external 
  inner = mean(tree$edge.length[tree$edge[,2] > Ntip(tree)])
  outer = mean(tree$edge.length[tree$edge[,2] <= Ntip(tree)])
  ratio = outer/inner
  
  ### Lineages surviving through times 
  time.names <- sort(as.numeric(as.character(unique(tree$tip.label))))
  
  taxa.sampling.time = as.numeric(as.character(tree$tip.label)) 
  taxa.sampling.times = sort(unique(taxa.sampling.time))
  
  n.samples <- length(taxa.sampling.times)
  
  lineages.per.sample <- c()
  for(i in 2:n.samples){
    current.sample <- which(taxa.sampling.time == taxa.sampling.times[i])
    current.MRCA <- getMRCA(tree, current.sample)
    current.tr <- extract.clade(tree, node=current.MRCA)
    
    current.tr <- drop.tip(phy=current.tr, tip = which(as.numeric(as.character(current.tr$tip.label)) > taxa.sampling.times[i]))
    
    current.runs <- grep(time.names[i], current.tr$tip.label)
    current.runs[length(current.runs)+1]<-9999
    lineages.per.sample[i] <- length(which(diff(current.runs)>1))
  }
  ltt = data.frame(cbind(lineages.per.sample, time = taxa.sampling.times, sample.times))
  
  mean.lin.time = mean(lineages.per.sample[-1]/(table(taxa.sampling.time)[-1])) 
  
  
  return(list(sackin = sackin.index, e.i.ratio = ratio, ltt = ltt, mean.lin.time = mean.lin.time))
}
time.sort = function(x, n){
  ret = c()
  for(i in unique(n)){
    for(j in unique(n)){
      if(j>=i){
        # print(i)
        a = unlist(c(x[which(n==i), which(n==j)]))
        ret = c(ret, sort(a))
      }
    }
  } 
  unname(ret)
}
time.sort.cv = function(x, n){
  
  ret = matrix(, ncol = 3)
  for(i in unique(n)){
    for(j in unique(n)){
      if(j>=i){
        # print(i)
        a = unlist(c(x[which(n==i), which(n==j)]))
        edges = matrix(nrow = length(a), ncol = 3)
        edges[,1] = i
        edges[,2] = j
        edges[,3] = sort(a)
        
        l <- list(ret, edges) 
        ret = do.call(rbind, l)
      }
    }
  } 
  unname(ret)
}
#y = time.sort(y, as.numeric(rownames(y)))
#s = time.sort(s, tmp)
fun = function(lam, y, s){
  val = sum((y-(s*lam))^2)
  #print(val)
  val  
}
find_pairwise_objective <- function(y, s) {
  
  tmp = as.numeric(sapply(strsplit(rownames(s), "_"), "[[", 2))#/30
  y.sort = time.sort(y, as.numeric(rownames(y)))
  s.sort = time.sort(s, tmp)
  
  optimal = optimize(fun, c(0,1), y=y.sort, s=s.sort,  tol = 1e-15)  
  return(list(objective = optimal$objective, lambda = optimal$minimum))
}
find_cv <- function(avg.dm){
  
  tmp = as.numeric(sapply(strsplit(rownames(avg.dm), "_"), "[[", 2))#/30
  s.sort = time.sort.cv(avg.dm, tmp)
  colnames(s.sort) = c("time1", "time2", "distance")
  s.sort = data.table(s.sort)
  cv = s.sort[time1==time2, .(mean = mean(distance), sd= sd(distance)), by = .(time1)][
    , .(CV = sd/mean), by = .(time1)]
  return(cv)
}
