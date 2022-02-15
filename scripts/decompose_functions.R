#; ARG Pipeline
#' Filename: decompose_functions.R
#' Written by: Ethan Romero Severson & Lauren Castro
#' Created:   August 11, 2019
#' Edited on: October 22, 2019

#  Description                   ----
#' How this script works (high-level):
#' Functions to efficiently decompose the ARG into binary trees at each residue

library(data.table)

get.arg = function(arg) {
  arg$arg$arg.list[[arg$arg$num.list]] <- arg$arg$current.arg
  arg = rbindlist(arg$arg$arg.list) 
  arg[time >0]
}
get.lineage = function(from, to, ln, tip.ind, res.ind, bp, root){
  if(tip.ind%%1!=0)return("Not a tip")
  #if(nt[from==tip.ind]!="T")return("Not a tip")
  cur.ind = tip.ind
  anc = ctimes = vector(mode = "numeric", length = length(unique(bp)))
  
  ctime = 0
  running = T
  counter = 1
  
  while(running) {
    cur.index.position = which(from%in%cur.ind)
    
    #parents = to[from %in% cur.ind]  
    parents = to[cur.index.position]  
    parents.length = length(parents)
    
    if(parents.length==1){ #is C type
      ctime = ctime + ln[cur.index.position]
      #anc = c(anc, parents)
      anc[counter] = parents
      #ctimes = c(ctimes, ctime)
      ctimes[counter] = ctime
      cur.ind = parents
      counter = counter + 1
    }
    if(parents.length==2){# is R type
      ctime = ctime + ln[cur.index.position][1]
      bpoint = bp[cur.index.position][1]
      cur.ind = ifelse(res.ind <= bpoint, parents[1], parents[2])
    }
    if(parents.length>2){
      
      if(parents.length%%2!=0){# sequence ends in C
        ctime = ctime + sum(ln[cur.index.position][c(seq(1,parents.length-1,2), parents.length)]) # this may be off 
        #anc = c(anc, parents[parents.length])
        anc[counter] = parents[parents.length]
        cur.ind = parents[parents.length]
        #ctimes = c(ctimes, ctime)
        ctimes[counter] = ctime
        counter = counter + 1
      }else{ # sequence ends in a RR
        ctime = ctime + sum(ln[cur.index.position][c(seq(1,parents.length-2,2), parents.length-1)])
        bpoint = bp[cur.index.position][parents.length]
        cur.ind = ifelse(res.ind <= bpoint, parents[parents.length-1], parents[parents.length])
      }
    }
    if(as.numeric(cur.ind)==root)running = F
  }
 
  list(anc[1:counter-1],ctimes[1:counter-1])
}
# get.lineage(from,to,ln,5,100)
get.residue.dmat = function(arg, res.ind) {
  arg = get.arg(arg)
  to = as.numeric(arg$to.node)
  root = max(to, na.rm=T)
  end = max(which(to==root))
  to = to[1:end]
  
  from = as.numeric(arg$from.node)[1:end]
  ln = as.numeric(arg$edge.length)[1:end]
  et = arg$edge.type[1:end]
  #node info
  bp = as.numeric(arg$recombo)[1:end]
  nt = arg$from.node.type[1:end]
  
  tips = unique(from[nt=="T"])
  dat = vector(mode = "list", length = length(tips))

  for(i in 1:length(tips))  {
     dat[[i]] = get.lineage(from,to,ln,tips[i],res.ind, bp = bp, root)
     names(dat)[i] = as.character(tips[i])
  }

  ret = matrix(0, nrow=length(tips), ncol=length(tips))
  #
  for(i in 1:length(tips)) {
    anc.i = dat[[as.character(tips[i])]][[1]] # anc.i -- nodes in i's path
    ct.i = dat[[as.character(tips[i])]][[2]]  # anc. i --- node height's for i
    for(j in 1:length(tips)){
      if(j>i){
        anc.j = dat[[as.character(tips[j])]][[1]] # anc j -- nodes in j's path
        ct.j = dat[[as.character(tips[j])]][[2]] #  anc j -- node height's in j's path
        mrca = intersect(anc.j, anc.i)[1] ## find first shared node
        tmrca = sum(ct.i[which(anc.i==mrca)], ct.j[which(anc.j==mrca)])
        ret[i,j] = tmrca
        ret[j,i] = tmrca
      }
    }
  }
 
  return(ret)
}
