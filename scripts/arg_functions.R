#' ARG Functions
#' Filename: arg_functions.R
#' Written by: Lauren Castro
#' Created:   Feb  7, 2019
#' Edited on:  Oct 23, 2020
#' Description: 
#' Functions that control the ARG simulation and 


suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(data.table))

generate_pophistory_constant       <- function(num.collapses, pop.sizes, collapse.length, final.time) {
  population.boundaries =  data.frame(start = integer(), end = integer(),
                                      pop.realm = character(), pop.size = integer())
  population.boundaries$pop.realm  = factor(population.boundaries$pop.realm, levels = c("low", "high"))
  collapse.moments = data.frame(start = sample(1:final.time, size = num.collapses)) # don't have to worry about overlapping if
  collapse.moments = collapse.moments %>% arrange(desc(start))
  
  collapse.moments %>%
    mutate(valid = start-lead(start, n = 1) > collapse.length) -> collapse.valid
  
  while(length(which(collapse.valid$valid==FALSE)) > 0) {
    collapse.moments = data.frame(start = sample(1:final.time, size = num.collapses))
    collapse.moments = collapse.moments %>% arrange(desc(start))
    
    collapse.moments %>%
      mutate(valid = start-lead(start, n = 1) > collapse.length) -> collapse.valid
  }
  
  population.boundaries[1,] = data.frame(start = as.integer(final.time),  end = collapse.moments$start[1],
                                         pop.realm = "high",  pop.size = as.integer(pop.sizes$high))
  
  
  for(i in 1:nrow(collapse.moments)) {
    population.boundaries = rbind(population.boundaries, data.frame(start = collapse.moments$start[i], end = collapse.moments$start[i] - collapse.length,
                                                                    pop.realm = "low", pop.size = as.integer(pop.sizes$low)))
    population.boundaries = rbind(population.boundaries, data.frame(start = population.boundaries$end[nrow(population.boundaries)], end = collapse.moments$start[i+1],
                                                                    pop.realm = "high", pop.size = as.integer(pop.sizes$high)))
  }
  population.boundaries$end[which(is.na(population.boundaries$end))] = 0
  return(population.boundaries)
}
generate_pophistory_linear_current <- function(num.collapses, pop.size, collapse.length, 
                                               final.time, selection.frequency) {
  population.boundaries =  data.frame(start = integer(length = num.collapses), end = integer(length = num.collapses), pop.size = integer(length = num.collapses))
  collapse.moments = data.frame(start = round(seq(from = selection.frequency , to = final.time-selection.frequency, length.out = num.collapses)))
  for(i in 1:nrow(collapse.moments)) {
    population.boundaries[i,] = data.frame(start = collapse.moments$start[i],  end = collapse.moments$start[i]-collapse.length, 
                                           pop.size = pop.size)
  }
  population.boundaries %>% 
    arrange(desc(start))
}
generate_pophistory_linear_old     <- function(num.collapses, pop.sizes, collapse.length, final.time) {
  ## Currently uses a weird hack that doesn't allow a sample time to begin to end with where the latent/boundary is 
  population.boundaries =  data.frame(start = integer(length = num.collapses), end = integer(length = num.collapses), pop.size = integer(length = num.collapses))
  collapse.moments = data.frame(start = seq(from = 90, to = final.time, by = 180))
  for(i in 1:nrow(collapse.moments)) {
    population.boundaries[i,] = data.frame(start =  collapse.moments$start[i],  end = collapse.moments$start[i]-collapse.length, 
                                           pop.size = pop.sizes)
  }
  population.boundaries %>% 
    arrange(desc(start))
}
select_possible_events             <- function(nxt.event.type) {
  if(nxt.event.type == "R") {
    possible.events = c("R", "C", "L")
  } else if(nxt.event.type == "L") {
    possible.events = c("A")
  } else if(nxt.event.type == "C") {
    possible.events = c("R", "C", "L")
  } else {
    possible.events = c("R", "C", "L")
  }
  return(possible.events)
} # DT checked 



## Rate functions for events in ARG simulation
exp_coalesent_pop <- function(k, pop.size)  { #Recombo rate per day 
  # Exponential coalescent waiting time for a population with constant size 
  # k -- extant lineages
  coal.rate =(k*(k-1))/(2*pop.size)
  coal.wait  = rexp(n = 1, rate = coal.rate)
  if(is.na(coal.wait)) {
    browser()
  }
  return(coal.wait)
} # Contains defaults 
inv_cumf_activate <- function(u, time,k_l,lambda_l, lambda_a, alpha, beta = 1) {
  # Inverse cumulative function for time-varying activation rate
  # u--random uniform
  # time -- current time
  # k -- extant lineages
  # alpha -- starting population size
  # beta -- growth rate
  (1-(1-u)^((beta*lambda_l)/(k_l*lambda_a)))*(alpha + beta*time)*(1/beta)
}
inv_cumf_coal     <- function(u, time, k, alpha, beta = 1) {
  # Inverse cumulative function for time-varying coalescent
  ## u -- random uniform
  ## time -- current time
  ## k -- extant lineages
  ## alpha -- starting population size 
  ## beta -- growth rate 
  
  (1-(1-u)^(beta/choose(k,2))) * (alpha + beta*time)*(1/beta)
  
} # Contains defaults alpha and beta 
inv_cumf_latent   <- function(u, time, keff, beta=1, rho, alpha) {
  # Inverse cumulative function for time-varying latent rate
  # u--random uniform
  # time -- current time
  # k -- extant lineages
  # alpha -- starting population size
  # beta -- growth rate
  (1-(1-u)^(beta/(keff*rho)))*(alpha + beta*time)*(1/beta)
}
exp_recombo       <- function(k, re.rate)  { #Recombo rate per day 
  # Inverse cumulative function for time to next recombination event
  # k -- extant lineages
  # re.rate -- recombo rate in days 
  rexp(n = 1, rate = k*re.rate)
  #if(is.nan(rate)) break()
} # Contains defaults 
exp_latent        <- function(keff, latent.rate) { 
  # keff -- extant active lineages
  # re.rate -- recombo rate in days 
  rexp(n = 1, rate = keff*latent.rate)
}
exp_activate      <- function(klat = (k-keff), activate.rate) {
  # klat -- extant latent lineages
  # activate.rate -- activate rate in days 
  rexp(n = 1, rate = klat*activate.rate)
}
exp_activate2     <- function(klat = kl, rho, N_l) {
  #klat -- extant latent lineages
  #rho -- lambda/activate rate
  #N_L -- population size 
  rexp(n=1, rate = (klat/N_l)*rho)
} 


## Functions for elements within simulation 
find_edge_length_m <- function(host.arg.list, list.entry, host.arg, time.index, time.event, selected.lineage, possible.events) {
  #Is the selected.lineage in currrent ? 
  if(selected.lineage == "1e5.2") {
    ##browser()
  }
  host.match = host.arg[to.node == as.numeric(selected.lineage)][edge.type %in% possible.events]
  # host.arg %>%
  #   filter(to.node == as.numeric(selected.lineage)) -> total.match
  # total.match %>%
  #   filter(edge.type %in% possible.events)-> host.match
  
  if(nrow(host.match) > 0) {
    coal.desc = host.match[order(-time)][.N] 
    edge.length = (coal.desc$time[1]-time.index + time.event)
  } else {
    while(list.entry > 1) { # While there are more than one list entry in arg; was >=
      older.arg = host.arg.list[[list.entry-1]]
      host.match = older.arg[to.node == as.numeric(selected.lineage)] 
      if(nrow(host.match)>0) {
        coal.desc = host.match[order(-time)][.N]
        edge.length = (coal.desc$time[1]-time.index + time.event)
        list.entry = 0 # break out of loop
      } else{list.entry = list.entry - 1} # If there is no match, try the previous list 
    }
  }
  return(edge.length)
}
determine_edge     <- function(selected.lineage, host.arg, nxt.event, edge.type, time, sample.time, tips, possible.events, tip.df) {
  # selected.lineage -- the from node
  # host.arg.list/host.arg/num. entry for finding previous  entry
  # nxt event -- the time of the event 
  # edge type -- recombination/coalescent/initiate latency/activate latent edge
  # time -- current time
  # sample.time -- when the lineages were first sampled 
  edge.type = edge.type
  prev.events = nrow(host.arg$current.arg["to.node" == selected.lineage])
  #if(selected.lineage %in% tips & (length(which(host.arg$current.arg$to.node == selected.lineage)) == 0)) {
  if(selected.lineage %in% tips & prev.events == 0) {
    from.node.type = "T"
    ##browser()
    if(nrow(tip.df[which(tip.df$tip.number==selected.lineage),]) > 0) { 
      sampled.time = as.numeric(as.character(tip.df[which(tip.df$tip.number==selected.lineage), "sample.time"]))
    }else{
        sampled.time = as.numeric(as.character(sample.time))
      }
    edge.length = sampled.time-time + nxt.event
    if(edge.length < 0) browser()#break() #browser()
  } else { 
    from.node.type = "I"
    edge.length = find_edge_length_m(host.arg.list = host.arg$arg.list, list.entry = host.arg$num.list, host.arg = host.arg$current.arg, 
                                     time.index = time, time.event = nxt.event, selected.lineage = selected.lineage, possible.events)
    if(edge.length < 0) browser() # break()#browser()
  }
  return(edge.data = list(edge.type = edge.type, from.node.type=from.node.type, edge.length = edge.length))
}
find_recombo_site  <- function(node, host.arg, list.entry, host.arg.list) {
  
  host.arg[to.node == as.numeric(node)] -> host.match
  if(nrow(host.match) > 0) {
    recombo.site = host.match$recombo[1]
  } else {
    while(list.entry >=1) { # While there are more than one list entry in arg 
      older.arg = host.arg.list[[list.entry-1]]
      host.match = older.arg[to.node == as.numeric(node)]
      if(nrow(host.match) > 0) {
        recombo.site = host.match$recombo[1]
        list.entry = 0 # break otu of loop
      } else{list.entry = list.entry - 1} # If there is no match, try the previous list 
    }
  }
  return(recombo.site)
}

create_arg_entry    <- function(from.node, to.node, edge.length, from.node.type, edge.type, recombo, time, nxt.event) {
  new.entry = data.table(from.node = from.node, to.node = to.node, edge.length = edge.length, from.node.type=from.node.type,
                         edge.type = edge.type, recombo = recombo, time = time-nxt.event)
  new.entry %>%
    mutate_if(is.factor, as.character) %>% as.data.table
}
check_lineage_tip_a <- function(starting.k, host.a.node, selected.lineage) {
  ((host.a.node) <= as.numeric(as.character(selected.lineage))) & (as.numeric(as.character(selected.lineage)) <= (host.a.node + starting.k))
}
check_lineage_tip_b <- function(starting.k, selected.lineage) {
  (1 <= as.numeric(as.character(selected.lineage)) & (as.numeric(as.character(selected.lineage)) <= starting.k))
}
add_arg_entry       <- function(host.arg, new.entry, N) {
  if(host.arg$arg.entries < (N-1)) {
    #host.arg$current.arg[host.arg$arg.entries,] = new.entry
    set(host.arg$current.arg, as.integer(host.arg$arg.entries), names(host.arg$current.arg), as.list(new.entry))
    host.arg$arg.entries = host.arg$arg.entries  + 1
  } else{
    host.arg$arg.list[[host.arg$num.list]] <- host.arg$current.arg
    host.arg$current.arg <- initialize_arg(N)
    host.arg$num.list = host.arg$num.list + 1 
    host.arg$arg.entries=1
    set(host.arg$current.arg, as.integer(host.arg$arg.entries), names(host.arg$current.arg), as.list(new.entry))
    host.arg$arg.entries = host.arg$arg.entries  + 1
    if(host.arg$num.list >= 21 && (host.arg$num.list %% 10)==1 ) {
      extra.list = vector("list", 20)
      host.arg$arg.list = do.call(c, list(host.arg$arg.list, extra.list))
      }
  }
  return(host.arg)
}
add_track_entry     <- function(track.list, new.entry, N) {
  if(track.list$track.entries < (N-1)) {
    set(track.list$track.df, as.integer(track.list$track.entries), names(track.list$track.df), as.list(new.entry))
    track.list$track.entries = track.list$track.entries  + 1
  } else{
    track.list$track.list[[track.list$num.list]] <- track.list$track.df
    track.list$track.df <- initialize_track_df(N)
    track.list$num.list = track.list$num.list + 1 
    track.list$track.entries=1
    set(track.list$track.df, as.integer(track.list$track.entries), names(track.list$track.df), as.list(new.entry))
    track.list$track.entries = track.list$track.entries + 1
    if(track.list$num.list >= 21 && (track.list$num.list%%10) == 1) {
      extra.list = vector("list", 10)
      track.list$track.list = do.call(c, list(track.list$track.list, extra.list))
    }
  }
  return(track.list)
}

initialize_track_list <- function(N) {
  track.list = list(
    track.df = initialize_track_df(N),
    track.entries = 1,
    num.list = 1,
    track.list = vector("list", 20)
  )
}
initialize_track_df   <- function(N) {
  track <- data.table(
    lineages.total = numeric(N),
    lineages.effective = numeric(N),
    event = character(N),
    time.lapse = numeric(N),
    sim.time = numeric(N),
    pop.size = numeric(N)
  )
  track %>% map_if(is.factor, as.character) %>% as.data.table
}
combine_track         <- function(host.a.track, host.b.track, source.pop.track) {
  list.arg = list(host.a.track, host.b.track, source.pop.track)
  arg.source = c("a", "b", "source")
  
  nrow.list = map(list.arg, nrow)
  
  list.arg[which(nrow.list>0)] %>%
    bind_rows()  -> arg.df
  
  arg.df %>%
    mutate(source = rep(arg.source, nrow.list))
}
print_statement       <- function(k, keff, time) {
  print(paste0("Time: ", time))
  print(paste0("Effective Lineages: ", keff))
  print(paste0("Total Lineages: ", k))
}

# Functions for organizing ARG 
initialize_arg     <- function(N) {
  ## Initialize data frame of N size for ARG
  arg <- data.table(
    from.node = character(N),
    to.node = character(N),
    edge.length = numeric(N),
    from.node.type = character(N),
    edge.type = character(N),
    recombo = integer(N),
    time = numeric(N)
  )
  arg%>%map_if(is.factor, as.character) %>% as.data.table
}
initialize_arglist <- function(N) {
  host.arg = list(
    current.arg = initialize_arg(N),
    arg.entries = 1,
    num.list = 1,
    arg.list = vector("list", 20)
  )
}
initialize_time    <- function(N) {
  time.track <- data.table( 
    coal = numeric(N),
    activate = numeric(N),
    latent = numeric(N),
    recombo = numeric(N),
    k = numeric(N),
    keff = numeric(N),
    k_l = numeric(N),
    time = numeric(N)
  )
}
add_latent_arg     <- function(lineage, shost.arg, new.entry, N) {
  
}



sim_arg_initial_full <- function(starting.k, starting.node, final.time, sample.time, nxt.sample.time, alpha, re.rate, latent.size, activate.rate,
                                 population.boundaries, tip.df) {

  # Simulates coalescent/recombination/latency within derived populations
  # Time starts when the host population is sampled 
  # Initialize node list and node number for book keeping
  N=10000
  track.list = initialize_track_list(N)
  host.arg = initialize_arglist(N)  #time.track = initialize_time(N)
  time.entry = 1
  
  extant.active.list = vector(mode = "character", length = 500) # Create active list size 
  extant.active.list[1:starting.k] = as.character(seq(from = starting.node, to = (starting.node+starting.k-1), by = 1)) # Populate starting
  extant.latent.list = vector(mode = "character", length = 500) # Create latent list size 
  
  tips = as.character(as.numeric(seq(from = starting.node, to = starting.node + starting.k-1, by = 1)))
  node.number = starting.node+starting.k-1 # Keeps track of how many nodes are in the ARG  
  k = starting.k # k is the number of extant lineages 
  keff = k # Effective K (only the active)
  time = sample.time 
  
  pop.realm = 0
  sweep = FALSE
  nxt.boundary = population.boundaries$start[1]
  #print(paste0("Time;", time, "---Next Sample Time: ", nxt.sample.time))
  latent.phase = 2
  # Time boundary between sample and next sample
  while (time > nxt.sample.time) { 
    
    #### Generate times to next event 
    ## Coalescent dependent on population size  
    if(sweep == FALSE) {
      if(keff > 1) nxt.coal = inv_cumf_coal(u=runif(1), time = time, alpha=alpha, k=keff) else nxt.coal = 10000
      pop.size = alpha + 1*time
    } else {
      if(keff > 1 ) nxt.coal = exp_coalesent_pop(k=keff, pop.size = population.boundaries$pop.size[pop.realm]) else nxt.coal = 10000
    }
    
    if(keff >0) nxt.recombo = exp_recombo(k=keff, re.rate = re.rate) else nxt.recombo = 10e12
    if(latent.phase == 2 & keff > 0) nxt.latent = inv_cumf_latent(runif(1), time = time, keff = keff, rho = activate.rate, alpha = alpha) else nxt.latent = 10e12
    if(k-keff > 0) nxt.activate = exp_activate2(klat = (k-keff), rho = activate.rate, N_l = latent.size) else nxt.activate = 10e12
    
    waiting.times = c(nxt.coal, nxt.recombo, nxt.latent, nxt.activate)
    waiting.times = waiting.times[which(!is.na(waiting.times))]
    
    ## Determine if next boundary switch is because of population sweep or sampling 
    if(latent.phase == 2) nxt.switch = max(nxt.boundary, nxt.sample.time, 21) else nxt.switch = max(nxt.boundary, nxt.sample.time)
    if(time - min(waiting.times) > nxt.switch) { # Will the next event take place before the time boundary for population size ? 
      if(min(waiting.times) == nxt.recombo) {
        ###### Recombination Event 
        selected.lineage = sample(extant.active.list[1:keff], size = 1) # select random lineage
        node.number = node.number + 1  # designate new node number 
        possible.events = select_possible_events(nxt.event.type = "R")
        graph.edge = determine_edge(selected.lineage, host.arg, nxt.recombo, edge.type = "R", time, sample.time, tips, possible.events, tip.df=tip.df)
        ## Creating new entry into the ARG 
        recombo.site = sample(1:700, size = 1)
        new.entry1 = create_arg_entry(from.node = selected.lineage, to.node =  paste0(node.number, ".1"), edge.length = graph.edge$edge.length,
                                      nxt.event = nxt.recombo, from.node.type = graph.edge$from.node.type, edge.type = graph.edge$edge.type, recombo = recombo.site, time = time)
        new.entry2 = create_arg_entry(from.node = selected.lineage, to.node = paste0(node.number, ".2"), edge.length = graph.edge$edge.length,
                                      nxt.event = nxt.recombo, from.node.type = graph.edge$from.node.type, edge.type = graph.edge$edge.type, recombo = recombo.site, time = time)
        
        host.arg = add_arg_entry(host.arg, new.entry1, N) 
        host.arg = add_arg_entry(host.arg, new.entry2, N)
        
        ## Add track entry
        track.entry=data.table(lineages.total = k, lineages.effective = keff, event = as.factor("R"), time.lapse=nxt.recombo, sim.time = time, pop.size = pop.size) %>%
          mutate_if(is.factor, as.character) %>% as.data.table
        track.list = add_track_entry(track.list, track.entry, N)
        
        # Recombination Event/Bookeeping 
        extant.active.list = extant.active.list[-match(selected.lineage, extant.active.list)] # Remove the recombinant 
        extant.active.list[c(keff,keff+1)] = c(paste0(node.number,".1"), paste0(node.number,".2"))
        keff = keff + 1 # Recombination adds a lineage to active lineages
        k = k + 1 #  Recombination adds a lineage to total 
        time = time-nxt.recombo  # Move backward in time 
      } else if(min(waiting.times) == nxt.coal) {
        ###### Coalescent Event
        selected.lineages = sample(extant.active.list[1:keff], size = 2, replace = FALSE) # Two lineages coalesce
        node.number = node.number + 1 # designate new node number 
        
        #First Entry -- check if coalescing lineage is a tip
        possible.events = select_possible_events(nxt.event.type = "C")
        graph.edge.1 = determine_edge(selected.lineage=selected.lineages[1], host.arg, nxt.coal,
                                      edge.type = "C", time, sample.time, tips, possible.events, tip.df)
        new.entry.1 = create_arg_entry(from.node = selected.lineages[1], to.node = as.factor(node.number), edge.length = graph.edge.1$edge.length,
                                       nxt.event = nxt.coal, from.node.type = graph.edge.1$from.node.type, edge.type = graph.edge.1$edge.type, recombo = NA, time = time)
        
        #Second Entry -- check if a coalescing lineage is a tip
        graph.edge.2 = determine_edge(selected.lineage=selected.lineages[2], host.arg, nxt.coal,
                                      edge.type = "C", time, sample.time, tips, possible.events, tip.df)
        new.entry.2 = create_arg_entry(from.node = selected.lineages[2], to.node = as.factor(node.number), edge.length = graph.edge.2$edge.length,
                                       nxt.event = nxt.coal, from.node.type = graph.edge.2$from.node.type, edge.type = graph.edge.2$edge.type, recombo = NA, time = time)
        
        host.arg = add_arg_entry(host.arg, new.entry.1, N)
        host.arg = add_arg_entry(host.arg, new.entry.2, N)
        
        ## Add track entry
        track.entry=data.table(lineages.total = k, lineages.effective = keff, event = "C", time.lapse=nxt.coal, sim.time = time,pop.size = pop.size) %>% 
          mutate_if(is.factor, as.character) %>% as.data.table
        track.list = add_track_entry(track.list, track.entry, N)
        
        # Coalescent Event/Bookkeeping 
        extant.active.list = extant.active.list[-match(selected.lineages[1], extant.active.list)] #Remove first 
        extant.active.list = extant.active.list[-match(selected.lineages[2], extant.active.list)] # Remove second
        extant.active.list[keff-1] = node.number # Add new nodes 
        
        keff = keff - 1 # Coalesence removes lineages
        k = k-1
        time  = time-nxt.coal
      } else if(min(waiting.times) == nxt.latent) {
        ######### Latent Event 
        latent.node = sample(extant.active.list[1:keff], size = 1)
        possible.events = select_possible_events(nxt.event.type = "A")
        graph.edge = determine_edge(latent.node, host.arg, nxt.latent, edge.type = "A", time, sample.time, tips, possible.events, tip.df)
        
        ## Creating new entry into the ARG
        if(graph.edge$from.node.type == "T") { 
          recombo = NA 
        } else { 
          recombo = find_recombo_site(node = latent.node, host.arg = host.arg$current.arg, list.entry = host.arg$num.list, host.arg.list = host.arg$arg.list)
        }
        new.entry = create_arg_entry(from.node = latent.node, to.node = latent.node, edge.length = graph.edge$edge.length,
                                     nxt.event = nxt.latent, from.node.type = graph.edge$from.node.type, edge.type = graph.edge$edge.type, recombo = recombo, time = time)
        host.arg = add_arg_entry(host.arg, new.entry, N) 
        
        # Update time 
        track.entry=data.table(lineages.total = k, lineages.effective = keff, event = "AL", time.lapse=nxt.latent, sim.time = time, pop.size = pop.size) %>%
          mutate_if(is.factor, as.character) %>% as.data.table
        track.list = add_track_entry(track.list, track.entry, N)
        
        # Latent Event/Bookkeeping 
        extant.latent.list[(k-keff + 1)] = latent.node
        extant.active.list = extant.active.list[-match(latent.node, extant.active.list)]
        keff = keff -  1 # Recombination adds a lineage to active lineages
        k = k # The total number doesn't change  
        time = time-nxt.latent  # Move backward in time 
      } else {
        # Activate latent node event 
        # Edge Length  - A latent already has an internal node)
        active.node = sample(extant.latent.list[(k-keff)], size = 1)
        possible.events = select_possible_events(nxt.event.type = "L")
        graph.edge = determine_edge(active.node, host.arg, nxt.activate, edge.type = "L", time, sample.time, tips, possible.events, tip.df = tip.df)
        
        ## Creating new entry into the ARG
        if(graph.edge$from.node.type == "T") { 
          recombo = NA 
        } else { 
          recombo = find_recombo_site(node = active.node, host.arg = host.arg$current.arg, list.entry = host.arg$num.list, host.arg.list = host.arg$arg.list)
        } 
        new.entry = create_arg_entry(from.node = active.node, to.node = active.node, edge.length = graph.edge$edge.length,
                                     nxt.event = nxt.activate, from.node.type = graph.edge$from.node.type,  edge.type = graph.edge$edge.type, recombo = recombo, time = time)
        host.arg = add_arg_entry(host.arg, new.entry, N) 
        
        ## Add track entry
        track.entry=data.table(lineages.total = k, lineages.effective = keff, event = "LA", time.lapse=nxt.coal, sim.time = time, pop.size = pop.size) %>%
          mutate_if(is.factor, as.character) %>% as.data.table
        track.list = add_track_entry(track.list, track.entry, N) 
        
        ## Activation event/Bookkeeping  
        extant.active.list[keff + 1] = active.node
        extant.latent.list = extant.latent.list[-match(active.node, extant.latent.list)]
        keff = keff +  1 # Activation adds a lineage to active lineages
        k = k # The total number doesn't change  
        time = time-nxt.activate 
      }
    } else if(nxt.switch == 21) {
      #browser()
      time = 21
      extant.active.list = c(extant.active.list, extant.latent.list)
      extant.active.list = extant.active.list[which(extant.active.list != "")]
      keff = length(extant.active.list)
      k = keff
      
      #### List to enter in the entries into the arg 
      extant.latent.list = extant.latent.list[which(extant.latent.list != "")]
      if(length(extant.latent.list) != 0) { 
        for(i in 1:length(extant.latent.list)) { 
          active.node = extant.latent.list[i]
          nxt.activate = 0
          possible.events = select_possible_events(nxt.event.type = "L")
          graph.edge = determine_edge(active.node, host.arg, nxt.activate, edge.type = "L", time, sample.time= final.time, tips,possible.events, tip.df = tip.df)
          
          ## Creating new entry into the ARG
          if(graph.edge$from.node.type == "T") { 
            recombo = NA 
          } else { 
            recombo = find_recombo_site(node = active.node, host.arg = host.arg$current.arg, list.entry = host.arg$num.list, host.arg.list = host.arg$arg.list)
          }
          new.entry = create_arg_entry(from.node = active.node, to.node = active.node, edge.length = graph.edge$edge.length,from.node.type = graph.edge$from.node.type,
                                       nxt.event = 0, edge.type = graph.edge$edge.type, recombo = recombo, time = time)
          host.arg = add_arg_entry(host.arg, new.entry, N) 
          ## Add track entry
          track.entry=data.table(lineages.total = k, lineages.effective = keff, event = "LA", time.lapse=nxt.coal, sim.time = time, pop.size = pop.size) %>%
            mutate_if(is.factor, as.character) %>% as.data.table
          track.list = add_track_entry(track.list, track.entry, N) 
        }
      }
      extant.latent.list = NULL
      latent.phase = 1 
      if (nxt.boundary == 21) nxt.boundary = 20.99999
    } else if (nxt.switch == 0 ) {
      time=-1
    } else if (nxt.switch == nxt.boundary) {
     # browser() 
      if(sweep == FALSE) {
        sweep = TRUE
        pop.realm = pop.realm + 1
        time = population.boundaries$start[pop.realm]
        if(is.na(time)) time = 0
        nxt.boundary = population.boundaries$end[pop.realm]
        if(is.na(nxt.boundary)) nxt.boundary = 0
        pop.size = population.boundaries$pop.size[pop.realm]
        #print(paste0("Switching Pop Realm -- Pop Size = ", pop.realm, " Time - ", time, " Nxt Boundary:", nxt.boundary))
      } else {
        sweep = FALSE
        time = population.boundaries$end[pop.realm]
        if(is.na(time)) time = 0
        nxt.boundary = population.boundaries$start[pop.realm + 1]
        if(is.na(nxt.boundary)) nxt.boundary = 0
        #print(paste0("Switching Pop Realm -- Pop Size = ", pop.realm, " Time - ", time, " Nxt Boundary:", nxt.boundary))
      }
    } else {
      time = -1
    }  # Break out of loop 
  }  
  
  track.entry=data.table(lineages.total = k, lineages.effective = keff, event = as.factor("F"), time.lapse=NA, sim.time = time, pop.size = pop.size) %>%
    mutate_if(is.factor, as.character) %>% as.data.table
  track.list = add_track_entry(track.list, track.entry, N)
  
  #print(paste0("Ending with k - ", k))
  #print(paste0("Ending with keff - ", keff))
  #browser()
  tip.df[1:length(tips), 1] = tips
  tip.df[1:length(tips), 2] =  final.time
  return(list(arg = host.arg, extant.active.lineages = extant.active.list, 
              extant.latent.lineages = extant.latent.list, track = track.list,
              node.number = node.number, tips = tips, pop.realm = pop.realm, 
              sweep.status = sweep, latent.phase=2, nxt.boundary = nxt.boundary, tip.df = tip.df))#, time.track = time.track))
}


sim_arg_addsample_full <- function(sample.k, sample.time, final.time, nxt.sample.time, extant.active,
                                   extant.latent, host.arg, track.list, prev.tips, starting.node,
                                   alpha, re.rate, latent.size, activate.rate, population.boundaries, 
                                   pop.realm, sweep.status, latent.phase, nxt.boundary, tip.df) {

  # Simulates coalescent/recombination/latency within previous sampled populations
  # Time starts at next time sample, adding the extant existing nodes and new ones 
  # Active lineages start with new samples plus those remaining from previous time point 
  # ##browser()
  N = 10000
  extant.active.list = vector(mode = "character", length = 500) # Create active list size
  extant.active.list[1:sample.k] = as.character(seq(from = starting.node, to = (starting.node+sample.k-1), by = 1))
  prev.active.length = length(which(extant.active !=""))
  if(prev.active.length!=0) extant.active.list[sample.k + 1:length(extant.active[which(extant.active != "")])] = extant.active[which(extant.active != "")]
  
  extant.latent.list = vector(mode = "character", length = 500) # Create latent list size
  if(!is.null(extant.latent)) extant.latent.list[1:length(extant.latent)] = extant.latent
  
  tips.now = as.character(as.numeric(seq(from = starting.node, to = starting.node + sample.k-1, by = 1)))
  tips = c(prev.tips, tips.now)
  
  node.number = starting.node+sample.k # Keeps track of how many nodes are in the ARG
  k = length(extant.active.list[which(extant.active.list != "")]) + length(extant.latent.list[which(extant.latent.list != "")]) # k is the number of extant lineages
  keff = length(extant.active.list[which(extant.active.list != "")]) # Effective K (only the active)
  
  if(latent.phase == 2) time = sample.time  else time = 21
  pop.realm = pop.realm
  sweep  = sweep.status
  pop.size = population.boundaries$pop.size[pop.realm]
  nxt.boundary = nxt.boundary
  
  #print(paste0("Time;", time, "---Next Sample Time: ", nxt.sample.time))
  while (time > nxt.sample.time) {
    
    ## Coalescent dependent on population size
    if(sweep == FALSE) {
      if(keff > 1) nxt.coal = inv_cumf_coal(u=runif(1), time = time, alpha=alpha, k=keff) else nxt.coal = 10000
      pop.size = alpha + 1*time
    } else {
      if(keff > 1 ) nxt.coal = exp_coalesent_pop(k=keff, pop.size = population.boundaries$pop.size[pop.realm]) else nxt.coal = 10000
    }
    
    if(keff > 0) nxt.recombo = exp_recombo(k=keff, re.rate = re.rate) else nxt.recombo = 10e12
    if(latent.phase == 2 & keff > 0) nxt.latent = inv_cumf_latent(runif(1), time = time, keff = keff, rho = activate.rate, alpha = alpha) else nxt.latent = 10e12
    if(k-keff > 0) nxt.activate = exp_activate2(klat = (k-keff), rho = activate.rate, N_l = latent.size) else nxt.activate = 10e12
    
    waiting.times = c(nxt.coal, nxt.recombo, nxt.latent, nxt.activate)
    waiting.times = waiting.times[which(!is.na(waiting.times))]
    
    ## Determine if next boundary switch is because of population sweep or sampling
    if(latent.phase == 2) nxt.switch = max(nxt.boundary, nxt.sample.time, 21) else nxt.switch = max(nxt.boundary, nxt.sample.time)
    if(time - min(waiting.times) > nxt.switch) { # Will the next event take place before the time boundary for population size ?
      if(min(waiting.times) == nxt.recombo) {
        ###### Recombination Event
        selected.lineage = sample(extant.active.list[1:keff], size = 1) # select random lineage
        node.number = node.number + 1  # designate new node number
        possible.events = select_possible_events(nxt.event.type = "R")
        graph.edge = determine_edge(selected.lineage, host.arg, nxt.recombo, edge.type = "R", time, sample.time, tips, possible.events, tip.df=tip.df)
        ## Creating new entry into the ARG
        recombo.site = sample(1:700, size = 1)
        new.entry1 = create_arg_entry(from.node = selected.lineage, to.node =  paste0(format(node.number, scientific = FALSE), ".1"), edge.length = graph.edge$edge.length,
                                      nxt.event = nxt.recombo, from.node.type = graph.edge$from.node.type, edge.type = graph.edge$edge.type, recombo = recombo.site, time = time)
        new.entry2 = create_arg_entry(from.node = selected.lineage, to.node = paste0(format(node.number,scientific = FALSE), ".2"), edge.length = graph.edge$edge.length,
                                      nxt.event = nxt.recombo, from.node.type = graph.edge$from.node.type, edge.type = graph.edge$edge.type, recombo = recombo.site, time = time)
        
        host.arg = add_arg_entry(host.arg, new.entry1, N)
        host.arg = add_arg_entry(host.arg, new.entry2, N)
        
        ## Add track entry
        track.entry=data.table(lineages.total = k, lineages.effective = keff, event = as.factor("R"), time.lapse=nxt.recombo, sim.time = time, pop.size = pop.size) %>%
          mutate_if(is.factor, as.character) %>% as.data.table
        track.list = add_track_entry(track.list, track.entry, N)
        
        # Recombination Event/Bookeeping
        extant.active.list = extant.active.list[-match(selected.lineage, extant.active.list)] # Remove the recombinant
        extant.active.list[c(keff,keff+1)] = c(paste0(format(node.number, scientific = FALSE),".1"), paste0(format(node.number,scientific= FALSE),".2"))
        keff = keff + 1 # Recombination adds a lineage to active lineages
        k = k + 1 #  Recombination adds a lineage to total
        time = time-nxt.recombo  # Move backward in time
      } else if(min(waiting.times) == nxt.coal) {
        ###### Coalescent Event
        selected.lineages = sample(extant.active.list[1:keff], size = 2, replace = FALSE) # Two lineages coalesce
        node.number = node.number + 1 # designate new node number
        
        #First Entry -- check if coalescing lineage is a tip
        possible.events = select_possible_events(nxt.event.type = "C")
        graph.edge.1 = determine_edge(selected.lineage=selected.lineages[1], host.arg, nxt.coal,
                                      edge.type = "C", time, sample.time, tips, possible.events, tip.df = tip.df)
        new.entry.1 = create_arg_entry(from.node = selected.lineages[1], to.node = as.factor(node.number), edge.length = graph.edge.1$edge.length,
                                       nxt.event = nxt.coal, from.node.type = graph.edge.1$from.node.type, edge.type = graph.edge.1$edge.type, recombo = NA, time = time)
        
        #Second Entry -- check if a coalescing lineage is a tip
        graph.edge.2 = determine_edge(selected.lineage=selected.lineages[2], host.arg, nxt.coal,
                                      edge.type = "C", time, sample.time, tips, possible.events, tip.df = tip.df)
        new.entry.2 = create_arg_entry(from.node = selected.lineages[2], to.node = as.factor(node.number), edge.length = graph.edge.2$edge.length,
                                       nxt.event = nxt.coal, from.node.type = graph.edge.2$from.node.type, edge.type = graph.edge.2$edge.type, recombo = NA, time = time)
        
        host.arg = add_arg_entry(host.arg, new.entry.1, N)
        host.arg = add_arg_entry(host.arg, new.entry.2, N)
        
        ## Add track entry
        track.entry=data.table(lineages.total = k, lineages.effective = keff, event = "C", time.lapse=nxt.coal, sim.time = time,pop.size = pop.size) %>%
          mutate_if(is.factor, as.character) %>% as.data.table
        track.list = add_track_entry(track.list, track.entry, N)
        
        # Coalescent Event/Bookkeeping
        extant.active.list = extant.active.list[-match(selected.lineages[1], extant.active.list)] #Remove first
        extant.active.list = extant.active.list[-match(selected.lineages[2], extant.active.list)] # Remove second
        extant.active.list[keff-1] = node.number # Add new nodes
        
        keff = keff - 1 # Coalesence removes lineages
        k = k-1
        time  = time-nxt.coal
      } else if(min(waiting.times) == nxt.latent) {
        
        ######### Latent Event
        latent.node = sample(extant.active.list[1:keff], size = 1)
        possible.events = select_possible_events(nxt.event.type = "A")
        graph.edge = determine_edge(latent.node, host.arg, nxt.latent, edge.type = "A", time, sample.time, tips, possible.events, tip.df = tip.df)
        
        ## Creating new entry into the ARG
        if(graph.edge$from.node.type == "T") {
          recombo = NA
        } else {
          recombo = find_recombo_site(node = latent.node, host.arg = host.arg$current.arg, list.entry = host.arg$num.list, host.arg.list = host.arg$arg.list)
        }
        new.entry = create_arg_entry(from.node = latent.node, to.node = latent.node, edge.length = graph.edge$edge.length,
                                     nxt.event = nxt.latent, from.node.type = graph.edge$from.node.type, edge.type = graph.edge$edge.type, recombo = recombo, time = time)
        host.arg = add_arg_entry(host.arg, new.entry, N)
        
        # Update time
        track.entry=data.table(lineages.total = k, lineages.effective = keff, event = "AL", time.lapse=nxt.latent, sim.time = time, pop.size = pop.size) %>%
          mutate_if(is.factor, as.character) %>% as.data.table
        track.list = add_track_entry(track.list, track.entry, N)
        
        # Latent Event/Bookkeeping
        extant.latent.list[(k-keff + 1)] = latent.node
        extant.active.list = extant.active.list[-match(latent.node, extant.active.list)]
        keff = keff -  1 # Recombination adds a lineage to active lineages
        k = k # The total number doesn't change
        time = time-nxt.latent  # Move backward in time
      } else {
        # Activate latent node event
        # Edge Length  - A latent already has an internal node)
        active.node = sample(extant.latent.list[(k-keff)], size = 1)
        possible.events = select_possible_events(nxt.event.type = "L")
        graph.edge = determine_edge(active.node, host.arg, nxt.activate, edge.type = "L", time, sample.time, tips, possible.events, tip.df = tip.df)
        
        ## Creating new entry into the ARG
        if(graph.edge$from.node.type == "T") {
          recombo = NA
        } else {
          recombo = find_recombo_site(node = active.node, host.arg = host.arg$current.arg, list.entry = host.arg$num.list, host.arg.list = host.arg$arg.list)
        }
        new.entry = create_arg_entry(from.node = active.node, to.node = active.node, edge.length = graph.edge$edge.length,
                                     nxt.event = nxt.activate, from.node.type = graph.edge$from.node.type,  edge.type = graph.edge$edge.type, recombo = recombo, time = time)
        host.arg = add_arg_entry(host.arg, new.entry, N)
        
        ## Add track entry
        track.entry=data.table(lineages.total = k, lineages.effective = keff, event = "LA", time.lapse=nxt.coal, sim.time = time, pop.size = pop.size) %>%
          mutate_if(is.factor, as.character) %>% as.data.table
        track.list = add_track_entry(track.list, track.entry, N)
        
        ## Activation event/Bookkeeping
        extant.active.list[keff + 1] = active.node
        extant.latent.list = extant.latent.list[-match(active.node, extant.latent.list)]
        keff = keff +  1 # Activation adds a lineage to active lineages
        k = k # The total number doesn't change
        time = time-nxt.activate
      }
    } else if(nxt.switch == 21) { # If next event is where latent cells are activated (going back in time)
      
      time = 21
      extant.active.list = c(extant.active.list, extant.latent.list)
      extant.active.list = extant.active.list[which(extant.active.list != "")]
      keff = length(extant.active.list)
      k = keff
      
      #### List to enter in the entries into the arg
      extant.latent.list = extant.latent.list[which(extant.latent.list != "")]
      if(length(extant.latent.list) != 0) {
        for(i in 1:length(extant.latent.list)) {
          active.node = extant.latent.list[i]
          nxt.activate = 0
          possible.events = select_possible_events(nxt.event.type = "L")
          graph.edge = determine_edge(active.node, host.arg, nxt.activate, edge.type = "L", time, sample.time= final.time, tips,possible.events, tip.df=tip.df)
          
          ## Creating new entry into the ARG
          if(graph.edge$from.node.type == "T") {
            recombo = NA
          } else {
            recombo = find_recombo_site(node = active.node, host.arg = host.arg$current.arg, list.entry = host.arg$num.list, host.arg.list = host.arg$arg.list)
          }
          new.entry = create_arg_entry(from.node = active.node, to.node = active.node, edge.length = graph.edge$edge.length,from.node.type = graph.edge$from.node.type,
                                       nxt.event = 0, edge.type = graph.edge$edge.type, recombo = recombo, time = time)
          host.arg = add_arg_entry(host.arg, new.entry, N)
          ## Add track entry
          track.entry=data.table(lineages.total = k, lineages.effective = keff, event = "LA", time.lapse=nxt.coal, sim.time = time, pop.size = pop.size) %>%
            mutate_if(is.factor, as.character) %>% as.data.table
          track.list = add_track_entry(track.list, track.entry, N)
        }
      }
      extant.latent.list = NULL
      latent.phase = 1
      if (nxt.boundary == 21) nxt.boundary = 20.99999
      #nxt.boundary = 0 ## trying this to see if it works .... 
    } else if (nxt.switch == 0 ) {
    
      time=-1
    } else if (nxt.switch == nxt.boundary) { 
    
      if(sweep == FALSE) {
        sweep = TRUE
        pop.realm = pop.realm + 1
        time = population.boundaries$start[pop.realm]
        nxt.boundary = population.boundaries$end[pop.realm]
        if(is.na(nxt.boundary)) nxt.boundary = 0
        pop.size = population.boundaries$pop.size[pop.realm]
        #print(paste0("Switching Pop Realm -- Pop Size = ", pop.realm, " Time - ", time, " Nxt Boundary:", nxt.boundary))
      } else {
        sweep = FALSE
        pop.realm = pop.realm
        time = population.boundaries$end[pop.realm]
        nxt.boundary = population.boundaries$start[pop.realm + 1]
        if(is.na(nxt.boundary)) nxt.boundary = 0
        #print(paste0("Switching Pop Realm -- ", pop.realm, " Time - ", time, " Nxt Boundary:", nxt.boundary))
      }
    } else {
      time = -1
    }  # Break out of loop
  }
  
  track.entry=data.table(lineages.total = k, lineages.effective = keff, event = as.factor("F"), time.lapse=NA, sim.time = time, pop.size = pop.size) %>%
    mutate_if(is.factor, as.character) %>% as.data.table
  track.list = add_track_entry(track.list, track.entry, N)
 
  tip.df[(length(tips)-length(tips.now)+1):length(tips),1]= tips.now
  tip.df[(length(tips)-length(tips.now)+1):length(tips), 2] = sample.time
  #print(paste0("Ending with k - ", k))
  #print(paste0("Ending with keff - ", keff))
  #print(paste0("Ending time", time))
  return(list(arg = host.arg, extant.active.lineages = extant.active.list, 
              extant.latent.lineages = extant.latent.list, track = track.list,
              node.number = node.number, tips = tips,pop.realm = pop.realm, sweep.status = sweep, 
              latent.phase = latent.phase,nxt.boundary =nxt.boundary, tip.df=tip.df))
}


sim_arg_full <- function(simulation_params, population.boundaries) {
  with(simulation_params, { 
   
    ## Simulation to generate an ARG of HIV within one host with longitudinal samples 
    ## Time starts when last sequences are sampled 
    node.list = vector(mode = "numeric", length = length(sample.times) + 1)
    
    
    tip.df = data.frame(tip.number = character(sum(sample.sizes)+starting.k),
                        sample.time = character(sum(sample.sizes)+starting.k), stringsAsFactors = FALSE)
    
    time.arg = sim_arg_initial_full(starting.k = starting.k, starting.node = 1, final.time = final.time, 
                                    sample.time = final.time, nxt.sample.time = sample.times[1],
                                    alpha = alpha, re.rate = re.rate, latent.size = latent.size, 
                                    activate.rate = activate.rate, population.boundaries, tip.df = tip.df) 
    
    extant.active = time.arg$extant.active.lineages
    extant.latent = time.arg$extant.latent.lineages
  
    
    node.list[1] = time.arg$node.number
    if(sample.times[1] !=0) { 
      if(length(sample.times)>=2) {
        for(i in 1:(length(sample.times)-1)) {
          if(!exists("time.arg")) break() 
          time.arg=sim_arg_addsample_full(sample.k = sample.sizes[i],
                                          sample.time = sample.times[i],
                                          final.time = final.time,
                                          nxt.sample.time = sample.times[i+1],
                                          extant.active = extant.active,
                                          extant.latent = extant.latent,
                                          host.arg = time.arg$arg,
                                          track.list = time.arg$track,
                                          starting.node = time.arg$node.number + 1,
                                          prev.tips = time.arg$tips,
                                          alpha = alpha,
                                          re.rate = re.rate,
                                          latent.size = latent.size,
                                          activate.rate = activate.rate,
                                          population.boundaries = population.boundaries,
                                          pop.realm = time.arg$pop.realm,
                                          sweep.status = time.arg$sweep.status,
                                          latent.phase = time.arg$latent.phase,
                                          nxt.boundary = time.arg$nxt.boundary,
                                          tip.df = time.arg$tip.df)
          node.list[i+1] = time.arg$node.number
          extant.active = time.arg$extant.active.lineages
         
          extant.latent = time.arg$extant.latent.lineages
        }
      }
   
      time.arg = sim_arg_addsample_full(sample.k = sample.sizes[length(sample.times)],
                                        sample.time = sample.times[length(sample.times)],
                                        final.time = final.time,
                                        nxt.sample.time = 0,
                                        extant.active = extant.active,
                                        extant.latent = extant.latent,
                                        host.arg = time.arg$arg,
                                        track.list = time.arg$track,
                                        starting.node = time.arg$node.number + 1,
                                        prev.tips = time.arg$tips,
                                        alpha = alpha,
                                        re.rate = re.rate,
                                        latent.size = latent.size,
                                        activate.rate = activate.rate,
                                        population.boundaries = population.boundaries,
                                        pop.realm = time.arg$pop.realm,
                                        sweep.status = time.arg$sweep.status,
                                        latent.phase = time.arg$latent.phase,
                                        nxt.boundary = time.arg$nxt.boundary,
                                        tip.df = time.arg$tip.df)
    }
    
    ### Return full ARG, track lineages,
    host.arg = time.arg$arg
    track.lineages = time.arg$track
    extant.latent = time.arg$extant.latent.lineages
    time.entry = time.arg$time.track
    
    if(length(extant.latent[which(extant.latent!="")])==0) valid = TRUE else valid = FALSE
    #print(valid)
    if (length(time.arg$extant.active.lineages) ==  1) {
      #browser()
      return(list(arg = host.arg, track.lineages = track.lineages, nodes = node.list,valid = valid, time = time.entry)) 
    } else {
      #browser()
      return(NULL)
    }
    
  })
}
default_params <- function(beta = 1,
                           alpha  = 0,
                           starting.k = 10, # Last sample of lineage size
                           final.time = 180*6, # Latest sample time 
                           sample.times = rev(seq(from = 180, to = 180*5, by = 180)), # This will be a vector of times,
                           sample.sizes = c(5, 5, 5, 5,5), # This is the vector of samples corresponding to the same times 
                           re.rate = .01,
                           activate.rate = .001,
                           latent.size = 1000,
                           num.collapses = 1,
                           pop.sizes = 5,
                           collapse.length = 1e-12)
return(as.list(environment()))



## Example for running ARG 
# final.time = 6*180
# selection.frequency = 200
# selection.strength = 50
# 
# num.collapses = length(seq(from = 0, to = final.time, by = selection.frequency))
# activate.rate = 1.5
# recombo.rate =.01
# 
# population.boundaries = generate_pophistory_linear_current(num.collapses,
#                                                            collapse.length = 5,
#                                                            pop.size = selection.strength,
#                                                            final.time = final.time,
#                                                            selection.frequency = selection.frequency )
# 
# sim_params <- default_params(final.time = final.time,
#                              collapse.length = 5,
#                              num.collapses = num.collapses,
#                              sample.sizes = rep(x = 20, final.time-180),
#                              sample.times = rev(seq(from = 180, to = (final.time-180), by = 180)),
#                              re.rate = recombo.rate,
#                              latent.size = 1000,
#                              activate.rate = activate.rate,
#                              pop.size = selection.strength,
#                              starting.k = 20)
# 
# 
# 
# trial <- sim_arg_full(sim_params, population.boundaries)
