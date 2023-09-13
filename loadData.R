library(RSiena)

adj.matrix.1 = as.matrix(read.table("data/kapfti1.dat"))
adj.matrix.2 = as.matrix(read.table("data/kapfti2.dat"))

job.status <- read.table("data/kapfa_stat.dat")[[1]]
high.status <- read.table("data/kapfat.dat")[[2]]

create.kapferer.data <- function() {
  dummy_actor  <- diag(39)
  for(i in 1:9)   assign(paste("dummy_actor_0", i, sep = ""), coCovar(dummy_actor[,i], centered=F))
  for(i in 10:39) assign(paste("dummy_actor_",  i, sep = ""), coCovar(dummy_actor[,i], centered=F))
  
  array.data = array(NA_integer_, dim=c(39,39,2))
  array.data[,,1] = adj.matrix.1; array.data[,,2] = adj.matrix.2
  network = sienaDependent(array.data)
  status = coCovar(high.status)
  
  sienaDataCreate( network, status,
                   dummy_actor_01, dummy_actor_02, dummy_actor_03, dummy_actor_04, dummy_actor_05, dummy_actor_06, dummy_actor_07, dummy_actor_08, dummy_actor_09, dummy_actor_10,
                   dummy_actor_11, dummy_actor_12, dummy_actor_13, dummy_actor_14, dummy_actor_15, dummy_actor_16, dummy_actor_17, dummy_actor_18, dummy_actor_19, dummy_actor_20,
                   dummy_actor_21, dummy_actor_22, dummy_actor_23, dummy_actor_24, dummy_actor_25, dummy_actor_26, dummy_actor_27, dummy_actor_28, dummy_actor_29, dummy_actor_30,
                   dummy_actor_31, dummy_actor_32, dummy_actor_33, dummy_actor_34, dummy_actor_35, dummy_actor_36, dummy_actor_37, dummy_actor_38, dummy_actor_39)
}

data <- create.kapferer.data()
rm(create.kapferer.data)


