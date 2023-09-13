algorithm.estimation <- sienaAlgorithmCreate(projname = 'estimation', cond=F) 

model.stand <- siena07(algorithm.estimation, data=data, effects=effects.fixef.stand)
model.stand

# full transitivity
model.fullt <- siena07(algorithm.estimation, data=data, effects=effects.fixef.fullt)
model.fullt

# no transitivity:
model.notrt <- siena07(algorithm.estimation, data=data, effects=effects.fixef.notrt)
model.notrt

# no status
model.nosts <- siena07(algorithm.estimation, data=data, effects=effects.fixef.nosts)
model.nosts





