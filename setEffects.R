effects.fixef.stand <- getEffects(data)

effects.fixef.stand <- includeEffects( effects.fixef.stand,  transTrip)
effects.fixef.stand <- includeEffects( effects.fixef.stand,  egoX, altX, simX,  interaction1 = "status" )

effects.ranef.stand <- effects.fixef.stand

for(i in 1:9)   effects.ranef.stand <- includeEffects( effects.ranef.stand,  egoX,  interaction1 = paste0("dummy_actor_0",i) )
for(i in 10:39) effects.ranef.stand <- includeEffects( effects.ranef.stand,  egoX,  interaction1 = paste0("dummy_actor_", i) )

# alternative for complete model of transitivity alternative to balance:
effects.fixef.fullt <- includeEffects( effects.fixef.stand,  outAct,  include = T ); effects.fixef.fullt
effects.ranef.fullt <- includeEffects( effects.ranef.stand,  outAct,  include = T )

# no transitivity:
effects.fixef.notrt <- includeEffects( effects.fixef.stand,  transTrip,  include = F); effects.fixef.notrt
effects.ranef.notrt <- includeEffects( effects.ranef.stand,  transTrip,  include = F)

# no status:
effects.fixef.nosts <- includeEffects( effects.fixef.stand,  egoX, altX, simX,  interaction1 = "status" , include=F); effects.fixef.nosts
effects.ranef.nosts <- includeEffects( effects.ranef.stand,  egoX, altX, simX,  interaction1 = "status" , include=F)






