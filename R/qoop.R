qoop <- function(fp, MODEL=1, VERSION=1, MIX=FALSE) {
  pop     <- popEPP$new(fp, MODEL, VERSION, MIX)
  for (i in 2:fp$SIM_YEARS) {
    pop$year <- i
    pop$aging()
    pop$add_entrants()
    pop$deaths()
    # pop$migration()
    if (exists("risk", where=fp))
      pop$disease_model()
    pop$adjust_pop()
    # pop$save_outputs() # save prevalence
  }
  class(pop) <- "qoop"
  return(pop)
}