setMembers <- function(a, where="public", x, y) {
    if (length(x) > 1)
        mapply(function(x, y) a$set(where, x, y, overwrite=T), x, y)
    else 
        a$set(where, x, y, overwrite=T)
    invisible()
}

# EPP parameter super class: later improve the generating of fp as well
#' @importFrom methods new
eppFP <- R6::R6Class("eppfp", class=F, cloneable=F, portable=F, lock_objects=F,
    public = list(
        p = NULL,
        initialize = function(fp) {
            p <<- fp[-which(names(fp)=="ss")]
            list2env(fp$ss, self)
        }
    )
)

# POP CLASS
# -----------------------------------------------------------------------------
# preferablly all fields should have a leading character
popEPP <- R6::R6Class("popepp", class=F, cloneable=F, portable=F, inherit=eppFP,
    lock_objects=F,
    public = list(
        data              = "array",
        VERSION           = "character",
        MODEL             = "numeric",
        MIX               = "logical",
        year              = "numeric",
        birth_age         = "numeric", birth_agrp        = "numeric",
        prev15to49        = "vector",  incid15to49       = "vector",
        prev              = 0,         incid             = 0,
        entrantprev       = "vector",  pregprevlag       = "vector",
        birthslag         = "array",   infections        = "array",
        hivdeaths         = "array",   natdeaths         = "array",
        popadjust         = "array",   hivp_entrants_out = "array",
        incrate15to49_ts  = "vector",  prev15to49_ts     = "vector", # can be array
        data_db           = "array",    # virgin population
        data_active       = "array",
        artcov            = numeric(2), # initially no one on treatment
        prev_last         = 0, # last time step prevalence
        hiv_sx_prob       = "array",
        hiv_mr_prob       = "array",
        adj_prob          = "array",
        rvec              = "vector")
)

setMembers(popEPP, "public", "initialize",
function(fp, MODEL=1, VERSION="R", MIX=F) {
    super$initialize(fp)
    # 'Initialize pop array'
    MODEL       <<- MODEL
    data        <<- array(0, c(pAG, NG, pDS, PROJ_YEARS))
    data[,,1,1] <<- p$basepop
    
    # Outputs
    entrantprev   <<- numeric(PROJ_YEARS)
    prev15to49    <<- incid15to49  <<- pregprevlag <<- entrantprev
    popadjust     <<- array(0, c(pAG, NG, PROJ_YEARS))
    infections    <<- hivdeaths <<- natdeaths <<- popadjust
    prev15to49_ts <<- incrate15to49_ts <<- rep(NA, length(p$rvec))
    hivp_entrants_out <<- array(0, c(NG, PROJ_YEARS))
    
    # use in model
    birthslag     <<- p$birthslag
})

# named list functions as methods
popFunc <- c(
set_data = function(FUN="+", x, AG=T, NG=T, DS=T, YEAR=T) {
    FUN <- match.fun(FUN)
    data[AG, NG, DS, YEAR] <<- FUN(data[AG, NG, DS, YEAR], x)
}, 

put = function(x, AG=T, NG=T, DS=T, YEAR=T) {
    'This can be done directly with `<-` using pop_object$data'
    data[AG, NG, DS, YEAR] <<- x
}, 

get = function(YEAR=NULL, AG=T, NG=T, DS=T) {
    'This can be done directly with `[` using pop_object$data'
    if (is.null(YEAR))
        data[AG, NG, DS, TRUE]
    else 
        data[AG, NG, DS, YEAR]
},

update_infection = function(infect) {
    infect <- infect * DT
    data[,, hivn.idx, year] <<- data[,, hivn.idx, year] - infect
    data[,, hivp.idx, year] <<- data[,, hivp.idx, year] + infect
    infections[,,year] <<- infections[,,year] + infect 
    incid15to49[year]  <<- incid15to49[year] + sum(infect[p.age15to49.idx,])  
},

remove_hiv_death = function(cd4_mx, hivpop, artpop) {
    # death by age group
    artD <- artpop$f_death
    hivD <- hivpop$f_death
    if (MODEL==2) { # add deaths from inactive population
      artD <- artD + artpop$f_death_db
      hivD <- hivD + hivpop$f_death_db
    }
    dbyAG <- DT * (colSums(hivD) + colSums(artD,,2))
    # deaths by single-year
    hiv_mx <- dbyAG / sumByAGs(data[,, hivp.idx, year], ag.idx)
    hiv_mx[is.nan(hiv_mx)] <- 0
    dbyA_pr <- apply(hiv_mx, 2, rep, h.ag.span)
    hivdeaths[,,year] <<- hivdeaths[,,year] + data[,, hivp.idx, year] * dbyA_pr
    data[,, hivp.idx, year] <<- data[,, hivp.idx, year] * (1 - dbyA_pr)
    if (MODEL == 2) {
      hivdeaths[1:pDB,,year] <<- hivdeaths[1:pDB,,year] +
        data_db[,, hivp.idx, year] * dbyA_pr[1:pDB, ]
      data_db[,,hivp.idx,year] <<- data_db[,,hivp.idx,year] * (1 - dbyA_pr[1:pDB,])
    }
},

adjust_pop = function() {
  popadjust[,,year] <<- p$targetpop[,,year] / rowSums(data[,,,year],,2)
  data[,,,year] <<- sweep(data[,,,year], 1:2, popadjust[,,year], "*")
},

aging = function() {
    data[-c(1,pAG),,,year] <<- data[-(pAG-1:0),,,year-1]
    data[pAG,,,year] <<- data[pAG,,,year-1] + data[pAG-1,,,year-1] # open ended
},

add_entrants = function() {
    ## Add lagged births into youngest age group and update population
    if (exists("popadjust", where=p) & p$popadjust) {
        if (MODEL==0)
            healthy <- p$entrantpop[,year-1]
        else {
            healthy <- p$entrantpop[,year-1]*(1-p$entrantprev[,year])
            positiv <- p$entrantpop[,year-1]*   p$entrantprev[,year]
        }
    } 
    else {
        if (MODEL==0)
            healthy <- birthslag[,year-1] * p$cumsurv[,year-1] / 
                       p$paedsurv_lag[year-1] + p$cumnetmigr[,year-1]
        else {
            survivedBirth <- birthslag[,year-1]*p$cumsurv[,year-1]
            prev_now <- p$entrantprev[,year]  # avoid repeat accesses
            imm_now  <- p$cumnetmigr[,year-1] 
            healthy <- survivedBirth * (1-prev_now/p$paedsurv_lag[year-1]) + 
                        imm_now*(1 - pregprevlag[year-1]*p$netmig_hivprob)
            positiv <- (survivedBirth + imm_now) * prev_now
        }
    }
    # save and update pop        
    data[1,, hivn.idx, year] <<- healthy
    if (MODEL != 0) {
        data[1,, hivn.idx, year] <<- healthy
        data[1,, hivp.idx, year] <<- positiv
    }
    if (MODEL==2) {  # this is unwise
        data_db[1,,hivn.idx,year] <<- healthy
        data_db[1,,hivp.idx,year] <<- positiv
    }
    if (MODEL!=0) {
        entrantprev[year]        <<- sum(positiv)/sum(healthy+positiv)
        hivp_entrants_out[,year] <<- positiv
    }
},

entrant_art = function() { # return these for updating HIV and ART pop
    yesno <- c(p$entrantartcov[,year], 1 - p$entrantartcov[,year])
    out <- data[1,, hivp.idx, year] * yesno
    out # 1:2 ART+, 3:4 ART-
},

deaths = function() {
    death_now <- sweep(data[,,,year], 1:2, 1 - p$Sx[,,year], "*")
    data[,,,year] <<- data[,,,year] - death_now
    natdeaths[,,year] <<- rowSums(death_now,,2)
},

migration = function() {
    netmigsurv <- p$netmigr[,,year] * (1 + p$Sx[,,year]) / 2
    mr.prob <- 1 + netmigsurv / rowSums(data[,,,year],,2)
    data[,,,year] <<- sweep(data[,,,year], 1:2, mr.prob, "*")
},

update_fertile = function() { # only on active pop
    update_active_pop_to(year)
    two_years  <- data_active + get_active_pop_in(year-1)
    birth_age  <<- rowSums(two_years[p.fert.idx, f.idx,])/2 * p$asfr[,year]
    birth_age  <<- rowSums(two_years[p.fert.idx, f.idx,])/2 * p$asfr[,year]
    birth_agrp <<- sumByAG(birth_age, ag.idx, TRUE, p.fert.idx)
    births     <- p$srb[,year] * sum(birth_agrp)
    if (year + AGE_START <= PROJ_YEARS)
      birthslag[, year + AGE_START - 1] <<- births
},

cal_prev_pregant = function(hivpop, artpop) { # only on active pop
    years   <- year - 1:0
    update_active_pop_to(year)
    two_years  <- data_active + get_active_pop_in(year-1)
    meanWomen <- two_years[p.fert.idx, f.idx, hivn.idx] / 2
    hivn <- sumByAG(meanWomen, ag.idx, TRUE, p.fert.idx)
    hivp <- rowMeans(hivpop$get(AG=h.fert.idx, NG=f.idx, YEAR=years),,2)
    art  <- rowMeans(artpop$get(AG=h.fert.idx, NG=f.idx, YEAR=years),,3)
    pregprev <- sum(birth_agrp * 
      (1 - hivn / (hivn + colSums(p$frr_cd4[,,year] * hivp) + 
      colSums(p$frr_art[,,,year] * art,,2)))) / sum(birth_age)
    pregprevlag[year + AGE_START - 1] <<- pregprev
},

save_prev_n_inc = function() {
    prev15to49[year]  <<- sum(data[p.age15to49.idx,, hivp.idx, year]) / 
                          sum(data[p.age15to49.idx,,         , year])
    update_active_pop_to(year-1)
    incid15to49[year] <<- incid15to49[year] /
                          sum(data_active[p.age15to49.idx,,hivn.idx])
    prev[year]  <<- sum(data[,,hivp.idx,year]) / sum(data[,,,year])
    incid[year] <<- incid15to49[year] / sum(data_active[,,hivn.idx]) # toBfixed
}, 

infect_mix = function(hivpop, artpop, ii) {
    ts <- (year-2)/DT + ii
    update_active_pop_to(year)
    hiv_treated       <- sweep(data_active[,,hivp.idx], 2, artcov, '*')
    hiv_not_treated   <- data_active[,,hivp.idx] - hiv_treated
    transm_prev <- (hiv_not_treated + hiv_treated * (1 - p$relinfectART)) / 
                    rowSums(data_active,,2) # prevalence adjusted for art
    # +intervention effects and time epidemic start
    w  <- p$iota * (p$proj.steps[ts] == p$tsEpidemicStart)
    transm_prev <- rvec[ts] * transm_prev + w
    # sweep over sexual mixing matrices
    ir_m <- rowSums(sweep(p$mat_m, 2, transm_prev[, f.idx], "*")) # male
    ir_f <- rowSums(sweep(p$mat_f, 2, transm_prev[, m.idx], "*")) # female
    irmf   <- cbind(ir_m, ir_f)
    # if (exists("f_fun", fp)) # that fun
    #   ir <- ir * fp$f_fun
    infections.ts <- irmf * data_active[,,hivn.idx]

    incrate15to49_ts[,,ts] <<- transm_prev
    prev15to49_ts[ts] <<- prevlast <<- sum(data[,,hivp.idx,year])/sum(data[,,,year])
    infections.ts
},

disease_model = function() {
    # To be conclude, just for playing
    new_obese <- sweep(data[,,1,year], 1:2, p$risk, "*")
    data[,,1,year] <<- data[,,1,year] - new_obese
    data[,,2,year] <<- data[,,2,year] + new_obese
}

)

setMembers(popEPP, "public", names(popFunc), popFunc)