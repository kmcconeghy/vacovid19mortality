
# Setup -------------------------------------------------------------------
  # important packages for the project
  library(data.table) # helps with large data manipulations
  library(splines) # for spline terms
  library(parglm) # glm model alternative, allows parallel computing 

  ## load data ----
    # d_ipw = [INSERT YOUR DATA HERE]
    # longitudinal dataset, person-period with baseline time fixed variables and time-varying
    # Outcome is censoring event (vaccination)
    # includes both baseline and time-fixed covariates

  ## sample ----
    # FOR TEST RUNS
    if (F) {
      d_ipw = d_ipw %>%
        inner_join(
        slice_sample(distinct(., PatientICN), prop=0.1),
        .,
        by = 'PatientICN')
      
      d_surv = inner_join(d_surv, distinct(d_ipw, PatientICN))
    }

# IPW weights ----
  ## Set up IPW data
    setDT(d_ipw, key=c('PatientICN'))

    # create some model vars, interactions etc.
    d_ipw = d_ipw[, `:=`(
      intercept = 1,
      day2 = day*day,
      dayXhr = day*dx_hr_1p, 
      day2Xhr = day*day*dx_hr_1p,
      can_t2 = can_t*can_t,
      can_b2 = can_b*can_b,
      dem_age2 = dem_age*dem_age
    )]

    # generate matrix with regression variables
    d_xmat = select(d_ipw, 
                    intercept, day, day2, dayXhr, day2Xhr, dx_hr_1p,
                    can_t, can_t2, can_b, can_b2,
                    hsr_ip, hsr_hsp_t, hsr_opt_t, 
                    hsr_rx, hsr_rx_t, dx_Mace, dx_Mace_t,
                    dem_age, dem_age2, vacc_flu, 
                    hsr_lb_any, hsr_lb_trpn, dx_UnstablyHoused,
                    dem_afam, dx_WeightLoss, dx_Mace, 
                    dx_DMany, dx_Renal, dx_NeuroOther, 
                    dem_afam) %>%
      data.matrix(.)

    # once I make matrix I don't need those variables in my dataset
    d_ipw = d_ipw[, .(PatientICN, index, day, vacc)]
    
  # Estimate model ----
    d_glm_ipw = parglm.fit(
      y = d_ipw$vacc, # vacc is vaccine event (censoring)
      x = d_xmat, # matrix
      family=binomial(),
      model=F,
      # control sets some regression step parameters
      control = parglm.control(method="FAST",
                               nthreads=10,
                               epsilon = 1e-7,
                               maxit=25)
    ) 

  ## compute weights ----
    # Pr (vacc=1) at each timepoint from the model
    d_ipw$pr_vacc = d_glm_ipw$fitted.values
    
    #d_surv = [INSERT YOUR DATA HERE]
    # Outcome now death event
    # longitudinal dataset, person-period
    # Only needs baseline characteristics

      # groups for merging datasets
        grp_nms = c('PatientICN', 'treat', 'index')
        setDT(d_surv, key=grp_nms)
    
      # Data.table left join for pr_vacc for speed
        d_surv[, 
               pr_vacc := d_ipw[d_surv, 
                                on = .(PatientICN, index, day), 
                                x.pr_vacc]
               ]
    
      # IF TREAT = 1; Pr(no censor) = pr_den; 
      # IF TREAT = 0; Pr(no cens) = 1 - pr_den;
        d_surv[, pr_den_a := treat*(pr_vacc) + (1-treat)*(1-pr_vacc)]
    
      # If assigned no vaccine, 1 / pr(no cens) for teach timepoint ----
      # if assign vaccine, 1 / 1 before end of grace (cannot censor)
      # you get vaccine at end of grace (day 7), 1 / pr(no cens)
      # you get vaccine before day 7, 1 / 1  thereafter (because cannot censor after you get vaccine)
        grace = 7
      
        d_surv[, `:=`(ipw_us = fcase(treat == 0, 1 / pr_den_a,
                                          treat == 1 & day < grace, 1 / 1,
                                          treat == 1 & day == grace & t_vacc==grace, 1 / pr_den_a,
                                          treat == 1 & day == grace & t_vacc<grace, 1 / 1,
                                          treat == 1 & day > grace, 1 / 1))] 
      
      # with ipw calculated, now take cumulative product
        d_surv[, `:=`(ipw_us_c = cumprod(ipw_us)), by=grp_nms] 
        quantile(d_surv$ipw_us_c, c(seq(0.97, 1, 0.0025)))

      # truncate IPW
        trunc = 0.001
        
        d_tau = quantile(d_surv$ipw_us_c, 1-trunc)
        d_surv$ipw_us_t = d_surv$ipw_us_c
        d_surv$ipw_us_t[d_surv$ipw_us_c >= d_tau] = d_tau
      
        quantile(d_surv$ipw_us_t, c(seq(0.97, 1, 0.0025)))
    
# Survival curves ----

    # generate spline for time ----
      d_bs_day = as.data.frame(bs(1:60, df=5)) %>%
        mutate(day=1:60)
      
      colnames(d_bs_day) = c('day1', 'day2', 'day3', 'day4', 'day5', 'day')
      setDT(d_bs_day)
      
      d_surv = d_bs_day[d_surv, on = c('day')]
      
      d_surv = d_surv[, `:=`(
        intercept = 1,
        can_b2 = can_b*can_b,
        dem_age2 = dem_age*dem_age
      )]
    
      d_xmat = select(d_surv, 
                      intercept, day1, day2, day3, day4, day5,
                      can_b, can_b2, dem_age, vacc_flu, hsr_ip, 
                      hsr_rx, hsr_lb_any, hsr_lb_trpn, dx_UnstablyHoused,
                      dx_WeightLoss, dx_Mace, dx_DMany, dx_Renal, 
                      dx_NeuroOther, dx_Anemia, dx_PVD, dx_hr_1p, 
                      dx_Pulmonary, dem_afam) %>%
        data.matrix(.)
      
      d_surv = setDF(d_surv) %>%
        select(event, day, treat, ipw_us_t)
  
  ## Estimate survival separately for each model ----
    d_glm_pe_t0 = parglm.fit(
      y = d_surv$event[d_surv$treat==0],
      x = d_xmat[d_surv$treat==0, ],
      family=binomial(),
      weights = d_surv$ipw_us_t[d_surv$treat==0],
      model=F,
      control = parglm.control(method="FAST",
                               nthreads=10,
                               epsilon = 1e-7,
                               maxit=25)
    ) 
    
    d_glm_pe_t1 = parglm.fit(
      y = d_surv$event[d_surv$treat==1],
      x = d_xmat[d_surv$treat==1, ],
      family=binomial(),
      weights = d_surv$ipw_us_t[d_surv$treat==1],
      model=F,
      control = parglm.control(method="FAST",
                               nthreads=10,
                               epsilon = 1e-7,
                               maxit=25)
    ) 
    
  ## Survival probabilities ----
    d_surv$pr_ev = NA_real_
    d_surv$pr_ev[d_surv$treat==1] = d_glm_pe_t1$fitted.values
    d_surv$pr_ev[d_surv$treat==0] = d_glm_pe_t0$fitted.values
    
    d_res = summ_survprob(d_surv)

# End ----
