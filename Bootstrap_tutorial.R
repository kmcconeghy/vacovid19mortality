
# Setup -------------------------------------------------------------------
  library(data.table)
  library(splines)
  library(parglm)
  options(dplyr.summarise.inform=F)
  
  #d_start = [INSERT YOUR starting dataset]

# Bootstrap function ----
  boot_it = function(d_ipw) {

      # Key by Person ID 
      setDT(d_ipw, key = 'ID')

      # A list of unique persons
      d_ids = unique(d_ipw[,.(ID)])

      # generate a random value for each person (~Poisson(1))
      d_ids[, freqwt := rpois(.N, 1)]

      # add back to dataset
      d_ipw[, freqwt := d_ids[d_ipw, 
                              on = .(ID),
                              x.freqwt]
            ]

    # make some covariates for regression
      d_ipw = d_ipw[, `:=`(
        intercept = 1,
        day2 = day*day,
        dayXhr = day*dx_hr_1p,
        day2Xhr = day*day*dx_hr_1p,
        can_t2 = can_t*can_t,
        can_b2 = can_b*can_b,
        dem_age2 = dem_age*dem_age
      )]

    # Matrix
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

      # dont need covariates (this is a big dataset so I drop, but step not necessary)
      #d_ipw = d_ipw[, .(ID, index, day, vacc, freqwt)]

      # fit censoring model
      d_glm_ipw = parglm.fit(
        y = d_ipw$vacc,
        x = d_xmat,
        family=binomial(),
        weights = d_ipw$freqwt,
        start = readRDS(here('output', 
                             'full_adjus_wt', 
                             'c19d_ipw_coefs.Rds')),
        model=F,
        control = parglm.control(method="FAST",
                                 nthreads=10,
                                 epsilon = 1e-7,
                                 maxit=25)
      ) 
      
    ## compute weights 
      d_ipw$pr_vacc = d_glm_ipw$fitted.values
      
      d_surv = read_fst(here('dbdf', 'c19d_chrt_lng_surv.fst'))
      
      setDT(d_surv)
      
      d_surv[, 
             pr_vacc := d_ipw[d_surv, 
                              on = .(PatientICN, index, day), 
                              x.pr_vacc]
      ]
      
      d_surv[, 
             freqwt := d_ids[d_surv, 
                              on = .(PatientICN), 
                              x.freqwt]
      ]
      
      rm(d_ipw, d_glm_ipw, d_ids)
      
      grp_nms = c('PatientICN', 'treat', 'index')
      
      setDT(d_surv, key=grp_nms)
    
      # IF TREAT = 1; Pr(no censor) = pr_den; 
      # IF TREAT = 0; Pr(no cens) = 1 - pr_den;
      
      d_surv[, pr_den_a := treat*(pr_vacc) + (1-treat)*(1-pr_vacc)]
    
      # If assigned no vaccine, 1 / pr(no cens) for teach timepoint 
      # if assign vaccine, 1 / 1 before end of grace (cannot censor)
      # you get vaccine at end of grace (day 7), 1 / pr(no cens)
      # you get vaccine before day 7, 1 / 1  thereafter (because cannot censor after you get vaccine)
      grace = 7
      
      d_surv[, `:=`(ipw_us = fcase(treat == 0, 1 / pr_den_a,
                                   treat == 1 & day < grace, 1 / 1,
                                   treat == 1 & day == grace & t_vacc==grace, 1 / pr_den_a,
                                   treat == 1 & day == grace & t_vacc<grace, 1 / 1,
                                   treat == 1 & day > grace, 1 / 1))] 
      
      
      d_surv[, `:=`(ipw_us_c = cumprod(ipw_us)), by=grp_nms] 
      
      trunc = 0.001
      d_tau = quantile(d_surv$ipw_us_c, 1-trunc)
      d_surv$ipw_us_t = d_surv$ipw_us_c
      d_surv$ipw_us_t[d_surv$ipw_us_c >= d_tau] = d_tau
      
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
      
      d_surv = select(d_surv, PatientICN, index, event, day, treat, freqwt, ipw_us_t)
    
      d_glm_pe_t0 = parglm.fit(
        y = d_surv$event[d_surv$treat==0],
        x = d_xmat[d_surv$treat==0, ],
        family=binomial(),
        start = readRDS(here('output','full_adjus_wt', 'c19d_surv_t0_coefs.Rds')),
        weights = d_surv$ipw_us_t[d_surv$treat==0]*d_surv$freqwt[d_surv$treat==0],
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
        start = readRDS(here('output', 'full_adjus_wt', 'c19d_surv_t1_coefs.Rds')),
        weights = d_surv$ipw_us_t[d_surv$treat==1]*d_surv$freqwt[d_surv$treat==1],
        model=F,
        control = parglm.control(method="FAST",
                                 nthreads=10,
                                 epsilon = 1e-7,
                                 maxit=25)
      ) 
      
      d_surv$pr_ev = NA_real_
      d_surv$pr_ev[d_surv$treat==1] = d_glm_pe_t1$fitted.values
      d_surv$pr_ev[d_surv$treat==0] = d_glm_pe_t0$fitted.values
      
      d_surv[, `:=`(pr_surv = cumprod(1 - pr_ev)), by=grp_nms] 
      
      d_res = d_surv %>%
        group_by(treat, day) %>%
          summarize(pr_ev = weighted.mean(1-pr_surv, freqwt)) %>%
        ungroup %>%
        pivot_wider(., id_cols =c('day'), 
                    names_from = treat, 
                    names_prefix = 'pr_ev_',
                    values_from = pr_ev
        ) %>%
        mutate(cid = pr_ev_1 - pr_ev_0,
               cir = pr_ev_1 / pr_ev_0)
      
      return(d_res)
  }
  
# Run boot ----
  #d_tst = boot_it(d_start)
  
  run_boots = function(i) {
    d_bres = boot_it(d_start)
    cat('.')
    
    saveRDS(d_bres, 
            here('output', 'full_adjus_wt', 'boot', 
                 '20230808', paste0('c19d_b.', i, '.Rds')))
  }
  
   set.seed(as.integer(ymd('2023-08-08')))
  # 
   map(1:250, ~run_boots(.))
  
  # set.seed(as.integer(ymd('2023-03-02')))
  # 
  # map(32:200, ~run_boots(.))
  
  # set.seed(as.integer(ymd('2023-03-05')))
  # map(201:250, ~run_boots(.))
  
