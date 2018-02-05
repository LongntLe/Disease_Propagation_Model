##################################################################################
##################################################################################
#Credit: Sherry Towers, ASU, Polymatheia
#A script emplying MC parameter sweeping algorithm
#Modified for visual proof of non-convexity of least square function
#and to fit for Congo Measles Data
#By Long Le
##################################################################################

##################################################################################
# par(pch=20) sets a solid round dot style for the plots
# The chron package contains utilities to calculate dates
##################################################################################
rm(list = ls(all = TRUE))  # resets R to fresh
require("chron")
require("sfsmisc")
par(pch=20)  
source("sir_func.R")
set.seed(732552)

x <- c(25,37,92,209,167,149,68,19)
t <- 30*c(1,2,3,4,5,6,7,8,9)
adat <- as.data.frame(cbind(t,x))
#plot(df, pch=21, bg="red", xlab="time(month)", ylab="number of reports")
colnames(adat) <- c("times_of_observed", "incidence_observed")

##################################################################################
# read in 2007 influenza surveillance data (confirmed cases) for
# the Midwest (CDC region 5) from
# http://www.cdc.gov/flu/weekly/regions2007-2008/datafinalHHS/whoregX.htm
# where X=5 is midwest
# X=1 is northeast
# X=2 is NY and NJ
# X=3 are eastern seabord states like PA, DE, etc
#
# the weeks are number of weeks relative to Jan 1, 2007
# week 1 in 2007 ended Jan 6, 2007
##################################################################################
#adat = read.table("midwest_influenza_2007_to_2008.dat",header=T,sep=",")
#cat("\n")
#cat("The data file contains: ",names(adat),"\n")
#cat("\n")

##################################################################################
# The CDC data is weekly, with the date of each point corresponding to the
# date of the end of the week over which the data were collected.
# Let's convert these dates to time in days, relative to Jan 1, 2007
# We will be using this vector of dates, vtime_data, to obtain the model estimates
# of the incidence at that time.
# adat$week is relative to Jan 1, 2007, with week #1 occuring the first week in
# January, 2007.
##################################################################################
#adat$time_in_days_rel_jan_1_2007 = julian(1,6,2007)+(adat$week-1)*7-julian(1,1,2007)

##################################################################################
# Specifically for this data:
# make sure we are far enough into the season that there is at least one case per week
##################################################################################
#adat=subset(adat,week>=47)  
incidence_observed = adat$incidence_observed #adat$B
times_of_observed = adat$times_of_observed #adat$time_in_days_rel_jan_1_2007
time_binning = min(diff(times_of_observed))

##################################################################################
##################################################################################
# now, set up the iterations of the Monte Carlo method
# At each iteration, we will randomly sample a hypothesis for the
# reproduction number, R0, and the time-of-introduction of the virus to the
# population, t0.
#
# With these hypotheses, we will solve for the model predicted incidence
# and calculate the least squares statistic comparing this to the observed
# incidence.  We store the hypotheses for R0 and t0 in the vectors vR0 and vt0,
# and the resulting least squares statistic in the vector vleastsq
#
# At each iteration, we'll check if the predicted incidence is the best fit
# so far, and if so, we'll store that in vbest_leastsq_fit_incidence_prediction
#
# best_leastsq_so_far keeps track of the best-fit least squares so far obtained
# in the iterations.
##################################################################################
vR0 = numeric(0)
vt0 = numeric(0)
vleastsq = numeric(0) 

best_leastsq_so_far = 1e10 
vbest_leastsq_fit_incidence_prediction = rep(0,length(incidence_observed))

niter = 10000  
for (iter in 1:niter){
  
  ###############################################################################
  # This process is computationally intensive, so once in a while during the
  # iterations it is nice to inform user the script is doing something, 
  # and not hung
  ###############################################################################
  if (iter%%100==0){
    cat("Doing iteration ",iter," out of ",niter,"\n")  
  }
  
  ###############################################################################
  ################################################################################
  # set up the model parameters
  # npop is approximately the population of IL IN MI MN OH WI (CDC region 5)
  # I_0 is the initial number infected
  # R_0 is the inital number recovered and immune
  # S_0 is the susceptibles
  #
  # 1/gamma is the average recovery period of the disease
  # R0      is the reproduction number
  # t0      is the time-of-introduction of the disease to the population, measured
  #         in days from Jan 1, 2007
  #
  # For the SIR model, R0=beta/gamma, thus given our hypotheses for gamma and R0,
  # we calculate beta as beta=R0*gamma (note that beta and gamma appear in the model
  # equations, not gamma and R0, which is why we need to calculate beta given R0
  # and gamma).
  ###############################################################################
  npop = 78000000 
  I_0 = 1      
  R_0 = 0      
  S_0 = npop-I_0-R_0
  
  gamma = 1/15
  
  ###############################################################################
  # randomly sample R0 and t0 uniformly
  ###############################################################################
  R0 = runif(1,3,6)             
  #t0 = as.integer(runif(1,160,(min(times_of_observed)-time_binning)))
  
  ###############################################################################
  # or, you can use the Normal distribution to preferentially sample close to
  # a particular value
  ###############################################################################
  #R0 = rnorm(1,1.20,0.10) 
  #t0 = as.integer(rnorm(1,200,20))
  
  ###############################################################################
  # calculate beta for the SIR model, and fill the vparameters and inits
  # vectors that get passed to the lsoda method that solves our system of
  # equations
  ###############################################################################
  beta  = R0*gamma
  
  vparameters = c(gamma=gamma,beta=beta)
  inits = c(S=S_0,I=I_0,R=R_0)
  
  ###############################################################################
  # We get the model solution for all days from t0 to the last week of the
  # data time series.  If t0 is greater than the minimum date in the data time
  # series, we need to print out a warning, because the time of introduction had
  # to be before we actually started observing cases in the data!
  ###############################################################################
  #t0 = adat$time[1]
  t0 = 0
  tmin = t0
  tmax = max(times_of_observed)
  
  if (tmin>(min(times_of_observed)-time_binning)){
    cat("\n")
    cat("**************************************************************************\n")
    cat("**************************************************************************\n")
    cat("The time-of-introduction is _after_ the first cases appeared!",t0,min(times_of_observed)-time_binning,"\n")
    cat("**************************************************************************\n")
    cat("**************************************************************************\n")
    cat("\n")
  }
  tmin = min(t0,min(times_of_observed)-time_binning)
  
  vt = seq(tmin,tmax)
  
  ###############################################################################
  # Now solve the system of differential equations numerically with lsoda in the 
  # deSolve package.  Put the results in solved_model
  # The derivative_calc_func for the SIR model is in the sir_func.R script
  ###############################################################################
  solved_model = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters)) 
  
  ###############################################################################
  # the B influenza data is incidence, not prevalence (actually, the B influenza
  # data is the true incidence times the fraction that the CDC actually confirms)
  #
  # The incidence over some time step is the difference in S over that time
  # step (in a model with no births or deaths, immigration or emigration, or 
  # recurring susceptibility).
  #
  # The time step are derived from the dates at the end of the week of each
  # data point (vtime_data)
  #
  # solved_model$time%in%vtime_data returns the indices of the elements of simodel$time that
  # are also found in vtime_data
  ###############################################################################
  tmin_data = min(times_of_observed)-time_binning
  tmax_data = max(times_of_observed)
  vtime_data = seq(tmin_data,tmax_data,time_binning)
  
  susceptible_predicted = solved_model$S[solved_model$time%in%vtime_data]  
  incidence_predicted = -diff(susceptible_predicted)
  
  ###############################################################################
  # from the model estimate of the incidence and the data, we can 
  # estimate the fraction of cases that were confirmed
  ###############################################################################
  frac_confirmed = sum(incidence_observed)/sum(incidence_predicted)
  
  ###############################################################################
  # normalize the model prediction so area under curve
  # equals the sum of the data incidence
  ###############################################################################
  incidence_predicted = incidence_predicted*frac_confirmed 
  
  ###############################################################################
  # now calculate the least-squares
  # statistic that compares the data to this model calculated
  # under a particular hypothesis of R0 and t0
  ###############################################################################
  
  if (length(incidence_predicted)==length(incidence_observed)
      &!is.na(sum(incidence_predicted))){
    
    ########################################################################### 
    # calculate the least squares statistic
    ########################################################################### 
    leastsquares = sum((incidence_predicted-incidence_observed)^2)
    vR0 = append(vR0,R0)
    vt0 = append(vt0,t0)
    vleastsq = append(vleastsq,leastsquares)
    
    if (leastsquares<best_leastsq_so_far){
      best_leastsq_so_far = leastsquares
      vbest_leastsq_fit_incidence_prediction = incidence_predicted
      R0_best = R0
      t0_best = t0
      cat("The best value of R0 so far is:",R0_best,"\n")
      cat("The best value of t0 so far is:",t0_best,"\n")
    }
    
    ######################################################################## 
    # plot the best-fit results every once in a while
    # cex is the point size
    ######################################################################## 
    if (iter%%100==0){
      text_main = paste("Confirmed Measles Cases, Congo, 2017 season:\n result of",iter,"Monte Carlo fit iterations")
      #text_main = paste("Least Square Calculation with varied reproduction number")
      mult.fig(4,main=text_main,oma=c(1,2,4,1))
      
      num_points_to_show = 250
      ymax = max(vleastsq)
      #if (length(vleastsq)>num_points_to_show) ymax = sort(vleastsq)[num_points_to_show]
      l = which(vleastsq<=ymax)
      
      lmin = which.min(vleastsq)
      
      plot(vR0[l]
           ,vleastsq[l]
           ,ylab="Least squares"
           ,xlab="R0 hypothesis"
           ,main=paste("Best-fit R0 so far:",round(R0_best,3))
           ,col.main=3)
      points(vR0[lmin],vleastsq[lmin],col=3,cex=2)
      
      plot(vt0[l]
           ,vleastsq[l]
           ,ylab="Least squares"
           ,xlab="t0 hypothesis (days rel Jan 1, 2017)"
           ,main=paste("Best-fit t0 so far:",t0_best)
           ,col.main=3)
      points(vt0[lmin],vleastsq[lmin],col=3,cex=2)

      ymax = max(c(incidence_observed,vbest_leastsq_fit_incidence_prediction))
      plot(times_of_observed
           ,incidence_observed
           ,ylim=c(0,1.2*ymax)
           ,xlab="Time, in weeks relative to Jan 1, 2017"
           ,ylab="Incidence"
           ,cex=2)
      lines(times_of_observed
            ,vbest_leastsq_fit_incidence_prediction
            ,col=2
            ,lwd=5)
    }
    
  } # end check that the predicted incidence vector is the same length as the observed and doesn't contain NA's
} # end loop over the Monte Carlo iterations

