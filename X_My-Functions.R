#################################
# Script for modifications to 
# Alex's L3 Assess functions
# Mainly to make the plots look
# prettier for publicatoin
################################

#### Generic functions from Alex ####
Get_xaxis_scale <- function(x_data) {
  
  # modified from code provided by
  # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/
  
  xmax_data = max(x_data)
  xmin_data = min(x_data)
  
  xpow = log10(xmax_data-xmin_data)
  xint = 10 ^ (xpow-round(xpow,0))
  
  if (xint>=0 & xint<2.5) {
    xint = 0.2
  }
  if (xint>=2.5 & xint<5) {
    xint = 0.5
  }
  if (xint>=5 & xint<7.5) {
    xint = 1
  }
  if (xint>=7.5) {
    xint = 2
  }
  
  xint = xint * 10^round(xpow,0) # major ticks
  xmin = xint * round(xmin_data/xint,0)
  xmax = xint * (round(xmax_data / xint,0) + 1)
  
  results = list(xmin = xmin,
                 xmax = xmax,
                 xint = xint)
  
}

Get_yaxis_scale <- function(y_data) {
  
  # modified from code provided by
  # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/
  
  ymax_data = max(1.1 * y_data)
  ymin_data = min(y_data)
  
  ypow = log10(ymax_data-ymin_data)
  yint = 10 ^ (ypow-round(ypow,0))
  
  if (yint>=0 & yint<2.5) {
    yint = 0.2
  }
  if (yint>=2.5 & yint<5) {
    yint = 0.5
  }
  if (yint>=5 & yint<7.5) {
    yint = 1
  }
  if (yint>=7.5) {
    yint = 2
  }
  
  yint = yint * 10^round(ypow,0) # major ticks
  ymin = yint * round(ymin_data/yint,0)
  ymax = yint * (round(ymax_data / yint,0) + 1)
  
  results = list(ymin = ymin,
                 ymax = ymax,
                 yint = yint)
  
}

AddAxesAndTickLabelsToPlot <- function(xmin=NA, xmax, xint, ymin=NA, ymax, yint, cexval=NA, cexaxisval=NA, lwdval=NA, lineval=NA, lasval=NA) {
  
  if (is.na(xmin)) xmin=0
  if (is.na(ymin)) ymin=0
  if (is.na(cexval)) cexval=1
  if (is.na(cexaxisval)) cexaxisval=1
  if (is.na(lwdval)) lwdval=1
  if (is.na(lasval)) lasval=1
  if (is.na(lineval)) lineval=0
  
  axis(1, at = seq(xmin, xmax, xint), line = lineval, labels = F)
  axis(2, at = seq(ymin, ymax, yint), line = lineval, labels = F)
  axis(1, at = seq(xmin, xmax, xint), lwd=lwdval, labels=T, line=lineval, cex=cexval, cex.axis=cexaxisval, las=lasval)
  axis(2, at = seq(ymin, ymax, yint), lwd=lwdval, labels=T, line=lineval, cex=cexval, cex.axis=cexaxisval, las=lasval)
  
}


#### Modified functions ####

X_PlotLengthBasedCatchCurve_RetCatch <- function(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                               SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                               MaxAge, NatMort, TimeStep, MainLabel, xaxis_lab, yaxis_lab, xmax, xint,
                                               ymax, yint, PlotCLs, FittedRes, nReps, Error.Colour) {
  
  .pardefault <- par(no.readonly = TRUE) # store current par settings
  
  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, SelectivityType, ObsRetCatchFreqAtLen,
                                        lbnd, ubnd, midpt, SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
  }
  
  params = res$params
  vcov.params = res$vcov.Params
  ExpCatchAtLen = res$ModelDiag$ExpRetCatchPropInLenClass
  ObsRelCatchAtLen = ObsRetCatchFreqAtLen/sum(ObsRetCatchFreqAtLen)
  
  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = nReps, params, vcov.params))
  EstPropAtLen.sim = as.matrix(data.frame(matrix(nrow = nReps, ncol = length(midpt))))
  
  for (j in 1:nReps) {
    params = unlist(sims[j,])
    CatchCurveType=1 #1=length-based, 2=age and length based
    Res=AgeAndLengthBasedCatchCurvesCalcs(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityVec, DiscMort, MaxAge, NatMort, TimeStep)
    
    EstPropAtLen.sim[j,] = unlist(Res$ExpRetCatchPropInLenClass)
    cat("j",j,'\n')
  }
  
  EstProp.sim = apply(EstPropAtLen.sim, 2, median)
  EstProp.sim_low = apply(EstPropAtLen.sim, 2, quantile, probs = 0.025)
  EstProp.sim_up = apply(EstPropAtLen.sim, 2, quantile, probs = 0.975)
  
  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "Prop. (retained catch)"
  xlims = Get_xaxis_scale(ubnd)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  ylims = Get_yaxis_scale(ObsRelCatchAtLen)
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint
  
  plot(midpt, ObsRelCatchAtLen, "p", main=MainLabel, cex.main=1.2, pch=16, cex=0.8, xaxt = "n", yaxt = "n",
       xlab=list(xaxis_lab,cex=1.2),ylab=list(yaxis_lab,cex=1.2), frame=F, xlim=c(0,xmax), ylim=c(0,ymax))
  if (is.data.frame(GrowthParams)) {
    lines(midpt, Res$ExpRetCatchPropInLenClass_Fem, lty="dotted", col="dark grey")
    lines(midpt, Res$ExpRetCatchPropInLenClass_Mal, lty="dotted", col="blue")
    legend("bottomright", legend=c("Female","Male"), y.intersp = 1.0, inset=c(0.05,0.05),
           cex = 0.8, bty="n", lty="dotted", col=c("dark grey","blue"))
  }
  if (PlotCLs == TRUE) {
    sm1 = spline(Res$midpt, EstProp.sim_low, n=100, method="natural")
    sm2 = spline(Res$midpt, EstProp.sim_up, n=100, method="natural")
    sm1$y[which(sm1$y<0)]=0
    sm2$y[which(sm2$y<0)]=0
    if (!is.na(MLL)) {
      sm1$y[which(sm1$x<MLL)]=0
      sm2$y[which(sm2$x<MLL)]=0
    }
    x = c(sm1$x, rev(sm2$x)) # using shading for 95% CLs
    y = c(sm1$y, rev(sm2$y))
    polygon(x,y, col=Error.Colour,border=NA)
    sm1_x <- sm1$x
    sm1_y <- sm1$y
    sm2_x <- sm2$x
    sm2_y <- sm2$y
  }
  points(midpt, ObsRelCatchAtLen, col="grey20", pch=16, cex=0.8)
  points(midpt, ExpCatchAtLen, col="grey20", pch=1, cex=0.8)
  axis(1, at = seq(0, xmax, xint), line = 0.2, labels = F)
  axis(2, at = seq(0, ymax, yint), line = 0.2, labels = F)
  axis(1, at = seq(0, xmax, xint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  axis(2, at = seq(0, ymax, yint), lwd = 0, labels = T, line = 0, cex.axis = 1, las = 1)
  
  # inverse logit transformed value
  params = res$params # from point estimate, not from resampled values
  Fval = round(1/(1+exp(-params[1])),2)
  uprbnd <- res$ParamEst[1,3]
  lwrbnd <- res$ParamEst[1,2]
  Fest = bquote("F =" ~ .(Fval) ~ y^-1 ~"("*.(lwrbnd)~"-"~.(uprbnd) ~ y^-1*")")
  if (SelectivityType==1) {
    legend("topright", pch=-1, legend=as.expression(Fest),
           lty="solid",col="black", bty='n', cex=0.8,lwd=-1, y.intersp=1.2, adj=0)
  }
  if (SelectivityType==2) {
    L50est=paste("L50 =",round(exp(params[2]),0),"mm")
    L95est=paste("L95 - L50 =",round(exp(params[3]),0),"mm")
    legend("topright", pch=-1, legend=c(as.expression(Fest), L50est, L95est),
           lty="solid",col="black", bty='n', cex=0.8,lwd=-1, y.intersp=1.2)
  }
  legend("topleft", legend=c("Observed","Estimated"), y.intersp = 1.0, inset=c(0.13,0),
         lty=1, cex = 0.8, bty="n", seg.len = 0, pch=c(16,1), col=c("grey20","grey20"))
  
  par(.pardefault)
  
  
  to.plot <- cbind(ExpCatchAtLen, ObsRelCatchAtLen) %>% 
    as.data.frame() %>% 
    pivot_longer(everything(.), names_to="Obs.Est", values_to="Prop") %>% 
    group_by(Obs.Est) %>% 
    mutate(midpoint = midpt) %>% 
    ungroup()
  
  params.for.plot <- list()
  
  params.for.plot[["ExpCatchAtLen"]] <- ExpCatchAtLen
  params.for.plot[["ObsRelCatchAtLen"]] <- ObsRelCatchAtLen
  params.for.plot[["plotting.points"]] <- to.plot
  params.for.plot[["x"]] <- x
  params.for.plot[["y"]] <- y
  params.for.plot[["Fest"]] <- Fest
  params.for.plot[["sm1_x"]] <- sm1_x
  params.for.plot[["sm1_y"]] <- sm1_y
  params.for.plot[["sm2_x"]] <- sm2_x
  params.for.plot[["sm2_y"]] <- sm2_y
  
  
  return(params.for.plot)
}

#### Plotting mortality at length ####
Get_xaxis_scale <- function(x_data) {
  
  # modified from code provided by
  # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/
  
  xmax_data = max(x_data)
  xmin_data = min(x_data)
  
  xpow = log10(xmax_data-xmin_data)
  xint = 10 ^ (xpow-round(xpow,0))
  
  if (xint>=0 & xint<2.5) xint = 0.2
  if (xint>=2.5 & xint<5) xint = 0.5
  if (xint>=5 & xint<7.5) xint = 1
  if (xint>=7.5) xint = 2
  
  xint = xint * 10^round(xpow,0) # major ticks
  xmin = xint * round(xmin_data/xint,0)
  xmax = xint * (round(xmax_data / xint,0) + 1)
  
  results = list(xmin = xmin,
                 xmax = xmax,
                 xint = xint)
  
}


X_PlotLengthBasedCatchCurve_Mortality <- function(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityAtLen,
                                                PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                                MaxAge, NatMort, TimeStep, xmax, xint, ymax, yint, FittedRes_1, FittedRes_2, colour_1, colour_2) {
  
  #.pardefault <- par(no.readonly = TRUE) # store current par settings
  
  # if model already fitted, can input results rather than refit

    res_1 = FittedRes_bruv
    res_2 = FittedRes_dpird
  
  
  CatchCurveType=1 #1=length-based, 2=age and length based
  
  params_1 = res_1$params
  res_1 = AgeAndLengthBasedCatchCurvesCalcs(params_1, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                          MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep)
  params_2 = res_2$params
  res_2 = AgeAndLengthBasedCatchCurvesCalcs(params_2, DistnType, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge, CatchCurveType,
                                            MLL, SelectivityType, lbnd, ubnd, midpt, SelectivityAtLen, DiscMort, MaxAge, NatMort, TimeStep)
  
  GrowthModelType = res_1$GrowthModelType
  
  xlims = Get_xaxis_scale(midpt)
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  
  yaxis_lab = expression(paste(italic("F") ~ (year^{-1})))
  xaxis_lab = "Length (mm)"
  
  # F at age
  # if (GrowthModelType ==  1 | GrowthModelType == 3) { # single sex
  #   ylims = Get_yaxis_scale(res$FAtLen)
  #   if (is.na(ymax)) ymax = ylims$ymax
  #   if (is.na(yint)) yint = ylims$yint
  #   plot(midpt, res$FAtLen, "l", main="F (retention + discard mortality)", cex.main=1.0, pch=1, cex=0.6,
  #        xaxt = "n", yaxt = "n", xlab=list(xaxis_lab,cex=1.2),ylab=NA, frame=F, xlim=c(0,xmax),
  #        ylim=c(0,ymax), col="red")
  #   AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint, cexval=NA, cexaxisval=1, lwdval=1.75, lineval=-0.2, lasval=1)
  #   mtext(expression(paste(italic("F") ~ (year^{-1}))),las=3,side=2,line=2.5,cex=1.0,lwd=1.75)
  # }
  
  plot_data_1 = cbind(midpt, res_1$FAtLen, res_1$ZAtLen) %>% 
    as.data.frame() %>% 
    mutate(Data_Type = "Stereo-BRUV")
  
  plot_data_2 = cbind(midpt, res_2$FAtLen, res_2$ZAtLen) %>% 
    as.data.frame() %>% 
    mutate(Data_Type = "Line fishing")
  
  plot_data <- rbind(plot_data_1, plot_data_2)
  
  ggplot()+
    geom_line(data=plot_data, aes(x=midpt, y=V2, colour=Data_Type, group=Data_Type))+
    geom_line(data=plot_data, aes(x=midpt, y=(V3-V2)), colour="grey40", lty="dashed")+
    scale_colour_manual(values = c("Stereo-BRUV"=colour_1, "Line fishing"=colour_2), name="Data Source")+
    xlab(xaxis_lab)+
    ylab(yaxis_lab)+
    scale_x_continuous(breaks=seq(0,xmax, xint))+
    theme_classic()
  
  #par(.pardefault)
}



