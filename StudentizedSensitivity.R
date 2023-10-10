#' StudentizedSensitivity
#'Function to perform a Studentized Sensitivity analysis on the sample average treatment
#'effect in a paired observational study
#'
#' @param PairedDiff: Vector of treated-minus-control paired differences.
#' @param null: Value of the sample average treatment effect under the null.
#' @param alpha: Desired Type I error rate.
#' @param alternative: Can be "less", "greater", or "two.sided".
#' @param Gamma: Vector of values for Gamma at which to perform the sensitivity
#'  analysis.
#' @param nperm: Number of permutations to perform permutation test.
#' @param Changepoint: If true, function returns the maximal Gamma for which the
#' test rejects at level alpha.
#' @param SensitivityInterval: If true, function returns (100-alpha) sensitivity
#' intervals. They will be one-sided if the alternative is less than or greater than,
#' and two-sided if the alternative is two-sided.
#'
#' @return Gamma: Vector of Gammas for which the sensitivity analysis was performed.
#' @return pval: P-values for each value of Gamma.
#' @return GammaPval: Matrix combining Gamma and pval.
#' @return Changepoint: Maximal Gamma for which the test rejected at level alpha.
#' @return SensitivityInterval: Upper and lower bounds for 100(1-alpha) sensitivity
#' intervals for each value of Gamma.
#'
#' @examples
#' #Reading in Data
#' data <- read.csv(file="/data/maffei.csv")
#' treatment <- data$Treatment
#' PairedDiff <- data$CNPlus[1:20]-data$CNPlus[21:40]
#' #Function call
#' StudentizedSensitivity(PairedDiff, null = 0, alpha = 0.05, alternative = "greater", Gamma = 1, nperm = 50000, Changepoint = T, SensitivityInterval = T)
#' @export


StudentizedSensitivity = function(PairedDiff, null = 0, alpha = 0.05, alternative = "greater", Gamma = 1, nperm = 50000, Changepoint = T, SensitivityInterval = T)
{
   if(any(Gamma < 1))
   {
   	stop("Values for Gamma must be >= 1")
   }
   if(alternative!="less" & alternative!= "greater" & alternative != "two.sided")
   {
   	stop("Values for alternative are `less', `greater', or `two.sided'")
   }
   if(length(null) > 1)
   {
   	stop("Value under the null must be a scalar")
   }
      if(alpha < 0 | alpha > 0.5)
   {
   	stop("alpha must be between 0 and 0.5")
      }

  PairedDifftrue <- PairedDiff
  alphatrue <- alpha
  I <- length(PairedDiff)
  Adjust <- PairedDiff - null

  if(alternative == "less")
  {
  	Adjust <- -Adjust
  }
  if(alternative == "two.sided")
  {
  	alpha <- alphatrue/2

  	if(mean(Adjust) < 0)
  	{
  		Adjust <- -Adjust
  	}
  }

  pval <- rep(0, length(Gamma))

  for(i in 1:length(Gamma))

  {
  D <- (Adjust) - (Gamma[i]-1)/(1+Gamma[i])*abs(Adjust)
  obs <- mean(D)/(sd(D)/sqrt(I))
  Adjmat <- matrix(abs(Adjust), I, nperm)
  Zmat <- matrix(runif(I*nperm) < Gamma[i]/(1+Gamma[i]), I, nperm)
  Dmat <- (2*Zmat-1)*(Adjmat) - (Gamma[i]-1)/(1+Gamma[i])*Adjmat
  perm <- colMeans(Dmat)/(sqrt(colVars(Dmat)/I))
  pval[i] <- (1+sum(perm>=obs))/(nperm + 1)
  }
  pvalret = pval
  if(alternative == "two.sided")
  {
    pvalret = 2*pval
  }
  Pmatrix <- cbind(Gamma, pvalret)
  colnames(Pmatrix) <- c("Gamma", "P-value")

  if(Changepoint == T)
  {
  	proceed <- StudentizedSensitivity(PairedDifftrue, null, alphatrue, alternative, Gamma=1, nperm,
  	                                  Changepoint = F, SensitivityInterval = F)$pval <= alphatrue

  	change <- 1

  	if(proceed)
  	{
  		change <- uniroot(StudentizedChangepoint, interval = c(1, 30), PairedDiff = PairedDifftrue, null = null,
  		                  alpha = alphatrue, alternative = alternative, nperm = nperm,
  		                  extendInt = "upX")$root
  	}
  }

  if(SensitivityInterval == T)
  	{
  		lb = rep(-Inf, length(Gamma))
  		ub = rep(Inf, length(Gamma))
  		for(i in 1:length(Gamma))
  		{
  		  # Warm Starts
  		UB = uniroot(BoundFinder, PairedDifftrue, Gamma[i],
  		             interval = c(mean(PairedDifftrue), mean(PairedDifftrue)+4*sd(PairedDifftrue)/sqrt(I)), extendInt = "yes")$root
  		LB = -uniroot(BoundFinder, -PairedDifftrue, Gamma[i],
  		              interval = c(-mean(PairedDifftrue)-4*sd(PairedDifftrue)/sqrt(I), -mean(PairedDifftrue)), extendInt = "yes")$root

  		SUB = Inf
  		SLB = -Inf

  		if(alternative == "greater")
  		{
  			SLB = uniroot(StudentizedSI, interval = c(UB-4*sd(PairedDifftrue)/sqrt(I), UB), extendInt = "yes",
  			              Gamma = Gamma[i], PairedDiff=PairedDifftrue, alternative = "greater", alpha = alpha, nperm = nperm)$root
  		}

  		if(alternative == "less")
  		{
  			SUB = uniroot(StudentizedSI, interval = c(LB, LB + 4*sd(PairedDifftrue)/sqrt(I)), extendInt = "yes",
  			              Gamma = Gamma[i], PairedDiff=PairedDifftrue, alternative = "less", alpha = alpha, nperm = nperm)$root
  		}

  		if(alternative == "two.sided")
  		{
  		 SLB = uniroot(StudentizedSI, interval = c(UB-4*sd(PairedDifftrue)/sqrt(I), UB), extendInt = "yes",
  		               Gamma = Gamma[i], PairedDiff=PairedDifftrue, alternative = "greater", alpha = alpha, nperm = nperm)$root
		   SUB = uniroot(StudentizedSI, interval = c(LB, LB+4*sd(PairedDifftrue)/sqrt(I)), extendInt = "yes",
		                 Gamma = Gamma[i], PairedDiff=PairedDifftrue, alternative = "less", alpha = alpha, nperm = nperm)$root
  		}

  		lb[i] = SLB
  		ub[i] = SUB
  		}

  	SImat = cbind(Gamma, lb, ub)
  	colnames(SImat) = c("Gamma", "Lower Bound", "Upper Bound")
  	}
  	if(Changepoint == F & SensitivityInterval == F)
  	{
  		return(list(Gamma=Gamma, pval = pvalret, GammaPval = Pmatrix))
  	}
  	if(Changepoint == F & SensitivityInterval == T)
  	{
  		return(list(Gamma = Gamma, pval = pvalret, GammaPval = Pmatrix, SensitivityInterval = SImat))
  	}
  	if(Changepoint == T & SensitivityInterval == F)
  	{
  		return(list(Gamma = Gamma, pval = pvalret, GammaPval = Pmatrix, Changepoint = change))
  	}
  	if(Changepoint == T & SensitivityInterval == T)
  	{
  		return(list(Gamma = Gamma, pval = pvalret, GammaPval = Pmatrix, Changepoint = change,
  		            SensitivityInterval = SImat))
  	}

}


####These are auxiliary functions used for root finding and for calculating columnwise variances in StudentizedSensitivity
StudentizedChangepoint = function(Gamma, PairedDiff, null, alternative, alpha, nperm)
{
  alphachange = alpha
  StudentizedSensitivity(PairedDiff, null, alpha, alternative, Gamma, nperm, Changepoint = F, SensitivityInterval = F)$pval - alphachange
}

StudentizedSI = function(null,  Gamma, PairedDiff,  alternative, alpha, nperm)
{
	StudentizedSensitivity(PairedDiff, null, alpha, alternative, Gamma, nperm, Changepoint = F, SensitivityInterval = F)$pval - alpha
}

BoundFinder = function(null,  PairedDiff, Gamma)
{
	mean(PairedDiff - null - (Gamma-1)/(1+Gamma)*abs(PairedDiff-null))
}

colVars <- function(x) {
  N = nrow(x)
  (colSums(x^2) - colSums(x)^2/N) / (N-1)
}


