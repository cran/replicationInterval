#' A common problem faced by journal reviewers and authors is the question of
#' whether the results of a replication study are consistent with the original
#' published study. One solution to this problem is to examine the effect size
#' from the published study and generate the range of effect sizes that could
#' reasonably be obtained (due to random sampling) in a replication attempt
#' (i.e., calculate a replication interval). If a replication effect size falls
#' outside the replication interval then that value could not have occurred to
#' due the effects of sampling error alone. Alternatively, if a replication
#' effect size falls within the replication interval then the replication
#' results could have reasonably occurred due to the effects of sampling error
#' alone. This package calculates the replication interval for two different types of effect sizes (i.e., correlation, Cohen's \emph{d} ).
#'\tabular{ll}{
#'Package: \tab replicationInterval\cr
#'Type: \tab Package\cr
#'Version: \tab 0.1\cr
#'Date: \tab 2014-07-19\cr
#'License: \tab Unlimited\cr
#'}
#'\code{\link{ri.d}} creates a replication interval for a Cohen's \emph{d} \cr
#'\code{\link{ri.r}} creates a replication interval for a correlation\cr
#'
#'@name replicationInterval-package
#'@aliases replicationInterval
#'@docType package
#'@title Replication Interval 
#'@author 
#'\tabular{ll}{
#'Author: \tab David J. Stanley \email{dstanley@@uoguelph.ca}\cr
#'Maintainer: \tab David J. Stanley \email{dstanley@@uoguelph.ca}
#'}
#'@references 
#'Spence, J.R., & Stanley, D.J.(2014). Replication Interval: What to expect when you're expecting ... a replication. [manuscript in preparation] \cr\cr
#'Stanley, D.J., & Spence, J.R.(2014). Expectations for replications: Are yours realistic? \emph{Perspectives on Psychological Science, 9}, 305-318.
#'@keywords package
#'@examples
#' ri.d(d=.65,n1=20,n2=20)
#' ri.d(d=.65,n1=20,n2=20,explain=TRUE)
#' ri.r(r=.35,n=40)
#' ri.r(r=.35,n=40,explain=TRUE)
NULL



quantile.d <- function (probability.in.tail,ncd,n1,n2,is.lower.tail) {
     nct.value=ncd/sqrt((1/n1)+(1/n2))
     df.t=n1+n2-2
     t.tail=qt(p=probability.in.tail, ncp=nct.value,df=df.t,lower.tail=is.lower.tail)
     d.tail=t.tail*sqrt((1/n1)+(1/n2))
     return(d.tail)
}


ri.d.output <- function(interval.output) {
     conf.level.percent=round(interval.output$conf.level*100,0)
     
     d <- interval.output$d
     lower.conf.limit <- interval.output$lower.conf.limit.d 
     upper.conf.limit <- interval.output$upper.conf.limit.d
     n1 <- interval.output$conf.limit.n1
     n2 <- interval.output$conf.limit.n2
     
     lower.population.lower.bound <- interval.output$lower.population.lower.bound.d
     lower.population.upper.bound <- interval.output$lower.population.upper.bound.d 
     
     upper.population.lower.bound <- interval.output$upper.population.lower.bound.d
     upper.population.upper.bound <- interval.output$upper.population.upper.bound.d 
     
     lower.replication.interval <- interval.output$lower.replication.interval.d
     upper.replication.interval <- interval.output$upper.replication.interval.d
     
     n1.replication <- interval.output$replication.interval.n1
     n2.replication <- interval.output$replication.interval.n2
     
     text.ci <- sprintf("Due to sampling error, the published/sample d = %1.2f (N1 = %d and N2 = %d) could have been created by a population d-value as low as %1.2f (lower-population d-value) or as high as %1.2f (upper-population d-value).",d,n1,n2,lower.conf.limit,upper.conf.limit)
     
     text.lower.population <- sprintf("If the lower-population d-value of %1.2f created the published d = %1.2f (due to sampling error) then %d%% of replication d-values (using N1 = %d and N2 = % d) will fall between %1.2f and %1.2f.",lower.conf.limit,d,conf.level.percent,n1.replication,n2.replication,lower.population.lower.bound,lower.population.upper.bound)
     
     
     text.upper.population <- sprintf("If the upper-population d-value of %1.2f created the published d = %1.2f (due to sampling error) then %d%% of replication d-values (using N1 = %d and N2 = % d) will fall between %1.2f and %1.2f.",upper.conf.limit,d,conf.level.percent,n1.replication,n2.replication,upper.population.lower.bound,upper.population.upper.bound)
     
     
     text.summary <- sprintf("That is, if the study was replicated (using N1 = %d and N2 = %d) the researcher could obtain any d-value in the %1.2f to %1.2f interval due to sampling error.",n1.replication,n2.replication,lower.replication.interval,upper.replication.interval)
     
     text.replication.interval.brief <- sprintf("Replication Interval for d = %1.2f is {%1.2f, %1.2f}.", d,lower.replication.interval,upper.replication.interval)
     
     text.replication.interval <- sprintf("Thus, for d = %1.2f (when replications are conducted with N1 = %d and N2 = %d) the replication interval is {%1.2f, %1.2f}.", d,n1.replication,n2.replication,lower.replication.interval,upper.replication.interval)
     
     
     text.explain <- paste("Explanation:",text.ci,text.lower.population,text.upper.population, text.summary,text.replication.interval)
     
     
     text.out <- list()
     text.out$replication.interval <- text.replication.interval.brief
     text.out$replication.interval.explanation <- text.explain
     return(text.out)
}


#' Creates a replication interval based on a published/sample Cohen's \emph{d} 
#' @param d A sample Cohen's \emph{d}-value (standardized mean difference) created with pooled variance denominator.
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param conf.level (optional 0 to 1 value) Confidence level desired (0 to 1). If not specified .95 (i.e, 95 percent) will be used.
#' @param explain (optional boolean) Default is FALSE. If TRUE, text output explaining the interval is provided. 
#' @param extended.output (optional boolean) Default is FALSE. If TRUE, additional details (e.g., confidence interval) provided in numeric return output.
#' @return A list of values (\code{lower.replication.interval.d, upper.replication.interval.d}) containing the replication interval (and related statistics if requested with the \code{extended.output} argument).
#' @examples
#' ri.d(d=.65,n1=20,n2=20)
#' ri.d(d=.65,n1=20,n2=20,explain=TRUE)
#' @export 
ri.d <- function (d,n1,n2,conf.level=.95,explain=FALSE,extended.output=FALSE) {
     
     d.value.ci=MBESS::ci.smd(smd=d,n.1=n1,n.2=n2,conf.level=conf.level)
     
     d <- d.value.ci$smd
     
     conf.lower.bound <- d.value.ci$Lower.Conf.Limit.smd
     conf.upper.bound <- d.value.ci$Upper.Conf.Limit.smd
     
     probability.in.tail=(1 - conf.level)/2
     
     
     lower.population.lower.bound <- quantile.d(probability.in.tail=probability.in.tail, ncd=conf.lower.bound,n1=n1,n2=n2,is.lower.tail=TRUE)
     lower.population.upper.bound <- quantile.d(probability.in.tail=probability.in.tail, ncd=conf.lower.bound,n1=n1,n2=n2,is.lower.tail=FALSE)
     
     upper.population.lower.bound <- quantile.d(probability.in.tail=probability.in.tail, ncd=conf.upper.bound,n1=n1,n2=n2,is.lower.tail=TRUE)
     upper.population.upper.bound <- quantile.d(probability.in.tail=probability.in.tail, ncd=conf.upper.bound,n1=n1,n2=n2,is.lower.tail=FALSE)
     
     interval.output <- list()
     interval.output$conf.level <- conf.level
     interval.output$d <- d
     interval.output$lower.conf.limit.d <- conf.lower.bound
     interval.output$upper.conf.limit.d <- conf.upper.bound
     interval.output$conf.limit.n1 <- n1
     interval.output$conf.limit.n2 <- n2
     interval.output$lower.population.lower.bound.d <- lower.population.lower.bound
     interval.output$lower.population.upper.bound.d <- lower.population.upper.bound
     interval.output$upper.population.lower.bound.d <- upper.population.lower.bound
     interval.output$upper.population.upper.bound.d <- upper.population.upper.bound
     interval.output$replication.interval.n1 <- n1
     interval.output$replication.interval.n2 <- n2
     interval.output$lower.replication.interval.d <- lower.population.lower.bound
     interval.output$upper.replication.interval.d <- upper.population.upper.bound

     interval.output.brief <- list()
     interval.output.brief$lower.replication.interval.d <- lower.population.lower.bound
     interval.output.brief$upper.replication.interval.d <- upper.population.upper.bound
     
     if (extended.output==TRUE) {
          interval.results <- interval.output
     } else {
          interval.results <- interval.output.brief
     }
     
     if (explain==TRUE) {
          cat("\n")
          rText<-ri.d.output(interval.output)
          cat(rText[[1]])
          cat("\n\n")
          cat(rText[[2]])
          cat("\n\n\n")
          
     }

     return(interval.results)
}


















quantile.r <- function (probability.in.tail,ncr,n,is.lower.tail) {
     ncz <- atanh(ncr)
     ncz.se <- 1 / sqrt(n-3)
     z.tail <- qnorm(p=probability.in.tail, mean=ncz, sd=ncz.se,lower.tail=is.lower.tail)
     r.tail <- tanh(z.tail)
     return(r.tail)
}

ci.r <- function(r,n,conf.level) {
     probability.in.tail <- (1- conf.level)/2
     
     obs.z <- atanh(r)
     obs.z.se <- 1/sqrt(n-3)
     
     high.z <- qnorm(p=probability.in.tail, mean=obs.z,sd=obs.z.se,lower.tail=FALSE)
     high.r <- tanh(high.z)
     
     low.z <- qnorm(p=probability.in.tail, mean=obs.z,sd=obs.z.se,lower.tail=TRUE)
     low.r <- tanh(low.z)
     
     
     conf.interval.output <- list()
     
     conf.interval.output$r <- r
     conf.interval.output$lower.conf.limit.r <- low.r
     conf.interval.output$upper.conf.limit.r <- high.r
     return(conf.interval.output)
}

ri.r.output <- function(interval.output) {
     conf.level.percent=round(interval.output$conf.level*100,0)
     
     r <- interval.output$r
     lower.conf.limit <- interval.output$lower.conf.limit.r 
     upper.conf.limit <- interval.output$upper.conf.limit.r
     n <- interval.output$upper.conf.limit.n
     
     lower.population.lower.bound <- interval.output$lower.population.lower.bound.r
     lower.population.upper.bound <- interval.output$lower.population.upper.bound.r 
     
     upper.population.lower.bound <- interval.output$upper.population.lower.bound.r
     upper.population.upper.bound <- interval.output$upper.population.upper.bound.r 
     
     lower.replication.interval <- interval.output$lower.replication.interval.r
     upper.replication.interval <- interval.output$upper.replication.interval.r
     
     n.replication <- interval.output$replication.interval.n
     
     text.ci <- sprintf("Due to sampling error, the published/sample r = %1.2f (N = %d) could have been created by a population correlation as low as %1.2f (lower-population correlation) or as high as %1.2f (upper-population correlation).",r,n,lower.conf.limit,upper.conf.limit)
     
     text.lower.population <- sprintf("If the lower-population correlation of %1.2f created the published r = %1.2f (due to sampling error) then %d%% of replication correlations (using N = %d) will fall between %1.2f and %1.2f.",lower.conf.limit,r,conf.level.percent,n.replication,lower.population.lower.bound,lower.population.upper.bound)
     
     
     text.upper.population <- sprintf("If the upper-population correlation of %1.2f created the published r = %1.2f (due to sampling error) then %d%% of replication correlations (using N = %d) will fall between %1.2f and %1.2f.",upper.conf.limit,r,conf.level.percent,n.replication,upper.population.lower.bound,upper.population.upper.bound)
     
     
     text.summary <- sprintf("That is, if the study was replicated (using N = %d) the researcher could obtain any correlation in the %1.2f to %1.2f interval due to sampling error.",n.replication,lower.replication.interval,upper.replication.interval)
     
     text.replication.interval.brief <- sprintf("Replication Interval for r = %1.2f is {%1.2f, %1.2f}.", r,lower.replication.interval,upper.replication.interval)
     
     text.replication.interval <- sprintf("Thus, for r = %1.2f (when replications are conducted with N = %d) the replication interval is {%1.2f, %1.2f}.", r,n.replication,lower.replication.interval,upper.replication.interval)
     
     
     text.explain <- paste("Explanation:",text.ci,text.lower.population,text.upper.population, text.replication.interval)
     
     text.out <- list()
     text.out$replication.interval <-text.replication.interval.brief
     text.out$replication.interval.explanation <- text.explain
     return(text.out)
}


#' Creates a replication interval based on a published/sample correlation. 
#' @param r A sample correlation
#' @param n Sample size 
#' @param conf.level (optional 0 to 1 value) Confidence level desired (0 to 1). If not specified .95 (i.e, 95 percent) will be used.
#' @param explain (optional boolean) Default is FALSE. If TRUE, text output explaining the interval is provided. 
#' @param extended.output (optional boolean) Default is FALSE. If TRUE, additional details (e.g., confidence interval) provided in numeric return output.
#' @return A list of values (\code{lower.replication.interval.r, upper.replication.interval.r}) containing the replication interval (and related statistics if requested with the \code{extended.output} argument).
#' @examples
#' ri.r(r=.35,n=40)
#' ri.r(r=.35,n=40,explain=TRUE)
#' @export
ri.r <- function (r,n,conf.level = .95,explain=FALSE,extended.output=FALSE) {

     
     r.value.ci=ci.r(r=r,n=n,conf.level=conf.level)
     r <- r.value.ci$r
     
     
     conf.lower.bound <- r.value.ci$lower.conf.limit.r
     conf.upper.bound <- r.value.ci$upper.conf.limit.r
     
     probability.in.tail=(1 - conf.level)/2
     
     lower.population.lower.bound <- quantile.r(probability.in.tail=probability.in.tail, ncr=conf.lower.bound,n=n,is.lower.tail=TRUE)
     lower.population.upper.bound <- quantile.r(probability.in.tail=probability.in.tail, ncr=conf.lower.bound,n=n,is.lower.tail=FALSE)
     
     
     upper.population.lower.bound <- quantile.r(probability.in.tail=probability.in.tail, ncr=conf.upper.bound,n=n,is.lower.tail=TRUE)
     upper.population.upper.bound <- quantile.r(probability.in.tail=probability.in.tail, ncr=conf.upper.bound,n=n,is.lower.tail=FALSE)
     
     interval.output <- list()
     interval.output$conf.level <- conf.level
     
     interval.output$r <- r
     interval.output$lower.conf.limit.r <- conf.lower.bound
     interval.output$upper.conf.limit.r <- conf.upper.bound
     interval.output$upper.conf.limit.n <- n
     
     interval.output$lower.population.lower.bound.r <- lower.population.lower.bound
     interval.output$lower.population.upper.bound.r <- lower.population.upper.bound
     
     interval.output$upper.population.lower.bound.r <- upper.population.lower.bound
     interval.output$upper.population.upper.bound.r <- upper.population.upper.bound
     
     interval.output$replication.interval.n <- n
     
     interval.output$lower.replication.interval.r <- lower.population.lower.bound
     interval.output$upper.replication.interval.r <- upper.population.upper.bound

     
     interval.output.brief <- list()
     interval.output.brief$lower.replication.interval.r <- lower.population.lower.bound
     interval.output.brief$upper.replication.interval.r <- upper.population.upper.bound
     
     if (extended.output==TRUE) {
          interval.results <- interval.output
     } else {
          interval.results <- interval.output.brief
     }
     
     if (explain==TRUE) {
          rText<-ri.r.output(interval.output)
          cat("\n")
          cat(rText[[1]])
          cat("\n\n")
          cat(rText[[2]])
          cat("\n\n\n")
          
     }
     
     
     return(interval.results)
}









