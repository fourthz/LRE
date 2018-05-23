gSIGN.test<-function (x, y = NULL, md = 0, alternative = "two.sided", conf.level = 0.95) 
{
  choices <- c("two.sided", "greater", "less")
  alt <- pmatch(alternative, choices)
  alternative <- choices[alt]
  if (length(alternative) > 1 || is.na(alternative)) 
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
  if (!missing(md)) 
    if (length(md) != 1 || is.na(md)) 
      stop("median must be a single number")
  if (!missing(conf.level)) 
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 
        0 || conf.level > 1) 
      stop("conf.level must be a number between 0 and 1")
  if (is.null(y)) {
    dname <- paste(deparse(substitute(x)))
    x <- sort(x)
    diff <- (x - md)
    n <- length(x)
    nt <- length(x) - sum(diff == 0)
    s <- sum(diff > 0)
    estimate <- median(x)
    method <- c("One-sample Sign-Test")
    names(estimate) <- c("median of x")
    names(md) <- "median"
    names(s) <- "s"
    CIS <- "Conf Intervals"
    if (alternative == "less") {
      pval <- sum(dbinom(0:s, nt, 0.5))
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)][1]
      if (k < 1) {
        conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
        xl <- -Inf
        xu <- x[n]
        ici <- c(xl, xu)
      }
      else {
        ci1 <- c(-Inf, x[n - k + 1])
        acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(-Inf, x[n - k])
        acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
        xl <- -Inf
        xu <- (((x[n - k + 1] - x[n - k]) * (conf.level - 
                                               acl2))/(acl1 - acl2)) + x[n - k]
        ici <- c(xl, xu)
      }
    }
    else if (alternative == "greater") {
      pval <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)][1]
      if (k < 1) {
        conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
        xl <- x[1]
        xu <- Inf
        ici <- c(xl, xu)
      }
      else {
        ci1 <- c(x[k], Inf)
        acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(x[k + 1], Inf)
        acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
        xl <- (((x[k] - x[k + 1]) * (conf.level - acl2))/(acl1 - 
                                                            acl2)) + x[k + 1]
        xu <- Inf
        ici <- c(xl, xu)
      }
    }
    else {
      p1 <- sum(dbinom(0:s, nt, 0.5))
      p2 <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
      pval <- min(2 * p1, 2 * p2, 1)
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)/2][1]
      if (k < 1) {
        conf.level <- (1 - 2 * (sum(dbinom(k, n, 0.5))))
        xl <- x[1]
        xu <- x[n]
        ici <- c(xl, xu)
      }
      else {
        ci1 <- c(x[k], x[n - k + 1])
        acl1 <- (1 - 2 * (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(x[k + 1], x[n - k])
        acl2 <- (1 - 2 * (sum(dbinom(0:k, n, 0.5))))
        xl <- (((x[k] - x[k + 1]) * (conf.level - acl2))/(acl1 - 
                                                            acl2)) + x[k + 1]
        xu <- (((x[n - k + 1] - x[n - k]) * (conf.level - 
                                               acl2))/(acl1 - acl2)) + x[n - k]
        ici <- c(xl, xu)
      }
    }
  }
  else {
    if (length(x) != length(y)) 
      stop("Length of x must equal length of y")
    xy <- sort(x - y)
    diff <- (xy - md)
    n <- length(xy)
    nt <- length(xy) - sum(diff == 0)
    s <- sum(diff > 0)
    dname <- paste(deparse(substitute(x)), " and ", deparse(substitute(y)), 
                   sep = "")
    estimate <- median(xy)
    method <- c("Dependent-samples Sign-Test")
    names(estimate) <- c("median of x-y")
    names(md) <- "median difference"
    names(s) <- "S"
    CIS <- "Conf Intervals"
    if (alternative == "less") {
      pval <- sum(dbinom(0:s, nt, 0.5))
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)][1]
      if (k < 1) {
        conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
        xl <- -Inf
        xu <- xy[n]
        ici <- c(xl, xu)
      }
      else {
        ci1 <- c(-Inf, xy[n - k + 1])
        acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(-Inf, xy[n - k])
        acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
        xl <- -Inf
        xu <- (((xy[n - k + 1] - xy[n - k]) * (conf.level - 
                                                 acl2))/(acl1 - acl2)) + xy[n - k]
        ici <- c(xl, xu)
      }
    }
    else if (alternative == "greater") {
      pval <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)][1]
      if (k < 1) {
        conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
        xl <- xy[1]
        xu <- Inf
        ici <- c(xl, xu)
      }
      else {
        ci1 <- c(xy[k], Inf)
        acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(xy[k + 1], Inf)
        acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
        xl <- (((xy[k] - xy[k + 1]) * (conf.level - acl2))/(acl1 - 
                                                              acl2)) + xy[k + 1]
        xu <- Inf
        ici <- c(xl, xu)
      }
    }
    else {
      p1 <- sum(dbinom(0:s, nt, 0.5))
      p2 <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
      pval <- min(2 * p1, 2 * p2, 1)
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)/2][1]
      if (k < 1) {
        conf.level <- (1 - 2 * (sum(dbinom(k, n, 0.5))))
        xl <- xy[1]
        xu <- xy[n]
        ici <- c(xl, xu)
      }
      else {
        ci1 <- c(xy[k], xy[n - k + 1])
        acl1 <- (1 - 2 * (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(xy[k + 1], xy[n - k])
        acl2 <- (1 - 2 * (sum(dbinom(0:k, n, 0.5))))
        xl <- (((xy[k] - xy[k + 1]) * (conf.level - acl2))/(acl1 - 
                                                              acl2)) + xy[k + 1]
        xu <- (((xy[n - k + 1] - xy[n - k]) * (conf.level - 
                                                 acl2))/(acl1 - acl2)) + xy[n - k]
        ici <- c(xl, xu)
      }
    }
  }
  if (k < 1) {
    cint <- ici
    attr(cint, "conf.level") <- conf.level
    rval <- structure(list(statistic = s, p.value = pval, 
                           estimate = estimate, null.value = md, alternative = alternative, 
                           method = method, data.name = dname, conf.int = cint))
    oldClass(rval) <- "htest"
    return(rval)
  }
  else {
    result1 <- c(acl2, ci2)
    result2 <- c(conf.level, ici)
    result3 <- c(acl1, ci1)
    Confidence.Intervals <- round(as.matrix(rbind(result1, 
                                                  result2, result3)), 4)
    cnames <- c("Conf.Level", "L.E.pt", "U.E.pt")
    rnames <- c("Lower Achieved CI", "Interpolated CI", "Upper Achieved CI")
    dimnames(Confidence.Intervals) <- list(rnames, cnames)
    cint <- ici
    attr(cint, "conf.level") <- conf.level
    rval <- structure(list(statistic = s, parameter = NULL, 
                           p.value = pval, conf.int = cint, estimate = estimate, 
                           null.value = md, alternative = alternative, method = method, 
                           data.name = dname))
    oldClass(rval) <- "htest"
    rval$Confidence.Intervals <- Confidence.Intervals
    return(rval)
  }
}



onesam.power<-function(dist, test, n, parm1, parm2, draws){
  draws <- draws
  n <- n
  if(dist=="normal"){
    mu <- parm1*parm2
    sigma <- parm2
    samp <- matrix(rnorm(n * draws, mu, sigma), byrow = draws, ncol = n)
  }
  if(dist=="logistic"){
    location <- parm1*parm2
    scale <- parm2*.5513288954
    samp <- matrix(rlogis(n * draws, location, scale), byrow = draws, ncol = n)
  }  
  if(dist=="dexp"){
    require(smoothmest)
    mu <- parm1*parm2
    lambda <- parm2/sqrt(2)
    samp <- matrix(rdoublex(n * draws, mu, lambda), byrow = draws, ncol = n)
  } 
  if(dist=="uniform"){      
    mu <- parm1*parm2
    temp <- ((parm2/sqrt(1/12)*.5)*-1)
    max <- -temp+mu
    min <- temp+mu
    samp <- matrix(runif(n * draws, min, max), byrow = draws, ncol = n)
  }  
  if(test=="sign"){
    p.vals<-apply(samp, 1, function(x)
      gSIGN.test(x, md = 0, alternative = "two.sided", conf.level = 0.95)$p.value)
    power <- sum(p.vals < .05)/draws
    return(power)
  }
  if(test=="t"){
    # The apply command is used to apply a function to rows or columns of a matrix.
    # Here I'm saying to apply function x to the sample matrix (samp) to each row (1).
    p.vals<-apply(samp, 1, function(x)
      t.test(x, mu=0, alternative = "two.sided")$p.value)
    power <- sum(p.vals < .05)/draws
    return(power)}
  if(test=="vdw"){
    require(snpar)
    p.vals<-apply(samp, 1, function(x)
      ns.test(x, q=0)$p.value)
    power <- sum(p.vals < .05)/draws
    return(power)}  
  if(test=="wilcoxon"){
    p.vals<-apply(samp, 1, function(x)
      wilcox.test(x, mu=0)$p.value)
    power <- sum(p.vals < .05)/draws
    return(power)}
} 




find.n<-function(x)abs(onesam.power(dist=distr, test=proc, n=x, parm1=es, parm2=parm_two, draws=10000)-.8)


distr<-"dexp"
proc<-"sign"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(1450,1900), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE)

distr<-"dexp"
proc<-"sign"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(420,530), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)






distr<-"dexp"
proc<-"sign"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(190,245), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(125,145), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(84,116), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(64,80), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(49,60), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(37,48), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(30,46), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(28,39), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(20,37), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(15,32), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(12,29), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(12,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(12,23), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(10,24), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(10,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(7,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(7,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"sign"
es<-1
parm_two<-1
result<-optimize(find.n, c(7,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(2750,3450), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(690,920), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(320,385), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(180,220), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(110,154), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(84,100), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(58,73), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(44,58), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(31,50), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(25,41), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(20,36), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(17,32), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(14,28), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(12,23), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(10,22), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(9,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(7,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(6,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(5,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"t"
es<-1
parm_two<-1
result<-optimize(find.n, c(4,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(2250,2800), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(570,710), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(250,325), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(150,175), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(94,120), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(70,86), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(51,61), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(40,50), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(30,45), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(23,37), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(19,34), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(15,29), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(12,26), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(12,24), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(10,22), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(9,20), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(7,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(7,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(6,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"vdw"
es<-1
parm_two<-1
result<-optimize(find.n, c(5,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(1900,2300), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(480,580), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(215,265), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(125,155), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(80,106), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(62,76), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(41,58), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(32,48), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(25,40), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(20,35), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(17,31), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(15,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(10,26), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(11,22), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(8,22), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(8,20), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(7,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(6,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(5,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"dexp"
proc<-"wilcoxon"
es<-1
parm_two<-1
result<-optimize(find.n, c(5,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(3300,4550), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(880,1170), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(375,500), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(225,290), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(150,178), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(110,140), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(79,102), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(58,76), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(44,62), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(37,55), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(30,48), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(25,41), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(20,37), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(18,37), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(15,29), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(12,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(12,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(10,24), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(10,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"sign"
es<-1
parm_two<-1
result<-optimize(find.n, c(7,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(2750,3550), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(710,880), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(310,395), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(180,220), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(110,142), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(84,100), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(60,72), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(47,59), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(33,50), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(26,42), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(21,36), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(17,31), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(13,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(12,24), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(10,23), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(8,20), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(7,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(6,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(5,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"t"
es<-1
parm_two<-1
result<-optimize(find.n, c(5,15), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(2650,3500), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(640,830), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(295,385), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(170,200), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(108,134), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(82,94), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(57,71), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(45,55), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(32,49), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(26,41), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(21,36), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(17,31), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(14,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(13,26), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(11,22), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(9,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(7,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(7,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(6,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"vdw"
es<-1
parm_two<-1
result<-optimize(find.n, c(5,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(2550,3150), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(630,830), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(280,355), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(155,205), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(106,132), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(74,90), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(56,68), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(39,57), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(30,48), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(24,39), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(19,34), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(17,30), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(15,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(12,24), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(10,22), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(9,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(8,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(7,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(6,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"logistic"
proc<-"wilcoxon"
es<-1
parm_two<-1
result<-optimize(find.n, c(5,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(4450,5650), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(1090,1480), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(520,640), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(295,350), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(186,232), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(132,168), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(94,122), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(77,96), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(63,80), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(49,62), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(42,53), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(35,46), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(25,44), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(23,39), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(20,37), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(18,32), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(15,29), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(12,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(13,27), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"sign"
es<-1
parm_two<-1
result<-optimize(find.n, c(10,24), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(2850,3550), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(700,890), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(320,390), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(185,225), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(116,140), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(82,100), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(59,74), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(45,55), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(34,50), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(27,42), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(22,37), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(18,30), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(15,28), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(12,25), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(10,22), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(9,20), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(7,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(6,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(6,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"t"
es<-1
parm_two<-1
result<-optimize(find.n, c(5,15), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(2900,3700), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(700,920), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(310,400), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(175,235), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(120,140), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(84,98), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(59,74), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(47,58), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(35,50), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(28,41), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(21,35), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(18,32), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(14,29), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(12,25), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(10,23), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(9,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(8,20), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(6,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(6,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"vdw"
es<-1
parm_two<-1
result<-optimize(find.n, c(6,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(3000,3750), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(720,910), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(330,410), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(195,220), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(120,158), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(86,102), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(62,79), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(50,61), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(35,51), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(27,45), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(21,38), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(18,32), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(15,28), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(13,26), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(11,23), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(9,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(8,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(7,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(6,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"normal"
proc<-"wilcoxon"
es<-1
parm_two<-1
result<-optimize(find.n, c(5,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(8300,10700), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(2200,2600), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(980,1140), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(550,680), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(350,440), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(250,305), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(180,225), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(140,170), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(110,135), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(88,106), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(72,96), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(58,78), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(54,64), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(44,57), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(37,48), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(30,46), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(25,44), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(20,39), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(18,37), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"sign"
es<-1
parm_two<-1
result<-optimize(find.n, c(18,32), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(2800,3600), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(750,900), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(320,380), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(190,220), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(125,135), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(85,95), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(60,75), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(45,60), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(35,50), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(27,39), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(23,35), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(17,31), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(15,25), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(12,24), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(10,22), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(9,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(8,20), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(6,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(6,16), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"t"
es<-1
parm_two<-1
result<-optimize(find.n, c(5,15), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(1400,1600), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(450,500), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(220,240), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(140,150), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(90,100), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(70,80), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(50,60), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(40,50), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(35,45), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(27,41), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(23,35), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(19,31), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(15,29), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(13,26), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(12,23), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(10,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(9,20), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(8,19), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(7,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"vdw"
es<-1
parm_two<-1
result<-optimize(find.n, c(6,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.05
parm_two<-1
result<-optimize(find.n, c(3100,3500), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.1
parm_two<-1
result<-optimize(find.n, c(750,950), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.15
parm_two<-1
result<-optimize(find.n, c(360,400), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.2
parm_two<-1
result<-optimize(find.n, c(200,240), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.25
parm_two<-1
result<-optimize(find.n, c(140,150), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.3
parm_two<-1
result<-optimize(find.n, c(95,110), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.35
parm_two<-1
result<-optimize(find.n, c(70,85), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.4
parm_two<-1
result<-optimize(find.n, c(55,65), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.45
parm_two<-1
result<-optimize(find.n, c(40,55), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.5
parm_two<-1
result<-optimize(find.n, c(33,47), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.55
parm_two<-1
result<-optimize(find.n, c(25,39), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.6
parm_two<-1
result<-optimize(find.n, c(23,35), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.65
parm_two<-1
result<-optimize(find.n, c(19,31), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.7
parm_two<-1
result<-optimize(find.n, c(15,28), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.75
parm_two<-1
result<-optimize(find.n, c(13,26), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.8
parm_two<-1
result<-optimize(find.n, c(12,23), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.85
parm_two<-1
result<-optimize(find.n, c(10,21), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.9
parm_two<-1
result<-optimize(find.n, c(8,20), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-0.95
parm_two<-1
result<-optimize(find.n, c(8,18), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
distr<-"uniform"
proc<-"wilcoxon"
es<-1
parm_two<-1
result<-optimize(find.n, c(7,17), tol = 0.0001)
dat<-as.data.frame(t(as.matrix(unlist(result))))
dat$dist<-distr
dat$test<-proc
dat$es<-es
write.table(dat, "C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\May 2018 - Rerunning Optimizer\\find_sample_sizes_21May2018.csv", append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)










onesam.power<-function(dist, n, parm1, parm2, draws){
  draws <- draws
  n <- n
  if(dist=="normal"){
    mu <- parm1*parm2
    sigma <- parm2
    samp <- matrix(rnorm(n * draws, mu, sigma), byrow = draws, ncol = n)
  }
  if(dist=="logistic"){
    location <- parm1*parm2
    scale <- parm2*.5513288954
    samp <- matrix(rlogis(n * draws, location, scale), byrow = draws, ncol = n)
  }  
  if(dist=="dexp"){
    require(smoothmest)
    mu <- parm1*parm2
    lambda <- parm2/sqrt(2)
    samp <- matrix(rdoublex(n * draws, mu, lambda), byrow = draws, ncol = n)
  } 
  if(dist=="uniform"){      
    mu <- parm1*parm2
    temp <- ((parm2/sqrt(1/12)*.5)*-1)
    max <- -temp+mu
    min <- temp+mu
    samp <- matrix(runif(n * draws, min, max), byrow = draws, ncol = n)
  }  
  
  p.vals<-apply(samp, 1, function(x)
    gSIGN.test(x, md = 0, alternative = "two.sided", conf.level = 0.95)$p.value)
  signpower <- sum(p.vals < .05)/draws
  #return(signpower)
  
  
  # The apply command is used to apply a function to rows or columns of a matrix.
  # Here I'm saying to apply function x to the sample matrix (samp) to each row (1).
  p.vals<-apply(samp, 1, function(x)
    t.test(x, mu=0, alternative = "two.sided")$p.value)
  tpower <- sum(p.vals < .05)/draws
  #return(power)
  
  
  require(snpar)
  p.vals<-apply(samp, 1, function(x)
    ns.test(x, q=0)$p.value)
  vdwpower <- sum(p.vals < .05)/draws
  #return(power)
  
  
  p.vals<-apply(samp, 1, function(x)
    wilcox.test(x, mu=0)$p.value)
  wilcoxonpower <- sum(p.vals < .05)/draws
  #return(power)
  
  allpower<-cbind(dist, n, parm1, parm2, signpower, tpower, vdwpower, wilcoxonpower)
  colnames(allpower)<-c("Dist", "Sample_Size", "ES", "SD", "Sign", "t", "VDW", "Wilc")
  return(allpower)
  #write.table(allpower, paste(c("C:\\Users\\grant_morgan\\Box Sync\\Research\\Mike Seaman\\allpower.csv")), append=TRUE, sep=",", row.names = FALSE)
} 

onesam.power("dexp", n=1728, parm1 = .05, parm2 = 1, draws = 100000)




