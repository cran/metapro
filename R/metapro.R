#' @title Beta probability
#' @param p p-value
#' @param i rank
#' @param n The number of inputs
#' @importFrom "stats" "pbeta" "qbeta"
F_i = function(p, i, n)
{
  a = i
  b = n-i+1
  res = pbeta(q = p, shape1 = a, shape2 = b, lower.tail = T)
  return(res)
}


#' @title wFisher
#' @description sample size-weighted Fisher's method
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weight or sample size for each experiment
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effects, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @return p : Combined p-value
#' @return overall.eff.direction : The direction of combined effects.
#' @examples wFisher(p=c(0.01,0.2,0.8), weight = c(50,60,100),is.onetail=FALSE, eff.sign=c(1,1,1))
#' @importFrom "stats" "pgamma" "qgamma"
#' @references Becker BJ (1994). “Combining significance levels.” In Cooper H, Hedges LV (eds.), A handbook of research synthesis, 215–230. Russell Sage, New York.
#' @references Fisher RA (1925). Statistical methods for research workers. Oliver and Boyd, Edinburgh.
#' @export

wFisher = function(p, weight = NULL, is.onetail = TRUE, eff.sign)
{
  if(is.null(weight)){weight = rep(1, length(p))}
  idx.na = which(is.na(p))
  if(length(idx.na)>0){
    p = p[-idx.na];
    weight = weight[-idx.na];
    if(!is.onetail)
    {
      eff.sign = eff.sign[-idx.na]
    }
  }
  NP = length(p)
  NS = length(weight)
  if(NP!=NS){stop("The length of p and weight vector must be identical.")}
  N = NS
  Ntotal = sum(weight)
  ratio = weight/Ntotal
  Ns = N*ratio
  G = c()

  if(is.onetail)
  {
    for(i in 1:length(p))
    {
      G = append(G, qgamma(p = p[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign > 0)
    idx_neg = which(eff.sign < 0)
    # positive direction
    G = c()
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1-p[idx_neg]/2
    for(i in 1:length(p1))
    {
      G = append(G, qgamma(p = p1[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP1 = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
    # negative direction
    G = c()
    p2[idx_pos] = 1-p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2
    for(i in 1:length(p2))
    {
      G = append(G, qgamma(p = p2[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP2 = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
    resultP = 2* min(resultP1, resultP2)
    if(resultP > 1.0){resultP = 1.0}
    overall.eff.direction = if(resultP1<=resultP2){"+"}else{"-"}
  }
  RES = if(is.onetail){list(p=min(1,resultP))}else{list(p=min(1,resultP), overall.eff.direction=overall.eff.direction)}
  return(RES)
}


#' @title wZ
#' @description P-value combination based on weighted Z-method
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weights (e.g., sample sizes)
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @return p : Combined p-value
#' @return overall.eff.direction : The direction of combined effects.
#' @return sumz : Sum of transformed z-score
#' @examples wZ(p=c(0.01,0.2,0.8), weight = c(20,10,40), is.onetail=FALSE, eff.sign=c(1,-1,1))
#' @import "metap"
#' @references Becker BJ (1994). “Combining significance levels.” In Cooper H, Hedges LV (eds.), A handbook of research synthesis, 215–230. Russell Sage, New York.
#' @references Stouffer SA, Suchman EA, DeVinney LC, Star SA, Williams RMJ (1949). The American soldier, vol 1: Adjustment during army life. Princeton University Press, Princeton.
#' @references Mosteller, F. & Bush, R.R. (1954). Selected quantitative techniques. In: Handbook of Social Psychology, Vol. 1 (G. Lindzey, ed.), pp. 289–334. Addison‐Wesley, Cambridge, Mass.
#' @export
wZ = function(p, weight = NULL, is.onetail = TRUE, eff.sign)
{
  if(length(p) == 0){stop("No input p-values exist.")}
  if(is.null(weight)){weight = rep(1, length(p))}
  idx.na = which(is.na(p))
  if(length(idx.na)>0){
    p = p[-idx.na];
    weight = weight[-idx.na]
    if(!is.onetail)
    {
      eff.sign = eff.sign[-idx.na]
    }
  }
  if(is.onetail){res = sumz(p = p, weights = weight)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign >= 0)
    idx_neg = which(eff.sign < 0)
    # positive direction
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1 - p[idx_neg]/2
    res1 = sumz(p = p1, weights = weight)
    # negative direction
    p2[idx_pos] = 1 - p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2
    res2 = sumz(p = p2, weights = weight)

    if(res1$p <= res2$p){res = res1;overall.eff.direction = "+"}else{res = res2; overall.eff.direction = "-"}
  }
  if(is.onetail)
  {
    RES = list(p = res$p[,1], sumz = res$z[,1])
  }else{
    RES = list(p = res$p[,1] * 2 , overall.eff.direction = overall.eff.direction, sumz = res$z[,1])
  }

  return(RES)
}

#' @title Lancaster
#' @description P-value combination based on Lancaster's procedure
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weights (e.g., samples sizes)
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes (1 or -1). It works when is.onetail = FALSE
#' @return p : Combined p-value
#' @return overall.eff.direction : The direction of combined effects.
#' @examples lancaster(p=c(0.01,0.2,0.8), weight=c(20,50,10), is.onetail=FALSE, eff.sign=c(1,1,1))
#' @importFrom "stats" "pchisq" "qchisq"
#' @references Becker BJ (1994). “Combining significance levels.” In Cooper H, Hedges LV (eds.), A handbook of research synthesis, 215–230. Russell Sage, New York.
#' @references Lancaster HO (1949). “Combination of probabilities arising from data in discrete distributions.” Biometrika, 36, 370–382.
#' @export

lancaster = function(p, weight, is.onetail=TRUE, eff.sign)
{
  if(is.null(p)){stop("Input p-values are required.")}
  if(is.null(weight)){stop("Weights are required.")}
  idx.na = which(is.na(p))
  if(length(idx.na)>0){
    p = p[-idx.na];
    weight = weight[-idx.na];
    if(!is.onetail)
    {
      eff.sign = eff.sign[-idx.na]
    }
  }

  NP = length(p)
  NS = length(weight)
  if(NP!=NS){stop("The length of pvalue and samplesize vector must be identical")}
  N = NP
  Ntotal = sum(weight)
  SS = weight
  #Ns = N*ratio
  if(is.onetail)
  {
    X = c()
    for(i in 1:N)
    {
      X = append(X, qchisq(p = p[i],df = SS[i], lower.tail=F))
    }
    Xsum = sum(X)
    resultP = pchisq(q=Xsum, df = Ntotal, lower.tail=F)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign > 0)
    idx_neg = which(eff.sign < 0)
    # positive direction
    X = c()
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1-p[idx_neg]/2
    for(i in 1:N)
    {
      X = append(X, qchisq(p = p1[i],df = SS[i], lower.tail=F))
    }
    Xsum = sum(X)
    resultP1 = pchisq(q=Xsum, df = Ntotal, lower.tail=F)
    # negative direction
    X = c()
    p2[idx_pos] = 1-p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2
    for(i in 1:N)
    {
      X = append(X, qchisq(p = p2[i],df = SS[i], lower.tail=F))
    }
    Xsum = sum(X)
    resultP2 = pchisq(q=Xsum, df = Ntotal, lower.tail=F)
    resultP = 2* min(resultP1, resultP2)
    if(resultP>1.0){resultP = 1.0}
    overall.eff.direction = if(resultP1<=resultP2){"+"}else{"-"}
  }
  RES = if(is.onetail){list(p=min(resultP,1))}else{list(p=min(resultP,1), overall.eff.direction=overall.eff.direction)}
  return(RES)
}
