if (FALSE){

  # explore V matrix
  library(spidadev)
  fit <- lme( mathach ~ ses*Sector , hs, random = ~ 1 + ses|school)
summary(fit)  
  str(fit$modelStruct)
  methods(class='lme')
  VarCorr(fit)
  getGroupsFormula(fit)
  getVarCov(fit)   # returns G matrix
  vcov(fit)
  head(getData(fit))
  nobs(fit)
  getGroups(fit)
  formula(fit)
  terms(fit)
  matchCoefs("ses",fit)
  diags(fit)
  
  methods(class='modelStruct')
  formula(fit$modelStruct)   ######  formula for G: can form Z
  
  str(fit$modelStruct)
  
  methods(class="reStruct")
  
  coef(fit$modelStruct$reStruct)
  formula(fit$modelStruct$reStruct)
  
  ?coef.reStruct
  
  rs1 <- reStruct(list(A = pdSymm(diag(1:3), form = ~Score),
                       B = pdDiag(2 * diag(4), form = ~Educ)))
  coef(rs1)
  
  str(summary(fit))
  corMatrix( fit$modelStruct)
  corMatrix( fit$modelStruct$reStruct)
  methods('corMatrix')
  getVarCov(fit)
  ?corMatrix
  
fit$modelStruct$varStruct
  
  
  
  
  
  }