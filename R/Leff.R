# # MOVED From RDC package 2011-02-28
#
# ALL <- FALSE
#
# ##
# ##  RDC
# ##  Some functions for longitudinal data analysis
# ##
# #
# #  Functions for effect plots
# #       Leff, Leff2  - generate L matrices for effect plots and second order effect plots
# #       Linfo     - information on variables in model
# #       Lallterms - show all terms that involve an effect
# #       Lcall     - generate code for manual editing for a cbind call to generate an L matrix
# #  Also in this script:
# #       onames - pretty labels for variable in xyplot
# #       allp   - recursive function to generate all permutations
# #
# #  USER'S GUIDE:
# #
# #  Plotting effects: (call the variable of interest 'x')
# #
# #  Effects can be plotted in a number of ways. To help a reader visualize
# #  an effect it is probably best to plot them in more than one way:
# #
# #  1) response plot: yhat is plotted against 'x' controlling for other variables
# #
# #     a) variables that have no interactions with 'x' can be set at a value, their
# #        only effect on the plot is to move the entire graph up or down but they do
# #        not affect the shape of the graph.
# #
# #     b) Variables that do interact significantly with x should be set at different
# #        levels and can be visualized with panels ( y ~ x | z) in xyplot or with
# #        'groups'. Different choices will bring different features into prominence.
# #
# #  2) effect plot: the effect of 'x' on 'y' is plotted against other variables
# #         that have interactions with 'x'.
# #         If 'x' is involved in 3-way or higher interactions, then the effect
# #         of 'x can be plotted against one of the interacting variables
# #         controlling for others using panels or groups in xyplot.
# #
# #  Response plots are produced in the familiar way: create a 'pred' data
# #  frame and then use 'predict' to add yhat to the data frame.
# #
# #  The following discusses the creation of effect plots using a number of
# #  tools:
# #
# #       Leff   generates an L matrix to estimate the effect of 'x' as a function
# #              of other variables
# #
# #       Linfo  extracts useful information about variables in a model and
# #              their interactions with 'x'
# #
# #       Lallterms  generates a list of terms that can be used with 'wald'
# #              [essentially it replaces wald(fit , allp3( 'x','y','z')) etc.]
# #              Lallterms can be used in a way to generate only the overall
# #              F tests without showing the individual terms although it is
# #              desirable to look at the individual terms.]
# #
# #       Lcall  generates a call to cbind with the code to generate each column
# #              of the model matrix. The code can be edited to, e.g., differentiate
# #              each term to produce and effect L matrix manually if the effect
# #              of a complex term is needed.
# #
# #  The following describes the use of Leff to generate effect plots
# #
# #  Leff generates an L matrix to estimate the effect of a variable as
# #  a function of the variables that it interacts with
# #
# #  Limitation: Leff works only for a numerical variable or
# #  a factor that enters linearly, e.g. no squares or polynomials or functions
# #  To work with a variable that enters in a polynomial, use 'Lcall' and
# #  then differentiate each term manually.
# #
# #  Using Leff
# #
# #  1) Choose a linear variable (let's call it x) whose effect you want to explore
# #
# #  2) Study how the effect of x depends on other variables, i.e. what variables
# #     have interactions with x in the model.
# #     Doing this is simplified by using:
# #
# #     >  lallt <- Lallterms ( fit, 'x' )   # generates all terms including x
# #     >  lallt
# #     >  wald( fit, Lallterms( fit, 'x' ))
# #
# #     To see just the summary F tests, use:
# #
# #     >  summ ( wald( fit, Lallterms( fit, 'x' )) )
# #     or
# #     >  round ( summ ( wald( fit, Lallterms( fit, 'x' )) ) , 5)
# #
# #     This shows whether the effect of 'x' depends on other variables
# #     individually or jointly and provides guidance regarding the effects
# #     that should be displayed. For example, if another variable 'z' only has
# #     weak interactions, then it may be sufficient to control for 'z' without
# #     displaying what happens when 'z' is varied but to simply hold z at a
# #     median value.
# #
# #     If a variable has no interactions with 'x', then the effect of 'x' is
# #     invariant wrt that variable and the value at which the variable is set
# #     is immaterial although one must choose a a value for the function 'Leff'
# #     to work.
# #
# #  3) Identify
# #           a) the basic variables in the model, and
# #           b) the basic variables on which the effect depends, i.e.
# #              the basic variables that are involved in interaction with x
# #     Doing this is simplified with:
# #
# #     > Linfo( fit, 'x' )
# #
# #  4) Generate response plots:
# #     a) generate a prediction data frame with all variables
# #     b) plot desired responses curves to show effects
# #
# #  5) Create an effect data frame for effects 'pred.eff' in which all variables in the model appear.
# #           a) the basic numerical variables in the model that do not interact with 'x'
# #              only need to be present at an arbitrary level (as mentioned above).
# #           b) Factors should be present at all levels. e.g. sex = levels(zdc$sex)
# #           c) Variables that interact with x should be present in the values
# #              you intend to plot
# #           d) If x is numerical: include only x = 1
# #              If x is a factor, include it at all its levels: e.g. sex = levels(zdc$sex)
# #     The code for the prediction data frame would look like:
# #
# #     pred.eff <- expand.grid( sex = levels(zdc$sex) , age = 12:80, incomez = seq(-4,4,2), x = 1)
# #
# #
# #  6) Generate the effect matrix:
# #
# #     > Lf <- Leff( fit, 'x', pred.eff)    # Note 'x' should be in quotes.
# #     and inspect it:
# #     > head( Lf )
# #     Make sure it makes sense ( this is an experimental function )
# #
# #  7) Estimated the effects using 'wald'
# #
# #     > ww <- as.data.frame( wald( fit, Lf))
# #
# #     Merge with pred.eff
# #
# #     > wc <- cbind( ww, pred.eff)
# #
# #  8) Plot the effect with xyplot.
# #
# #     If 'x' is a factor, you need to select levels of the factor. The effect
# #     plot will show the difference between the selected level and the reference level
# #     To compare two levels that are NOT reference levels is IN DEVELOPMENT
# #
# #  9) To easily compare the effect at different levels of an interacting variable
# #     is laborious and can be done manually (see example below). A more effective
# #     way is IN DEVELOPMENT.
# #
# #  9) DONE!
# #     Use 'Leff2'
# #     See 'tutorial' for example using ch.abused and sex
# #     The idea is to set both variable at the value 1
# #
# #
# #
#
# ALL <- FALSE
#
# #' @export
# Anova2.lme <-
# function (mod, error, singular.ok = TRUE, ...)
# {
#     if (!missing(error)) {
#         sumry <- summary(error, corr = FALSE)
#         s2 <- sumry$sigma^2
#         error.df <- error$df.residual
#         error.SS <- s2 * error.df
#     }
#     SS.term <- function(term) {
#         which.term <- which(term == names)
#         subs.term <- which(assign == which.term)
#         relatives <- relatives(term, names, fac)
#         subs.relatives <- NULL
#         for (relative in relatives) subs.relatives <- c(subs.relatives,
#             which(assign == relative))
#         hyp.matrix.1 <- I.p[subs.relatives, , drop = FALSE]
#         hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop = FALSE]
#         hyp.matrix.2 <- I.p[c(subs.relatives, subs.term), , drop = FALSE]
#         hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop = FALSE]
#         hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0)
#             hyp.matrix.2
#         else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov(mod)))
#         hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term,
#             1, function(x) all(x == 0)), , drop = FALSE]
#         if (nrow(hyp.matrix.term) == 0)
#             return(c(SS = NA, df = 0))
#         lh <- linearHypothesis(mod, hyp.matrix.term, singular.ok = singular.ok,
#             ...)
#         abs(c(SS = lh$"Sum of Sq"[2], df = lh$Df[2]))
#     }
#     not.aliased <- !is.na(fixef(mod))
#     if (!singular.ok && !all(not.aliased))
#         stop("there are aliased coefficients in the model")
#     fac <- attr(mod$terms, "factors")
#     intercept <- has.intercept(mod)
#     I.p <- diag(length(fixef(mod)))
#     assign <- mod$assign
#     assign[!not.aliased] <- NA
#     names <- term.names(mod)
#     if (intercept)
#         names <- names[-1]
#     n.terms <- length(names)
#     p <- df <- f <- SS <- rep(0, n.terms + 1)
#     sumry <- summary(mod, corr = FALSE)
#     SS[n.terms + 1] <- if (missing(error))
#         sumry$sigma^2 * mod$df.residual
#     else error.SS
#     df[n.terms + 1] <- if (missing(error))
#         mod$df.residual
#     else error.df
#     p[n.terms + 1] <- f[n.terms + 1] <- NA
#     for (i in 1:n.terms) {
#         ss <- SS.term(names[i])
#         SS[i] <- ss["SS"]
#         df[i] <- ss["df"]
#         f[i] <- df[n.terms + 1] * SS[i]/(df[i] * SS[n.terms +
#             1])
#         p[i] <- pf(f[i], df[i], df[n.terms + 1], lower.tail = FALSE)
#     }
#     result <- data.frame(SS, df, f, p)
#     row.names(result) <- c(names, "Residuals")
#     names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
#     class(result) <- c("anova", "data.frame")
#     attr(result, "heading") <- c("Anova Table (Type II tests)\n",
#         paste("Response:", responseName(mod)))
#     result
# }
# environment(Anova2.lme) <- environment(Anova)
#
# if(ALL){
#     spses <- function( x ) gsp( x, c(-1,0,1), 2,1)[,-1]
#     fitmm <- lme( mathach ~ Sex*Sector*(ses+spses(ses)), hs, random = ~ 1 | school)
#     summary(fitmm)
#     wald( fitmm, 'spses')
#     wald( fitmm, Lallterms(fitmm,''))
#     Anova2.lme(fitmm,1)
# }
#
#
#
# #' Hypothesis matrices for first and higher-order effects
# #'
# #' Easier visualization, estimation, testing and graphing of specific effects
# #' in models with complex interactions.
# #'
# #'
# #' Functions for effect plots Leff, Leff2 - generate L matrices for effect
# #' plots and second order effect plots Linfo - information on variables in
# #' model Lallterms - show all terms that involve an effect Lcall - generate
# #' code for manual editing for a cbind call to generate an L matrix Also in
# #' this script: onames - pretty labels for variable in xyplot allp - recursive
# #' function to generate all permutations
# #'
# #' USER'S GUIDE:
# #'
# #' Plotting effects: (call the variable of interest 'x')
# #'
# #' Effects can be plotted in a number of ways. To help a reader visualize an
# #' effect it is probably best to plot them in more than one way:
# #'
# #' 1) response plot: yhat is plotted against 'x' controlling for other
# #' variables
# #'
# #' a) variables that have no interactions with 'x' can be set at a value, their
# #' only effect on the plot is to move the entire graph up or down but they do
# #' not affect the shape of the graph.
# #'
# #' b) Variables that do interact significantly with x should be set at
# #' different levels and can be visualized with panels ( y ~ x | z) in xyplot or
# #' with 'groups'. Different choices will bring different features into
# #' prominence.
# #'
# #' 2) effect plot: the effect of 'x' on 'y' is plotted against other variables
# #' that have interactions with 'x'.  If 'x' is involved in 3-way or higher
# #' interactions, then the effect of 'x can be plotted against one of the
# #' interacting variables controlling for others using panels or groups in
# #' xyplot.
# #'
# #' Response plots are produced in the familiar way: create a 'pred' data frame
# #' and then use 'predict' to add yhat to the data frame.
# #'
# #' The following discusses the creation of effect plots using a number of
# #' tools:
# #'
# #' Leff generates an L matrix to estimate the effect of 'x' as a function of
# #' other variables
# #'
# #' Linfo extracts useful information about variables in a model and their
# #' interactions with 'x'
# #'
# #' Lallterms generates a list of terms that can be used with 'wald'
# #' [essentially it replaces wald(fit , allp3( 'x','y','z')) etc.] Lallterms can
# #' be used in a way to generate only the overall F tests without showing the
# #' individual terms although it is desirable to look at the individual terms.]
# #'
# #' Lcall generates a call to cbind with the code to generate each column of the
# #' model matrix. The code can be edited to, e.g., differentiate each term to
# #' produce and effect L matrix manually if the effect of a complex term is
# #' needed.
# #'
# #' The following describes the use of Leff to generate effect plots
# #'
# #' Leff generates an L matrix to estimate the effect of a variable as a
# #' function of the variables that it interacts with
# #'
# #' Limitation: Leff works only for a numerical variable or a factor that enters
# #' linearly, e.g. no squares or polynomials or functions To work with a
# #' variable that enters in a polynomial, use 'Lcall' and then differentiate
# #' each term manually.
# #'
# #' Using Leff
# #'
# #' 1) Choose a linear variable (let's call it x) whose effect you want to
# #' explore
# #'
# #' 2) Study how the effect of x depends on other variables, i.e. what variables
# #' have interactions with x in the model.  Doing this is simplified by using:
# #'
# #' > lallt <- Lallterms ( fit, 'x' ) # generates all terms including x > lallt
# #' > wald( fit, Lallterms( fit, 'x' ))
# #'
# #' To see just the summary F tests, use:
# #'
# #' > summ ( wald( fit, Lallterms( fit, 'x' )) ) or > round ( summ ( wald( fit,
# #' Lallterms( fit, 'x' )) ) , 5)
# #'
# #' This shows whether the effect of 'x' depends on other variables individually
# #' or jointly and provides guidance regarding the effects that should be
# #' displayed. For example, if another variable 'z' only has weak interactions,
# #' then it may be sufficient to control for 'z' without displaying what happens
# #' when 'z' is varied but to simply hold z at a median value.
# #'
# #' If a variable has no interactions with 'x', then the effect of 'x' is
# #' invariant wrt that variable and the value at which the variable is set is
# #' immaterial although one must choose a a value for the function 'Leff' to
# #' work.
# #'
# #' 3) Identify a) the basic variables in the model, and b) the basic variables
# #' on which the effect depends, i.e.  the basic variables that are involved in
# #' interaction with x Doing this is simplified with:
# #'
# #' > Linfo( fit, 'x' )
# #'
# #' 4) Generate response plots: a) generate a prediction data frame with all
# #' variables b) plot desired responses curves to show effects
# #'
# #' 5) Create an effect data frame for effects 'pred.eff' in which all variables
# #' in the model appear.  a) the basic numerical variables in the model that do
# #' not interact with 'x' only need to be present at an arbitrary level (as
# #' mentioned above).  b) Factors should be present at all levels. e.g. sex =
# #' levels(zdc$sex) c) Variables that interact with x should be present in the
# #' values you intend to plot d) If x is numerical: include only x = 1 If x is a
# #' factor, include it at all its levels: e.g. sex = levels(zdc$sex) The code
# #' for the prediction data frame would look like:
# #'
# #' pred.eff <- expand.grid( sex = levels(zdc$sex) , age = 12:80, incomez =
# #' seq(-4,4,2), x = 1)
# #'
# #' 6) Generate the effect matrix:
# #'
# #' > Lf <- Leff( fit, 'x', pred.eff) # Note 'x' should be in quotes.  and
# #' inspect it: > head( Lf ) Make sure it makes sense ( this is an experimental
# #' function )
# #'
# #' 7) Estimated the effects using 'wald'
# #'
# #' > ww <- as.data.frame( wald( fit, Lf))
# #'
# #' Merge with pred.eff
# #'
# #' > wc <- cbind( ww, pred.eff)
# #'
# #' 8) Plot the effect with xyplot.
# #'
# #' If 'x' is a factor, you need to select levels of the factor. The effect plot
# #' will show the difference between the selected level and the reference level
# #' To compare two levels that are NOT reference levels is IN DEVELOPMENT
# #'
# #' 9) To easily compare the effect at different levels of an interacting
# #' variable is laborious and can be done manually (see example below). A more
# #' effective way is IN DEVELOPMENT.
# #'
# #' 9) DONE!  Use 'Leff2' See 'tutorial' for example using ch.abused and sex The
# #' idea is to set both variable at the value 1
# #'
# #' @aliases Leff fdiff flevs
# #' @param fit a model with a linear formula.
# #' @param x the name of the variable(s) with respect to which specific effects
# #' are needed
# #' @param data a data frame of values where specific effects are estimated %%
# #' ~~Describe \code{data} here~~
# #' @param debug %% ~~Describe \code{debug} here~~
# #' @note %% ~~further notes~~
# #' @author %% ~~who you are~~
# #' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
# #' @references %% ~put references to the literature/web site here ~
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' % Add one or more standard keywords, see file 'KEYWORDS' in the
# #' % R documentation directory.
# #'
# Leff  <- function( fit ,x = NULL, data,  debug = F){
#       if (length(x) != 1) stop( 'x must be of length 1')
#       #warning("Leff works only if x is linear or a factor")
#       mm <- model.matrix( formula(fit)[-2], data)
#       if( debug ) disp(terms(fit))
#       nams <- colnames(mm)
#       nams <- gsub( "^", ":", nams)
#       hasx <- grepl(paste(":",x,sep=''),nams)
#       if( debug ) disp(hasx)
#       mm[,!hasx] <- 0
#       attr(mm, 'coefs') <- colnames(mm)[hasx]
#       mm
# }
#
# Leff2  <- function( fit ,x = NULL, data,  debug = F){
#       if (length(x) != 2) stop( 'x must be of length 2')
#       #warning("Leff works only if x is linear or a factor")
#       mm <- model.matrix( formula(fit)[-2], data)
#       if( debug ) disp(terms(fit))
#       nams <- colnames(mm)
#       nams <- gsub( "^", ":", nams)
#       hasx1 <- grepl(paste(":",x[1],sep=''),nams)
#       hasx2 <- grepl(paste(":",x[2],sep=''),nams)
#       hasx12 <- hasx1 & hasx2
#       if( debug ) disp(hasx)
#       mm[,!hasx12] <- 0
#       attr(mm, 'coefs') <- colnames(mm)[hasx12]
#       mm
# }
#
#
# # GM: 2011-03-10: Lcall has been changed and moved to wald.R
# #     It can be invoked with 'Lform(fit)
# #
#
# #' @export
# Linfo <- function( fit , x = "", data = fit$data, debug = F){
#       #works only for lme objects at this time
#       nams <- names(fixef(fit))
#       tt <- terms(fit)
#       #vnams <- sapply(attr(tt, 'variables'),function(x) as.character(as.name(x))) [-(1:2)]
#       ff <- attr(tt,'factors')
#       vnams <- rownames(ff)[-1]
#       sel <- grep(x,colnames(ff))
#
#       vnams.dep <- rownames(ff)[ apply( cbind(0,ff[ , sel]),1,sum) >0]
#       ret <- list( vnams, vnams.dep)
#       names(ret) <- c("Variables in model:", paste("Variables in model that interact with \'",x,"\':",sep=""))
#       ret
# }
#
# if (ALL) {
#      Linfo(fit,'ses.m')
#
#      fit2 <- update(fit, .~ Sex*(ses.m+ses.d) + Sector)
#         Linfo(fit2,'ses.m')
#
#         sp <- function(x) gsp(x, c(-1,0,1), 2, 1)
#
#      fit3 <- update( fit, . ~   Sex*(ses.m+sp(ses.d)) + Sector)
#      Linfo(fit3,"ses.d")
#        terms(fit3)
#
#
#
#       summary(fit3)
#       ret <- list()
#
#       nams <- gsub( "^", ":", nams)   # delineate terms
#       nams <- gsub( "$", ":", nams)   # delineate terms
#       for ( ff in factors)   {
#           ff.string <- paste( ff, "([^:]*)" , sep = '')
#           disp( ff.string)
#           ff.rep <- paste(ff, " == \\'\\1\\'", sep = '')
#           disp(ff.rep)
#           nams <- gsub( ff.string, ff.rep, nams)
#       }
#      # for ( ii in seq_along(matrix)) {
#      #     mm.all   <- paste( "(:",names(matrix)[ii], "[^\\)]*\\))",sep='')
#      #     mm.match <- paste( "(",names(matrix)[ii], "[^\\)]*\\))",matrix[ii], sep ='')
#      #     mm.rep   <- paste( "\\1")
#      #     which.null <- grepl( mm.all, nams) mm.null  <-
#      #
#      # }
#       nams <- sub("(Intercept)", 1, nams)
#       nams <- gsub( "^:","(",nams)
#       nams <- gsub( ":$",")",nams)
#       nams <- gsub( ":", ") * (", nams)
#       nams <- paste( "with (data, cbind(", paste( nams, collapse = ",\n"), ")\n)\n", collapse = "")
#       nams
# }   # END
#
#
#
#
# if(ALL) {
#  cat(Lcall( fit, factors= c("Sex","Sector"))  )
#
#
# Lx <-  with (hs, cbind( ((1)),
# (Sex == 'Male'),
# (Sector == 'Public'),
# (ses.m),
# (ses.d),
# (Sex == 'Male') * (Sector == 'Public'),
# (Sex == 'Male') * (ses.m),
# (Sector == 'Public') * (ses.m),
# (Sex == 'Male') * (ses.d),
# (Sector == 'Public') * (ses.d),
# (ses.m) * (ses.d),
# (Sex == 'Male') * (Sector == 'Public') * (ses.m),
# (Sex == 'Male') * (Sector == 'Public') * (ses.d),
# (Sex == 'Male') * (ses.m) * (ses.d),
# (Sector == 'Public') * (ses.m) * (ses.d),
# (Sex == 'Male') * (Sector == 'Public') * (ses.m) * (ses.d) )
# )
#
# } # END
#
#
# ## Generate all terms that include a particular term to test and lapply to wald
# ## to get a test for everything
# ##
# ## Generalize allperms so it can handle a colon separated string
# ##
# ##   Recursive perm generator
# ##   Add new element to each position within previous permutations
# #  list ( c(1,2), c(2,1))
# # becomes
# #  list( c(3,1,2), c(1,3,2),c(1,2,3), c(3,2,1), c(2,3,1), c(2,1,3))
# #  etc
# #
#
# #' @export
# allp <- function( x ) {  # list of all permutations of x
#   insertall <- function( p, new ){
#      n <- length(p) + 1
#      ret <- rep(new,(length(p)+1)^2)
#      ret[- seq(1,n^2,n+1)]  <- rep(p, n)
#      split(ret,rep(1:n,each=n))
#   }
#   if (length( x ) == 1) return ( list(x))
#   prev <- Recall( x[-1])
#   new <- x[1]
#   ret <- list()
#   for ( p in prev ) {
#       ret <- c(ret,insertall(p,new))
#   }
#   ret
# }
# if(ALL)  system.time( allp(1:5))
#
#   # Generate a list for wald that contains all terms that contain a given term
#
#
#
# #' @export
# Lallterms <- function( fit, x ) {
#   # creates a list with all terms in a model that contain an effect
#   # The list can be used with 'wald' to test all effects
#   warning("If a variable occurs in a function, e.g. a spline, the test for the variable will include all relevant terms only if all functions of the variable appear with the same orders as the variable in interaction terms")
#   allt <- labels(terms(fit))
#   zterms <- grepv(x, allt)
#   znames <- gsub( ":", " * ", zterms)
#   zterms <- as.list(gsub( ":", ".*", zterms))
#   zterms <- gsub("\\(", "\\\\(", zterms)
#   zterms <- gsub("\\)", "\\\\)", zterms)
#   names( zterms ) <- znames
#   zlist <- as.list(zterms)
#   zlist
# }
#
#
#
# #' Hypothesis matrix for lmer objects
# #'
# #' %% ~~ A concise (1-5 lines) description of what the function does. ~~
# #'
# #' %% ~~ If necessary, more details than the description above ~~
# #'
# #' @param fit %% ~~Describe \code{fit} here~~
# #' @param nam %% ~~Describe \code{nam} here~~
# #' @note %% ~~further notes~~
# #' @author G. Monette
# #' @export
# Lall <-
# function (fit, nam)
# {
#     if (class(fit) != "lmer")
#         stop("only implemented for lmer")
#     v <- fit@frame[[nam]]
#     if (!is.factor(v))
#         stop("nam needs to specify the name of a factor")
#     lev0 <- levels(v)[1]
#     ret <- list()
#     namf <- nam
#     if (substring(namf, 1, 1) != "^")
#         namf <- paste("^", namf, sep = "")
#     ret[[nam]] <- Lmat(fit, namf)
#     ret[[paste(nam, "mu", sep = ".")]] <- Lmu(fit, nam)
#     ret[[paste(nam, "diff", sep = ".")]] <- Ldiff(fit, nam)
#     ret
# }
#
#
# if(ALL){
# ww <- wald( fit.pruned.abused.rsa, Lallterms( fit.pruned.abused.rsa , "ch.abused") )
# class(ww)
#     ww
#    Lallterms( fit.pruned.abused.rsa, "ch\\.abused")
#  Lallterms( fit.pruned.abused.rsa, "ch.abused")
#
# summ.wald <- function( x )  t(sapply(x , function(x) unlist( x[[1]])))
# round(summ(ww),5)
# }
#
#
# ### tests
#
# # Idea for a simple effect generator:
# if (ALL) {
# "
# NOTES:
# Goal:  Generate an L matrix to generate estimated effects for a linear variable (no function, no powers, just interactions)
# -- perhaps best to have an expression that can be used with with( pred, do.call(cbind,....))
# -- pred should perhaps be generated externally with 1's for a numverical variable of interest, or a non-reference level
# -- Note that code of the form  Sex == 'Male'  will work with a character variable if one modifies columns of the model matrix
#    however if one evaluates the model matrix, then one needs full factor variables. Perhaps 'factor' variables need to be
#    subsetted in the resulting data frame
# -- Possible approach:
# 1) Generate a model matrix over a data frame with
#    a) cross product of factor levels and
#    b) selected levels of numerical variables
#    c) differential variable set to 1 if numerical
# 2) Any column in which the diff. var. does not appear is set to 0
#
# Lcall: generates code to reproduce the model.matrix column by column
# -- specify factor names to allow expansion of conventional factor labels
#
# Linfo
#    Provide information: terms in model and which are involved with interactions
#
# "
# } #END
#
#
# if (ALL) {
# library(spida)
# hs$ses.m <- with(hs, capply( ses, school, mean, na.rm = T))
# hs$ses.d <- hs$ses - hs$ses.m
# fit <- lme( mathach ~ Sex*Sector*ses.m*ses.d, hs, random= ~1 | school, na.action = na.omit)
# summary(fit)
# wald(fit, 'ses.d')
# wald(fit, 'ses.m')
# wald(fit, 'Sector')
# wald(fit, 'Sex')
# } #END
#
# if(ALL) {
# pred  <- expand.grid( Sex = levels(hs$Sex), Sector = levels(hs$Sector), ses.d = c(-1,0,1), ses.m = c(-1,0,1))
#
#
#
# pred.ses.d <- subset( pred, ses.d == 1)
# pred.ses.d $ L <- Leff(fit ,'ses.d', pred.ses.d)
#  pred.ses.d $ L
#  class(pred.ses.d)
#  subset( pred.ses.d, Sex == "Male")$L
#
# str(terms(fit))
#
# attributes(model.matrix(fit))
#
# } #END
#
#
#
#
# #' @export
# onames <- function( x, pre = '', post = '', pre.pos = 'all', post.pos = 'all' , reverse = F) {
#    # makes prettier labels for strips
#     if ( is.null( names(x))) names(x) <- x
#     ret <- reorder(factor(names(x)), x)
#     if (reverse) ret <- reorder(factor(names(x)), -x)
#     nl <- length( levels(ret))
#     if( 'all' %in% pre.pos ) levels(ret) <- paste(pre, levels(ret), sep = '')
#     if( 'all' %in% post.pos ) levels(ret) <- paste(levels(ret),post, sep = '')
#     if( 'top' %in% pre.pos ) levels(ret)[nl] <- paste(pre, levels(ret)[nl], sep = '')
#     if( 'top' %in% post.pos ) levels(ret)[nl] <- paste(levels(ret)[nl], post, sep = '')
#     if( 'bottom' %in% pre.pos ) levels(ret)[1] <- paste(pre, levels(ret)[1], sep = '')
#     if( 'bottom' %in% post.pos ) levels(ret)[1] <- paste(levels(ret)[1], post, sep = '')
#     ret
# }
#
#
# #' @export
# fdiff <- function(f, sep = '-', all = F) {
#     # See similar function dlevels in wald.R
#       # creates a factor with all pairwise differences from an original factor or levels of a factor
#       # NO longer assumes default contrasts
#       # should work with model matrix to create all pairwise differences
#       # assumes no hyphens in original factor -> change sep
#       if( !is.factor (f))   f <- factor( f, levels = f)
#       cmat <- contrasts(f)
#       nams <- levels(f)
#       n <- length(nams)
#       if ( all ) { # all comparisons
#           plus <- rep( 1:n, n)
#           minus <- rep (1:n, each = n)
#           drop <- plus == minus
#           plus <- plus[!drop]
#           minus <- minus[!drop]
#       } else {  # no duplicates
#           zm <- matrix(1:n, nrow = n, ncol = n)
#           plus <- zm[col(zm) < row(zm)]
#           minus <- rep(1:(n - 1), (n - 1):1)
#       }
#       nrows <- length(plus)
#       Pmat <- matrix( 0, nrows, n)
#       Pmat[ cbind(1:nrows,minus) ] <- -1
#       Pmat[ cbind(1:nrows,plus) ] <- 1
#       rownames(Pmat) <- paste( nams[plus],sep, nams[minus], sep='')
#       Cmat <- Pmat %*% cmat
#       #return( Cmat)
#       fr <- factor(rownames(Pmat), levels=rownames(Pmat))
#       contrasts(fr,n-1) <- Cmat
#       fr
# }
#
# #' @export
# flevs <- function(f, sel = levels(f)) {
# # see similar function xlevels in wald.R
#         # allows selection of levels from a factor with correct contrasts
#         levs <- levels(f)
#         if ( is.numeric(sel)) sel <- levs[sel]
#         fr <- factor( sel, levels = levs)
#         ctf <- contrasts(f)
#         contrasts(fr,ncol(ctf)) <- ctf
#         fr
# }
#
# #' @export
# fmean <- function( f, type = c('II','III')) {
#     type = match.arg(type)
#     cmat <- contrasts(f)
#     comb <- switch( type,
#         II = table(f) ,
#         III = rbind( rep(1, length( levels(f))))
#     )
#     nam <-   switch( type,
#         II = "mean" ,
#         III = "center"
#     )
#     comb <- comb/sum(comb)
#     Cmat <- comb %*% cmat
#
#     fr <- factor(nam, c(nam,levels(f)[-1]))
#     contrasts(fr) <- Cmat[rep(1,length(levels(fr))),]
#     fr
# }
#
# #' @export
# fcenter <- function(f, type = "III") fmean( f, type = type)
#
#
#
#
#
#
#
#
#
#
#
#
#
# if(ALL) {
#       f   <- fitmm.mp
#       Lallterms(f)
#       zd <- expand.grid( x = 1:3, f = c('a','b','c'))
#       contrasts(zd$f)
#       model.matrix( ~ x + f, zd)
#
#       pdiff( zd$f, all = T)
#       pdiff( letters[1:4])
#       pdiff( Prestige$type)
#       pdiff( dl$EyeballMethod)
#
# # Experimentation on predict and model.matrix with incomplete factors
#
#  zd <- expand.grid( x =1:3, f = letters[1:3])
#  zd$y <- (zd$x) * 10*as.numeric(zd$f)
#  ff <- lm( y ~ x*f, zd)
#  predict(ff)
#  zp <- expand.grid( x = 1, f = c('a','b'), y = 1)
#  predict( ff, zp)        # seems to use ff$xlevels to get it right
#  model.matrix( ff, data = zp)   # doesn't work properly although this is what will
#                                 # allow Leff, and fdiff to fool it with the difference factor
#                                 # the problem is that factor need to be 'full' or have
#                                 # correct contrast matrices
#
#  model.matrix( ~x*f, zp)       # doesn't work properly  although this is what
#   Leff( ff, 'x', zp)
#
#   zd
#
#   zp <- expand.grid( x=1, f = flevs( zd$f))
#   zp
#   Leff( ff , 'x', zp)
#   wald( ff, Leff( ff, 'x', zp))
#   zp
#   zp2 <- expand.grid( x = 1:3, f = fdiff(zd$f))
#   zp2
#   Leff( ff , 'f', zp2)
#   wald( ff, Leff( ff, 'f', zp2))
# as.data.frame(  wald( ff, Leff( ff, 'f', zp2)))
#
#    rownames(zp2)
#
#   zpp <- expand.grid( x=1:3, f=flevs(zd$f))
#    zpp
#    zpp$y <- predict( ff, zpp)
#    xyplot( y ~ x ,zpp, groups = f, type = 'b')
#
# } # end of ALL
#
#
#
#
#
#
