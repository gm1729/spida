\name{Leff}
\alias{Leff}
\alias{fdiff}
\alias{flevs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypothesis matrices for first and higher-order effects
}
\description{
Easier visualization, estimation, testing and graphing of specific effects in models with complex interactions.
}
\usage{
Leff(fit, x = NULL, data, debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
a model with a linear formula.
}
  \item{x}{
the name of the variable(s) with respect to which specific effects are needed
}
  \item{data}{
a data frame of values where specific effects are estimated
%%     ~~Describe \code{data} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
}
\details{

Functions for effect plots
       Leff, Leff2  - generate L matrices for effect plots and second order effect plots
       Linfo     - information on variables in model
       Lallterms - show all terms that involve an effect
       Lcall     - generate code for manual editing for a cbind call to generate an L matrix
Also in this script:
       onames - pretty labels for variable in xyplot
       allp   - recursive function to generate all permutations

USER'S GUIDE:

Plotting effects: (call the variable of interest 'x')

Effects can be plotted in a number of ways. To help a reader visualize
an effect it is probably best to plot them in more than one way:

  1) response plot: yhat is plotted against 'x' controlling for other variables

     a) variables that have no interactions with 'x' can be set at a value, their
        only effect on the plot is to move the entire graph up or down but they do
        not affect the shape of the graph.

     b) Variables that do interact significantly with x should be set at different
        levels and can be visualized with panels ( y ~ x | z) in xyplot or with
        'groups'. Different choices will bring different features into prominence.

  2) effect plot: the effect of 'x' on 'y' is plotted against other variables
         that have interactions with 'x'.
         If 'x' is involved in 3-way or higher interactions, then the effect
         of 'x can be plotted against one of the interacting variables
         controlling for others using panels or groups in xyplot.

  Response plots are produced in the familiar way: create a 'pred' data
  frame and then use 'predict' to add yhat to the data frame.

  The following discusses the creation of effect plots using a number of
  tools:

       Leff   generates an L matrix to estimate the effect of 'x' as a function
              of other variables

       Linfo  extracts useful information about variables in a model and
              their interactions with 'x'

       Lallterms  generates a list of terms that can be used with 'wald'
              [essentially it replaces wald(fit , allp3( 'x','y','z')) etc.]
              Lallterms can be used in a way to generate only the overall
              F tests without showing the individual terms although it is
              desirable to look at the individual terms.]

       Lcall  generates a call to cbind with the code to generate each column
              of the model matrix. The code can be edited to, e.g., differentiate
              each term to produce and effect L matrix manually if the effect
              of a complex term is needed.

  The following describes the use of Leff to generate effect plots

  Leff generates an L matrix to estimate the effect of a variable as
  a function of the variables that it interacts with

  Limitation: Leff works only for a numerical variable or
  a factor that enters linearly, e.g. no squares or polynomials or functions
  To work with a variable that enters in a polynomial, use 'Lcall' and
  then differentiate each term manually.

  Using Leff

  1) Choose a linear variable (let's call it x) whose effect you want to explore

  2) Study how the effect of x depends on other variables, i.e. what variables
     have interactions with x in the model.
     Doing this is simplified by using:

     >  lallt <- Lallterms ( fit, 'x' )   # generates all terms including x
     >  lallt
     >  wald( fit, Lallterms( fit, 'x' ))

     To see just the summary F tests, use:

     >  summ ( wald( fit, Lallterms( fit, 'x' )) )
     or
     >  round ( summ ( wald( fit, Lallterms( fit, 'x' )) ) , 5)

     This shows whether the effect of 'x' depends on other variables
     individually or jointly and provides guidance regarding the effects
     that should be displayed. For example, if another variable 'z' only has
     weak interactions, then it may be sufficient to control for 'z' without
     displaying what happens when 'z' is varied but to simply hold z at a
     median value.

     If a variable has no interactions with 'x', then the effect of 'x' is
     invariant wrt that variable and the value at which the variable is set
     is immaterial although one must choose a a value for the function 'Leff'
     to work.

  3) Identify
           a) the basic variables in the model, and
           b) the basic variables on which the effect depends, i.e.
              the basic variables that are involved in interaction with x
     Doing this is simplified with:

     > Linfo( fit, 'x' )

  4) Generate response plots:
     a) generate a prediction data frame with all variables
     b) plot desired responses curves to show effects

  5) Create an effect data frame for effects 'pred.eff' in which all variables in the model appear.
           a) the basic numerical variables in the model that do not interact with 'x'
              only need to be present at an arbitrary level (as mentioned above).
           b) Factors should be present at all levels. e.g. sex = levels(zdc$sex)
           c) Variables that interact with x should be present in the values
              you intend to plot
           d) If x is numerical: include only x = 1
              If x is a factor, include it at all its levels: e.g. sex = levels(zdc$sex)
     The code for the prediction data frame would look like:

     pred.eff <- expand.grid( sex = levels(zdc$sex) , age = 12:80, incomez = seq(-4,4,2), x = 1)


  6) Generate the effect matrix:

     > Lf <- Leff( fit, 'x', pred.eff)    # Note 'x' should be in quotes.
     and inspect it:
     > head( Lf )
     Make sure it makes sense ( this is an experimental function )

  7) Estimated the effects using 'wald'

     > ww <- as.data.frame( wald( fit, Lf))

     Merge with pred.eff

     > wc <- cbind( ww, pred.eff)

  8) Plot the effect with xyplot.

     If 'x' is a factor, you need to select levels of the factor. The effect
     plot will show the difference between the selected level and the reference level
     To compare two levels that are NOT reference levels is IN DEVELOPMENT

  9) To easily compare the effect at different levels of an interacting variable
     is laborious and can be done manually (see example below). A more effective
     way is IN DEVELOPMENT.

  9) DONE!
     Use 'Leff2'
     See 'tutorial' for example using ch.abused and sex
     The idea is to set both variable at the value 1




}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
