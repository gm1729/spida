\name{adult}
\alias{adult}
\docType{data}
\title{
The Adult dataset
}
\description{
This dataset contains census data and is available from the UCI Machine Learning Repository. There are three variables with missing values.
}
\usage{data(adult)}
\format{
  A data frame with 48842 observations on the following 14 variables.
  \describe{
    \item{\code{V1}}{age: continuous variable}
    \item{\code{V2}}{workclass: nominal variable}
    \item{\code{V3}}{fnlwgt: continuous variable}
    \item{\code{V4}}{education: nominal variable}
    \item{\code{V5}}{marital-status: nominal variable }
    \item{\code{V6}}{occupation: nominal variable}
    \item{\code{V7}}{relationship: nominal variable }
    \item{\code{V8}}{race: nominal variable}
    \item{\code{V9}}{sex: nominal variable}
    \item{\code{V10}}{capital-gain: continuous variable}
    \item{\code{V11}}{capital-loss: continuous variable}
    \item{\code{V12}}{hours-per-week: continuous variable}
    \item{\code{V13}}{native-country: nominal variable}
    \item{\code{V14}}{Income class: nominal variable}
  }
}
\details{
The values of the nominal attributes are as follows: \cr

workclass: Private, Self-emp-not-inc, Self-emp-inc, Federal-gov, Local-gov, State-gov, Without-pay, Never-worked. \cr 
education: Preschool, 1st-4th, 5th-6th, 7th-8th, 9th, 10th, 11th, 12th, HS-grad, Some-college, Assoc-voc, Assoc-acdm, Bachelors, Masters, Prof-school, Doctorate \cr
marital-status: Married-civ-spouse, Divorced, Never-married, Separated, Widowed, Married-spouse-absent, Married-AF-spouse. \cr
occupation: Tech-support, Craft-repair, Other-service, Sales, Exec-managerial, Prof-specialty, Handlers -cleaners, Machine-op-inspct, Adm-clerical, Farming-fishing, Transport-moving, Priv-house-serv, Protective-serv, Armed-Forces. \cr
relationship: Wife, Own-child, Husband, Not-in-family, Other-relative, Unmarried. \cr
race: White, Asian-Pac-Islander, Amer-Indian-Eskimo, Other, Black. \cr
sex: Female, Male. \cr
native-country: United-States, Cambodia, England, Puerto-Rico, Canada, Germany, Outlying- US (Guam-USVI-etc), India, Japan, Greece, South, China, Cuba, Iran, Honduras, Philippines, Italy, Poland, Jamaica, Vietnam, Mexico, Portugal, Ireland, France, Dominican-Republic, Laos, Ecuador, Taiwan, Haiti, Columbia, Hungary, Guatemala, Nicaragua, Scotland, Thailand, Yugoslavia, El-Salvador, Trinadad&Tobago, Peru, Hong, Holand-Netherlands. \cr
Income class:  >50K, <=50K. \cr
}
\source{
The UCI Machine Learning Database Repository at: \cr
\url{http://archive.ics.uci.edu/ml/}

}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(adult)
summary(adult)
%%## maybe str(adult) ; plot(adult) ...
}
\keyword{datasets}
