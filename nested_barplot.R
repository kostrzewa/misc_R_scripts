# this function uses the "stacked barplot" functionality of R to
# create a barplot representing nested data
# the first bar is the "total data"
# the second bar is a sub-part of the total data
# and so on

# the function required arguments: datafile and namefile
# the function has optional arguments: debug, normalize

# the datafile is composed as follows:

# name barnum value
# "time spent in function a" 1 4e4
# "time spent in sub-part a1 inside function a" 1 2e4
# "time spent in sub-part a2 inside function a" 1 1e4
# "time spent in a1" 2 2e3
# "time spent in sub-part a11 inside function a1" 2 8e3
# "time spent in sub-part a12 inside function a1" 2 1e4

# the first three values make up the total height of the first bar 
#     (total time spent in function a)
# the last three values make up the total height of the second bar
#     (total time spent in sub-function a1)

# notice how the last three values add up to the value in the second line
# this is essential to get the correct proportions
# here, "time spent in a1" is thus not the total time spent
# in function a1 but the part of the time that is not spent
# in either a11 or a22

# further sub-parts can be added by using further parent values
# the data is grouped by the "barnum" column and the bars
# are plotted from left to right with increasing "barnum" number

# while double-nesting is possible in principle, it will not look very convincing

# the largest recommended number of divisions is 12 (i.e. 12 bar elements
# in total) because of the colour selection

# the namefile contains one line per nesting level which will serve as the
# label of the corresponding bar on the plot. For the example above, this would be:
# "a"
# "a1"

# the "normalize" argument normalizes all values by the total
# height of the first bar

library(reshape)

nested_barplot <- function(datafile,namefile,debug=T,normalize=T) {
  rawdat <- read.table(datafile,header=T)

  # order the raw data file by the bar number
  # this is important for the legend to work correctly
  ordat <- rawdat[sort(rawdat$barnum,index.return=T)$ix,]

  rm(rawdat)

  if(normalize) {
    norm <- sum(ordat[ordat$barnum==1,]$value)
    ordat$value <- ordat$value/norm
  }

  namedat <- read.table(namefile,header=T)

  if(debug) { 
    print(namedat)
    print(ordat) 
  }

  # the reshaping gives the correct format to the data
  wide <- reshape(data=ordat,timevar="name",idvar="barnum",direction="wide")

  # the reshaping will generate NAs because all the names are unique
  # by replacing these NAs with zeroes we make the bars stackable 
  wide[ is.na(wide) ] <- 0

  # remove the "barnum" row and transpose the data

  wide$barnum <- NULL
  twide <- t(wide)

  if(debug) { 
    print(wide) 
    print(twide)
    }

  # generate some colours for the plot using RColorBrewer
  n <- length(ordat$barnum)
  colours <- NULL
  if(n <= 12) {
    library(RColorBrewer)
    colours <- brewer.pal(n=n,name="Spectral")
  } else {
    rm(colours)
  }

  if(exists("colours")){
    barplot(twide,names=namedat$name,legend=ordat$name),col=colours)
  } else {
    barplot(twide,names=namedat$name,legend=ordat$name))
  }

}
