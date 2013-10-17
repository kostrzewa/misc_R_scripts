# this function uses the "stacked barplot" functionality of R to
# create a barplot representing nested data

# the first bar is the "total data"
# the second bar is a sub-part of the total data
# and so on, at least that's the intended effect

# the function required arguments: datafile and namefile
# the function has optional arguments: debug, normalize
# requirements: RColorBrewer (optional), reshape

# the datafile is composed as follows (the header is NOT optional!):

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

# further sub-parts can be added by using more bars
# the data is grouped by the "barnum" column and the bars
# are plotted from left to right with increasing "barnum" number

# the largest recommended number of divisions is 12 (i.e. 12 bar elements
# in total) because of the colour selection, but the code falls
# back onto "rainbow()" for more than 12 colours

# the namefile contains one line per nesting level which will serve as the
# label of the corresponding bar on the plot. For the example above, this file would look as follows:

# name
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
  
  
  norm <- sum(ordat[ordat$barnum==1,]$value)
  if(normalize) {
    ordat$value <- ordat$value/norm
    norm <- 1
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
  cols <- NULL
  # the largest set of colours returned by RColorBrewer has twelve colours
  # so we fall back to rainbow if we have more than 12 elements in total
  # we also fall back on rainbow if RColorBrewer is not available
  if(n <= 12 && require(RColorBrewer) ) {
    cols <- brewer.pal(n=n,name="Spectral")
  } else {
    cols <- rainbow(n=n)
  }

  # we want the legend to appear next to the plot, so we need to adjust the plot margins
  # and turn clipping off
  par(mar=c(8, 3, 4, 13), xpd=TRUE, family="Palatino")

  bp <- barplot(twide,names=namedat$name,xlab="",xaxt="n",col=cols)

  # add legends (one for each bar)
  # this x coordinate seems to work well for varying numbers of bars (tested 3 and 7)
  lxcoord <- max(bp)+max(bp)/(2*length(bp))
  # we start paining the legend at the very top of the page and then go down
  ypos <- 1.3*norm
  for( i in unique(ordat$barnum) ) {
    # number of elements in bar number i
    elements <- length(ordat[ordat$barnum==i,]$name)
    # we will also use this criterion to choose the colours from the cols vector using which
    # the legend writes top to bottom so we reverse the order, so that it matches
    # the way the stacks are drawn
    legend(x=lxcoord,y=ypos,legend=rev(ordat[ordat$barnum==i,]$name),bty="n",fill=rev(cols[which(ordat$barnum==i)]))
    # 0.11 seems to be a good increment, draw the next legend further down
    ypos <- ypos - norm*0.11*elements
  }

  # add slanted labels so that they fit and draw an axis
  text(cex=1, x=bp+0.1, y=0, namedat$name, xpd=TRUE, srt=45,adj=c(1.1,1.3))
  axis(1,labels=FALSE,tick=TRUE,tck=-0.01,at=bp)

}
