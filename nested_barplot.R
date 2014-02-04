# this function uses the "stacked barplot" functionality of R to
# create a barplot representing nested data

# the first bar is the "total data"
# the second bar is a sub-part of the total data
# and so on, at least that's the intended effect

# the function's required arguments: datafile and namefile
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

# the first three lines make up the total height of the first bar 
#     (total time spent in function a)
# the last three lines make up the total height of the second bar
#     (total time spent in sub-function a1)

# notice how the last three values add up to the value in the second line
# this is essential to get the correct proportions
# here, "time spent in a1" is thus not the total time spent
# in function a1 but the part of the time that is not spent
# in either a11 or a22

# further sub-parts can be added by using more bars
# the data is grouped by the "barnum" column and the bars
# are plotted from left to right with increasing "barnum" number

# the namefile contains one line per nesting level which will serve as the
# label of the corresponding bar on the plot. For the example above, this file would look as follows:

# name
# "a"
# "a1"

# the "normalize" argument normalizes all values by the total
# height of the first bar

library(reshape)

nested_barplot <- function(datafile,namefile,title="",debug=T,normalize=T,legend=T) {
  rawdat <- read.table(datafile,header=T)

  # order the raw data file by the bar number
  # this is important for the legend to work correctly
  ordat <- rawdat[sort(rawdat$barnum,index.return=T)$ix,]
  namedat <- read.table(namefile,header=T)

  # if there were offsets defined in the namedat file, parse them here
  # otherwise just generate a vector of zeroes
  offsets <- c(0)
  if( with(namedat, exists("offset")) ) { 
    for( i in 1:length(namedat$offset) ) {
      offsets[i] <- eval(parse(text=namedat$offset[i]))
    }
  } else {
    offsets <- rep(x=0,times=length(namedat$name)) 
  }
  
  rm(rawdat)
  
  norm <- sum(ordat[ordat$barnum==1,]$value)
  if(normalize) {
    ordat$value <- ordat$value/norm
    offsets <- offsets/norm
    norm <- 1
  }


  if(debug) { 
    print(namedat)
    print(ordat)
    print(offsets) 
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

  # generate 9 colours for the plot using RColorBrewer
  # fall back on rainbow if RColorBrewer is not available
  if( require(RColorBrewer) ) {
    cols <- brewer.pal(n=9,name="Set1")
  } else {
    cols <- rainbow(n=9)
  }

  # we want the legend to appear next to the plot, so we need to adjust the plot margins
  # and turn clipping off
  rightmar <- 0
  if(legend) {
    rightmar <- 13
  }
  
  par(mar=c(12, 4, 2, rightmar), xpd=TRUE, family="Palatino")

  bp <- barplot(twide,ylim=c(0,norm),names=namedat$name,xlab="",xaxt="n",col=cols,border=NA,main=title,offset=offsets,las=1)

  if(debug) {
    print("barplot bar midpoints")
    print(bp)
  }
  
  # add slanted labels so that they fit and draw an axis
  text(cex=1, x=bp+0.1, y=0, namedat$name, xpd=TRUE, srt=45,adj=c(1.1,1.3))
  axis(1,labels=FALSE,tick=TRUE,tck=-0.03,at=bp)

  # switch to relative coordinates so we can draw the legend correctly on any size plot area
  par(usr=c(0,1,0,1))

  if(legend) {
    # add legends (one for each bar)
    lxcoord <- 1.1
    # we start paining the legend at the very top of the page and then go down
    ypos <- 1.03
    increment_per_element <- 0.007
    increment_per_bunch <- 0.15
    for( i in unique(ordat$barnum) ) {
      # number of elements in bar number i
      elements <- length(ordat[ordat$barnum==i,]$name)
      # we will also use this criterion to choose the colours from the cols vector using which
      # the legend writes top to bottom so we reverse the order, so that it matches
      # the way the stacks are drawn
      # the barplot function wraps around the colour vector, so we do the same here
      legend(x=lxcoord,y=ypos,legend=rev(ordat[ordat$barnum==i,]$name),bty="n",fill=rev(cols[(which(ordat$barnum==i)-1)%%9+1]))
      # 0.069 seems to be a good height for each line, draw the next legend further down by 0.005
      ypos <- ypos - increment_per_element*elements - increment_per_bunch
    }
  }
}
