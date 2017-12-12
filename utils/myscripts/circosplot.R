###########################
#### Circos Plots in R ####
###########################

# mrs72@cam.ac.uk
# ms37@sanger.ac.uk

library(circlize)

circos.SVs <- function(x, gx, title, colour){
  width <-5000
  height <-5500
  chrnum <- length(gx$chrname)

  D1color <- heat.colors(chrnum, alpha=0.8)
  D2color <- rev(heat.colors(chrnum, alpha=0.8))
  Dcolor <- append(D1color, D2color)         # major_chrs
  Scolor <- heat.colors(chrnum, alpha=0.6)   # re-arrengment/misjoints

  chromosome.ranges <- matrix(0, ncol=2, nrow=chrnum*2)
  chromosome.ranges[,1] <- 1
  v1 <-  as.vector(gx$chrlen)
  v2 <-  rev(as.vector(gx$syglen))
  
  chromosome.ranges[,2] <-cbind(v1,v2)
  rn1 <- as.vector(gx$chrname)
  rn2 <- rev(as.vector(gx$sygname))
  rownames(chromosome.ranges) <- append(rn1,rn2)
  
  if( length(unique(rn1)) != length(unique(rn2))){  #anyDuplicated(rn1) included in this case already
    #print("Warning: some chrs are mapped to multiple scaffolds..")
    print("Error: possibly unplaced scaffolds are still in the alignment files?\n ..you might need to change the \"unplaced\" scaffolds name in your  mysettings.sh scripts to match the name of you unplaced scaffolds OR remove any scaffolds that wrongly map to a chr from the alignment list and then try again")
    return(1)
    #stop("possibly unplaced scaffolds are still in the alignment files?\n ..you might need to change the \"unplaced\" scaffolds name in your  mysettings.sh scripts to match the name of you unplaced scaffolds OR remove any scaffolds that wrongly map to a chr from the alignment list and then try again")
  }

  
  chrnames <-   sub('.*\\_', '',  rn1)
  sygnames <- sub('.*\\_', '',  rn2)

  gnames <- append(chrnames,sygnames)

  # fix gaps size between bands
  basic_gap=1
  large_gap=20

  chr_len=sum(as.numeric(v1))
  scaff_len=sum(as.numeric(v2))

  #gaps <- as.vector(large_gap)
  #gaps <- append(gaps,rep(basic_gap,chrnum-1))

  gaps <- as.vector(rep(basic_gap,chrnum-1))
  gaps <- append(gaps, large_gap)  # between last chr and scaffold
  gaps <- append(gaps,rep(basic_gap,chrnum-1))
  gaps <- append(gaps, large_gap)  # between 0 and 0


  # 2. Start plotting
  png(paste0("my_plot", '.png'), width = width, height = height)

  par(mar=c(3,3,3,3))
  circos.par("track.height" = 0.15, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.1,0.1), 
             start.degree = 80, gap.degree = gaps)
  circos.initialize(factors = rownames(chromosome.ranges), 
                    xlim=chromosome.ranges)
  circos.trackPlotRegion(ylim = c(-1, 1),   #c(-1, 1), 
                         bg.col=Dcolor, bg.border='black',
                         bg.lwd=3)

  #Add chromosome names
  for (i in 1:length(rownames(chromosome.ranges))){
    circos.axis(h=0,sector.index = rownames(chromosome.ranges)[i], 
                major.at = chromosome.ranges[i,2]/2, 
                labels = gnames[i], 
                direction = "outside", labels.cex=8,
                lwd=9, labels.niceFacing=T)
  }

  
  gx$sygi = 1
  gx$chri = 1
  gx <- gx[,c("sygname","sygi", "syglen","chrname", "chri", "chrlen")] # add initial position

  circos.genomicLink(region1=gx[,1:3],region2=gx[,4:6], 
                     col=Scolor, lwd=5,
                     rou=0.75)

 
  cols <- c(2,3:4)
  circos.genomicLink(region1=x[,cols],region2=x[,5:7], 
                     col="dodgerblue1", lwd=5, rou=0.75)
   
  text(0,1.1,title, cex=15)
  text(-0.75,0.85,"Scaffolds", cex=20)
  text(0.75,0.85,"Chrs", cex=20)
  circos.clear()
  invisible(garbage <- dev.off())
  return(0)
}


circos.noscaffs <- function(x, gx, title, colour){

  width <-5000
  height <-5500
  chrnum <- length(gx$chrname)

  D1color <- heat.colors(chrnum, alpha=0.8)
  D2color <- rev(heat.colors(chrnum, alpha=0.8))
  Dcolor <- append(D1color, D2color)         # major_chrs
  Scolor <- heat.colors(chrnum, alpha=0.6)   # re-arrengment/misjoints

  chromosome.ranges <- matrix(0, ncol=2, nrow=chrnum)
  chromosome.ranges[,1] <- 1
  rn1 <- as.vector(gx$chrname)
  rn2 <- rev(as.vector(gx$sygname))
  chromosome.ranges[,2] <- as.vector(gx$chrlen)
  rownames(chromosome.ranges) <- gx$chrname

  if( length(unique(rn1)) != length(unique(rn2))){  
    stop("possibly unplaced scaffolds are still in the alignment files?\n ..you might need to change the \"unplaced\" scaffolds name in your  mysettings.sh scripts to match the name of you unplaced scaffolds OR remove any scaffolds that wrongly map to a chr from the alignment list and then try again")
  }
   
  chrnames <-   sub('.*\\_', '',  rn1)
  gnames <- chrnames

  basic_gap=1
  large_gap=5

 
  gaps <- as.vector(rep(basic_gap,chrnum-1))
  gaps <- append(gaps, large_gap)  # between last chr and scaffold
  
  # 2. Start plotting
  png(paste0("noscaff", '.png'), width = width, height = height)
  
  par(mar=c(3,3,3,3))
  circos.par("track.height" = 0.15, cell.padding = c(0, 0, 0, 0), track.margin = c(0.1,0.1), 
             start.degree = 85, gap.degree = gaps)
  circos.initialize(factors = rownames(chromosome.ranges), 
                    xlim=chromosome.ranges)
  circos.trackPlotRegion(ylim = c(-1, 1), 
                         bg.col=Dcolor, bg.border='black',
                         bg.lwd=3)
  
  # Add chromosome names
  for (i in 1:length(rownames(chromosome.ranges))){
      circos.axis(h=0,sector.index = rownames(chromosome.ranges)[i], 
                major.at = chromosome.ranges[i,2]/2, 
                labels = gnames[i], 
                direction = "outside", labels.cex=9,
                lwd=9,
                labels.niceFacing=T)
  }

  cols <- c(1,3:4)
  circos.genomicLink(region1=x[,cols],region2=x[,5:7], 
                     col="dodgerblue1", lwd=5, rou=0.75)
   
  text(0,1.1,title, cex=15)
  circos.clear()
  invisible(garbage <- dev.off())
}

# Inputs:
# x = a 6-column Matrix with structural variants: i.e.
# x = [1. Chromosome - 1. Start pos. - 1. End pos. - 2. Chromosome - 2. Start pos. - 3. End pos.]

# some settings
title = 'fAnaTes1_caus vs. Seabass'
mycolour = 'lightblue'

args <- commandArgs(TRUE)
if (length(args) > 0){
  title <- args[1]
}else{
  title <- "no title given"
}

# read mis-joints
file<-"als.file"
x<-read.table(file, header=FALSE)


# read globals
file<-"chrs.file" 
colnames=c("sygname","syglen","chrname", "chrlen")
gx<-read.table(file, col.names=colnames,header=FALSE)

exit=circos.SVs(x, gx, title, mycolour)
if(!exit){
  exit2=circos.noscaffs(x, gx, title, mycolour)
}else{
  q("no", 1, FALSE)  # does not give the error code really.... useless
}


#if(0)tryCatch(circos.SVs(x, gx, title, mycolour),
#         error= function(err) { print(paste("Error :", err))},
#         finally = function(out) {  if(0)print("all done")}  )     
#if(0)tryCatch(circos.noscaffs(x, gx, title, mycolour),
#         error= function(err) { stop(paste("Error :", err)) },
#         finally = function(out) {  if(0)print("all done")}  )
