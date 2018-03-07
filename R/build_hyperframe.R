library(spatstat)
library(parallel)
library(future) # This is to process several datasets in parallel

#This dource file needs the svg2psp package to work, to convert SVG images of root cell walls to psp spatstat objects
# The package can be installed using by exectutung the command devtools::install_github("xraynaud/svg2psp") (needs the devtools package)
library(svg2psp)

plan(sequential) # Work sequentially. Use plan(tweak(multiprocess,workers=detectCores()-1)) to analyse datasets in parallel



dir = getwd()
setwd(paste(dir,"Data",sep=""))

coord_df_to_psp = function(data,box3 = box3(),array_res = c(1024,1024,1024), clean=T) {
 # data = data[data$I>0,]
  ignore = c()
  txy= 2* diff(box3$xrange)/array_res[2]
  tz = 2*diff(box3$zrange)/array_res[3]
  counter=0
  df = data.frame(X=NULL,Y=NULL,Z=NULL,I=NULL,id = NULL, id2 = NULL)
  for (c in 1:dim(data)[1]) {
    if (!(c %in% ignore)) {
      test = which(data[c,]$X<data$X+txy & data[c,]$X>data$X-txy & data[c,]$Y<data$Y+txy & data[c,]$Y>data$Y-txy & data[c,]$Z<data$Z+tz)
      
      nkeep = test[which.max( data[test,]$I)][1] #Si jamais plusieurs points ont la même intensité, ne concserve que le premier.
      test = c(test,which(data[nkeep,]$X<data$X+txy & data[nkeep,]$X>data$X-txy & data[nkeep,]$Y<data$Y+txy & data[nkeep,]$Y>data$Y-txy & data[nkeep,]$Z<data$Z+tz))
      keep = data[nkeep,] # point que l'on garde
      keep$Z = sum(data[test,]$I*data[test,]$Z)/sum(data[test,]$I)
      keep$id = nkeep
      keep$id2 = i
      
      ignore = c(ignore,test)
      df = rbind(df,keep)
    }
  }
  
  return(pp3(df$X,df$Y,diff(box3$zrange) + -df$Z, box3))
}

files = list.files(paste(dir,"Images",sep=""))

filename = unlist(strsplit(files,".tif"))

compteur = 1

imgs=list()
for (i in 1:length(filename)) {
  
  imgs[[i]] = future({
    
    image = filename[i]
    root = substr(unlist(strsplit(filename[i]," - "))[1],1,nchar(unlist(strsplit(filename[i]," - "))[1])-4)
    date = as.numeric(substr(strsplit(filename[i],"_")[[1]][2],2,3))
    rep_root = unlist(strsplit(substr(unlist(strsplit(filename[i]," - "))[1],1,nchar(unlist(strsplit(filename[i]," - "))[1])-4),"rep"))[2]
    set = substr(unlist(strsplit(filename[i]," - "))[2],1,1)
    rep_obs = substr(unlist(strsplit(filename[i]," - "))[2],2,2)
    
    print(paste("Processing",rep_root,set,rep_obs))
    
    source(paste(filename[i],"_props.R",sep=""))
    V=box3(xrange,yrange,zrange,unitname="um")
    
    read.csv(paste(filename[i],"_coords.csv",sep=""),h=T)->data
    data = data[data$I >0,]
    coord_df_to_psp(data,box3 = V, array_res=array,clean=T)-> X3
    
    ppp(X3$data$x,X3$data$y,window = owin(V$xrange,V$yrange)) -> Xproj
    Xproj <- Xproj[!duplicated(Xproj)]
    
    if (F) {
      read.table(paste(filename[i],"_fov.txt",sep=""),h=T)->tmp
      
      mask=list()
      
      for (l in unique(tmp$id)) {
        mask[[l]] = tmp[tmp$id==l,]
      }
      W = owin(poly=mask)
      W = dilation.owin(W,0.5)
      W = simplify.owin(W,1)
      W = W[owin(xrange,yrange)]
    }
    
    W = owin(xrange,yrange)
    
    Xproj = Xproj[W] 
    
    svg2psp(paste(filename[i],"_segments.svg",sep=""),owin=W,marks=2,bezier=2,upward=T,rescale=T,clip=T)->walls
    
    return(hyperframe(image = image, date = date, root = root, rep_root=rep_root,set=set,rep_obs=rep_obs, X3 =X3,Xproj = Xproj,walls=walls))
  })
}

imgdata = do.call("rbind.hyperframe",lapply(imgs,value))
save(imgdata,file = "../imgf-cardfish.R")
