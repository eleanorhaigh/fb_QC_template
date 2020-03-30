
#checks a ferrybox data table (fb) for stuck values in a given telid
#' Title
#'
#' @param w ferrybox data table 
#' @param tid telid of instrument to check
#' @param repeats minimum number of repeats to flag
#'
#' @return
#' @export
#'
#' @examples
stuck_value<-function(w, tid, repeats){
  setorder(w,dateTime) #order data table by time
  y<<-w[telid==tid] #extract telid of interest
  z=y[,.(Mean)] #extract mean column
  z=as.numeric(unlist(z)) #remove from data table format
  rled=(rle(z)) #test for repeats
  tester<<-as.data.table(rep(rled$lengths>=repeats,times=rled$lengths)) #generate logical vector of where repeats occur
  if (any(grepl("stuck", colnames(w)))){
    colnames(tester)<-c("telidspecreps") #rename stuck column 
    y<<-cbind(y,tester) #combine stuck column and data of relevant telid
    w<-merge(w,y[,.(dateTime, telid, serial, Cruise, telidspecreps)], by= c("dateTime", "telid", "serial", "Cruise"),all.x=T) #merge to original fb data
    w<<-as.data.table(w %>% mutate(stuck = coalesce(stuck,telidspecreps))) #add stuck values for this variable to main stuck value column
  }else{
    colnames(tester)<-c("stuck") #rename stuck column
    y<<-cbind(y,tester) #combine stuck column and data of relevant telid
    w<<-merge(w,y[,.(dateTime, telid, stuck, Cruise, serial)], by= c("dateTime", "telid", "Cruise", "serial"),all.x=T) #merge to original fb data
  }
}

#########################################################################

#generates linear model equation
#' Title
#'
#' @param x x data
#' @param y y data
#'
#' @return
#' @export
#'
#' @examples
lm_eqn = function(x, y) {
  m = lm(y ~ x)
  eq <- substitute(italic(y) == a + b %.% italic(x) * "," ~ 
                     ~italic(r)^2 ~ "=" ~ r2, list(a = format(coef(m)[1], 
                                                              digits = 2), b = format(coef(m)[2], digits = 2), 
                                                   r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
}

##########################################################################

substrRight<- function(y, n){
  substr(y, nchar(y)-n+1, nchar(y))
} #function to extract cruise name from directory

###########################################################################

# function to read in qcd ctd data from folder x #and name the data after the cruise on which it was taken-then get positiion -relies on ctd folder structure remaining the same
read_ctd=function (x) {
  setwd(ctdfolders[x]) 
  load("CTDQC.rdata",.GlobalEnv)
  newname=cruises[x]
  assign(newname,session,envir=.GlobalEnv)
  ctdpositions[x]<<-list(session$positions)
  ctddata[x]<<-list(session$data)
  
  finalstationdata<<-list()
  
  #loop through each station of cruise and extract data into list finalstationdata
  for (y in 1:length(ctddata[[x]])){
    stationnames<<-names(ctddata[[x]])
    stationdata=slot(ctddata[[x]][[y]],"data")
    stationdataf=data.table(t(Reduce(rbind, stationdata))) #convert data to data table and tranpose
    colnames(stationdataf)=names(stationdata) #rename columns
    stationdataf=stationdataf[depth > 3.5 & depth < 4.5] #extract depths of interest
    stationdataf[,time:=time+as.numeric(ctdpositions[[x]]$startTime[[y]]) ] 
    stationdataf[,time:=as.POSIXct(time,origin="1970-01-01", tz="UTC")] #put back into date format
    names(stationdataf)[names(stationdataf) == "time"] <- "dateTime"  #rename time column
    stationdataf[, `:=`(dateTime, as.POSIXct(round(as.numeric(dateTime)/60) * 
                                               60, origin = "1970-01-01", tz = "UTC"))] #round to nearest minute
    finalstationdata[[y]]<<-stationdataf
    
  }
  finalctddata[[x]]<<-rbindlist(finalstationdata,fill=T)
  rm(finalstationdata)
  
}

#############################################################################
#performs same task as read.ferrybox.10min on cruises from 2019 onwards with bad timestamp
#bad timestamped data must be in own seperate folder
#corects timestamp sing file name


fix_2019_timestamps_fb<-function (folder, recursive = F, print_file = F, debug = F) {
  read_10min <- function(f, print_file = F) {
    if (grepl("10minfiles", f, ignore.case = F)) {
      filename=f
      f = paste0(folder, f)
      ln = readLines(f)
      if (print_file) {
        print(f)
      }
      
      dateLine = grep("Date[ /\t]Time", ln, perl = T)
      d = read.table(f, sep = "\t", header = F, skip = dateLine + 
                       1, fill = T)
      cruise = strsplit(ln[2], "\t")[[1]][2]
      SIC = strsplit(ln[3], "\t")[[1]][2]
      comment = strsplit(ln[6], "\t")[[1]][2]
      header1 = unlist(strsplit(ln[dateLine], "\t"))
      header1[header1 == ""] = NA
      header1 = zoo::na.locf(header1)
      header2 = unlist(strsplit(ln[dateLine + 1], "\t"))
      header = paste(header1, header2, sep = "~~")
      header = gsub("Date/Time", "DateTime", header)
      colnames(d) = gsub("[^[:alnum:]~/]", "", header[1:length(colnames(d))])
      d = data.table(d)
      
      if ("Date~~" %in% header) {
        d[, `:=`("DateTime~~", paste(`Date~~`, `Time~~`))]
      }
      
      d[, `:=`(dateTime, as.character(d$"DateTime~~"))]
      
      year=substr(filename,12, 15)
      month=substr(filename,16,17)
      day=substr(filename,18,19)
      time=substr(d$dateTime, 12,19)
      time=as.POSIXct(time,format="%H:%M:%OS")
      time=substr(d$dateTime, 12,19)
      
      dattim=list()
      for(n in 1:10){
        if (any(time<="00:01:00")){
          if ((time[n]<="00:20:00" & time[n]>="00:00:00") & (month=="08" | month=="10") & day=="31"){ #aug oct month change
            newday=1
            dattim[n]=as.POSIXct(paste0(year,"-",as.numeric(month)+1,"-",newday," ",time[n]), tz="UTC") 
          }else if ((time[n]<="00:20:00" & time[n]>="00:00:00") & (month=="09") & day=="30"){
            newday=1
            dattim[n]=as.POSIXct(paste0(year,"-",as.numeric(month)+1,"-",newday," ",time[n]), tz="UTC") 
          }else if(time[n]<="00:20:00" & time[n]>="00:00:00"){
            dattim[n]=as.POSIXct(paste0(year,"-",month,"-",as.numeric(day)+1," ",time[n]), tz="UTC")
          }else{
            dattim[n]=as.POSIXct(paste0(year,"-",month,"-",day," ",time[n]), tz="UTC")
          }
        }else{
          dattim[n]=as.POSIXct(paste0(year,"-",month,"-",day," ",time[n]), tz="UTC") 
        }
      }
      
      datim=as.POSIXct(unlist(dattim),tz="UTC", origin= "1970-01-01 00:00.00 UTC")
      
      d[,dateTime:=datim]
      
      d = suppressWarnings(melt.data.table(d, id.vars = grep("Course|Long|Lat|Satellite|Speed|Heading|dateTime", 
                                                             colnames(d))))
      d[, `:=`(c("variable", "unit", "telid", "serial", 
                 "stat"), tstrsplit(variable, "~~"))]
      
      colnames(d) = gsub("~~[[:alnum:]]*", "", colnames(d))
      d = na.omit(d)
      
      if (anyDuplicated(d) > 0) {
        warning(paste("duplicates found and removed", 
                      f))
      }
      
      d[,value:=as.numeric(value)]
      d = dcast.data.table(unique(d), ... ~ stat, value.var = "value", 
                           fun.aggregate = median)
      d[, `:=`(Quality, as.character(Quality))]
      d[, `:=`(Cruise, cruise)]
      d[, `:=`(SIC, SIC)]
      d[, `:=`(Comment, comment)]
      
      if (debug) {
        d[, `:=`(filename, f)]
      }
      return(data.frame(d))
    }
  }
  if (print_file) {
    recursive = T
    dat = lapply(list.files(folder, recursive = recursive), 
                 read_10min, print_file = T)
    dat = rbindlist(dat, fill = T)
  }
  else {
    dat = pbapply::pblapply(list.files(folder, recursive = recursive), 
                            read_10min)
    dat = rbindlist(dat, fill = T)
  }
  return(dat[order(dateTime)])
}