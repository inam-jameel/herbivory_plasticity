setwd("~/Desktop/Anderson_data/herbivory/data/extracting_data")
getwd()

library(writexl)

##Data
GH<-read.delim("2genGH_012822.txt", header=T)
sapply(GH,class)

head(GH,20)
tail(GH)

GH$Master_notes<-as.character(GH$Master_notes)

reproductive<-subset(GH, phenology=="reproductive")
head(reproductive,20)




##Empty data set that we will populate with summary data
Empty <-read.csv("2genGH_empty.csv", header=T)

#Empty_input<-read.csv("NPEmpty.csv", header=T)
#Empty <- Empty_input[order(Empty_input$PlantID),]
#row.names(Empty) <- 1:nrow(Empty)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(GH $exp_ID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(GH $exp_ID)[i]
  sub<-subset(GH, exp_ID==target_ID)
  
  #        ##If the subset has less than 14 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==14)==F)
  { print(sub$exp_ID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(GH$exp_ID))){ 
  
  ##Subset the data for plantID i
  target_ID<-unique(GH$exp_ID)[i]
sub<-subset(GH, exp_ID==target_ID)


##Identify the row number in Empty that corresponds to this ID number 
target_row<-as.numeric(rownames(Empty[which(Empty$exp_ID==target_ID),]))


#exp_survival - Census 14 was the last full census, alive=1 dead=0
ifelse((sub$status[14]=="alive"), Empty$survival[target_row]<-1, Empty$survival[target_row]<-0)


##I measured leaf damage at censuses 1,3,5(only counted leaves for control),13	

Empty$num_total1[target_row]<-sub$num_total[1] 
Empty$num_damaged1[target_row]<-sub$num_damaged[1] 
Empty$total_damage1[target_row]<-sub$total_damage[1]
Empty$avg_damage1[target_row]<-sub$avg_damage[1]
Empty$LAR1[target_row]<-sub$LAR[1] 

Empty$num_total3[target_row]<-sub$num_total[3] 
Empty$num_damaged3[target_row]<-sub$num_damaged[3] 
Empty$total_damage3[target_row]<-sub$total_damage[3]
Empty$avg_damage3[target_row]<-sub$avg_damage[3]
Empty$LAR3[target_row]<-sub$LAR[3] 

Empty$num_total5[target_row]<-sub$num_total[5] 
Empty$num_damaged5[target_row]<-sub$num_damaged[5] 
Empty$total_damage5[target_row]<-sub$total_damage[5]
Empty$avg_damage5[target_row]<-sub$avg_damage[5]
Empty$LAR5[target_row]<-sub$LAR[5] 

Empty$num_total13[target_row]<-sub$num_total[13] 
Empty$num_damaged13[target_row]<-sub$num_damaged[13] 
Empty$total_damage13[target_row]<-sub$total_damage[13]
Empty$avg_damage13[target_row]<-sub$avg_damage[13]
Empty$LAR13[target_row]<-sub$LAR[13] 


## Master notes
Empty$Master_notes[target_row]<-sub$Master_notes[13] 

#Bolted and Flowered, yes=1 no=0
#to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
ifelse((sum(sub$status=="reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$bolted[target_row]<-1, Empty$bolted[target_row]<-0)

#First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$flowered[target_row]<-1, Empty$flowered[target_row]<-0)


#Mature_silique_number
Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)


#Failed_silique_number 
Empty$Failed_silique_number[target_row]<-sum(sub$Failed_silique_number_removed, na.rm=T)

#Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
#if (sum(sub$Silique_number, na.rm=T)>=1) { 
  
#  Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
  
  ##If the number of failed siliques is negative, set to zero
#  if(Empty$Failed_silique_number[target_row]<0) {
    
#    Empty$Failed_silique_number[target_row]<-0
    
  } #closes inner if statement



GH100 <-subset(Empty, exp_ID=="100")
GH711 <-subset(Empty, exp_ID=="711")

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "GHSummary2021.xlsx", sheetName="Summary2021", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 
