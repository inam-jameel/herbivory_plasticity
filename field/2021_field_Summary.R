
setwd("~/Desktop/Anderson_data/herbivory/data/field/2022_herbivory/census_files")
getwd()

library(writexl)

#### summarizing estess 2021 field cohort in 2021 field season ####

#### Gothic 2021 cohort ####
g2021 <- read.csv("2021_gInam_03oct22.csv")

sapply(g2021,class)

head(g2021,20)
tail(g2021)

g2021$Master_notes<-as.character(g2021$Master_notes)

reproductive<-subset(g2021, Reproductive=="1")
head(reproductive,20)




##Empty data set that we will populate with summary data
Empty <-read.csv("Field2022_gothic21_empty.csv", header=T)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)
Empty $survival<-as.numeric(Empty $survival)
Empty $Bolted<-as.numeric(Empty $Bolted)
Empty $Flowered<-as.numeric(Empty $Flowered)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(g2021$PlantID))) {

#        ##Subset for each plant
target_ID<-unique(g2021 $PlantID)[i]
sub<-subset(g2021, PlantID==target_ID)

#        ##If the subset has less than 10 rows (number of censuses), tell me the plant ID number
if((nrow(sub)==10)==F)
{ print(sub$PlantID[1]) }

}



##The loop begins - 
for (i in 1:length(unique(g2021$PlantID))) {
	
##Subset the data for plantID i
	target_ID<-unique(g2021$PlantID)[i]
	sub<-subset(g2021, PlantID==target_ID)
	
	
##Identify the row number in Empty that corresponds to this ID number 
	target_row<-as.numeric(rownames(Empty[which(Empty$PlantID==target_ID),]))
	
	
#exp_survival - Census 10 was the last full census, alive=1 dead=0
	ifelse((sub$Status[10]=="alive"), Empty$survival[target_row]<-1, Empty$survival[target_row]<-0)
		
## SKIP ####
##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census or if they had leaves counted during the damage census (note that we can only ask if there are >= 1 leaf if the element is not NA)
#	if(Empty$Season_survival[target_row]==1 | (Empty$Leaf_Number[target_row]>=1 & is.na(Empty$Leaf_Number[target_row])==F)) {
#		Empty$Overwinter_survival[target_row]<-1
#	}

##I measured leaf damage most censuses, so we will put them all here	

Empty$LAR_1[target_row]<-sub$LAR[1] 
Empty$LAR_2[target_row]<-sub$LAR[2] 
Empty$LAR_3[target_row]<-sub$LAR[3] 
Empty$LAR_4[target_row]<-sub$LAR[4] 
Empty$LAR_5[target_row]<-sub$LAR[5] 
Empty$LAR_6[target_row]<-sub$LAR[6] 
Empty$LAR_7[target_row]<-sub$LAR[7] 
Empty$LAR_8[target_row]<-sub$LAR[8] 
Empty$LAR_9[target_row]<-sub$LAR[9] 
Empty$LAR_10[target_row]<-sub$LAR[10] 




## Master notes
Empty$Master_notes[target_row]<-sub$Master_notes[10] 

## vegetative cover
Empty$vegetative_cover[target_row]<-sub$Vegetative_cover[9] 


#Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	

	#if(Empty$survival[target_row]==0) { #Outer if statement
	
	#Subset data where Status==Dead	
	#First_death<-subset(sub, Status=="dead,check again")
	
	#Extract row name of last census where status==Alive
	#Last_alive<-rownames(tail(sub[which(sub$Status=="alive"),], n=1))
	
	##Adjust First_death data set to include only data 	after the last Alive status was observed
	#if(length(Last_alive)>=1) { #Inner if statement
	
	##Adjust subset of First_death
	#First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
		
	#} #Closes inner if statement
	
	##Assign the date of first death
#	Empty$Date_first_death_observed[target_row]<-First_death$Date[1] 
	
#} #Closes outer if statement

	
###SKIP###	
##Gopher damage or Censor - was coded as GopherDamage or Censor in the Exclude column - sometimes a plant marked as being exposed to gopher damage doesn't die - for that reason, this has to meet conditions that it was censused as gopher damaged at least once and that we didn't observe it as alive in the last census
#	if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive")) {
#	Empty$Exclude[target_row]<-"GopherDamage"
#	}

#	if(any(sub$Exclude=="Censor") & (sub$Status[8]!="Alive")) {
#	Empty$Exclude[target_row]<-"Censor"
#	}
	
#	if(is.na(Empty$Exclude[target_row])==T) {
#	Empty$Exclude[target_row]<-"Include"
#	}

	
##Date_first_gopher, again only noted if the plant was not marked as alive in the last census
##This code is a little crazy because we are also calculating the date a plant was last alive prior to gopher damage and whether the living plant died prior to gopher damage
#if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive" )) {
		
#	Empty$Date_first_gopher[target_row]<-sub$Ordinal_date[which(sub$Exclude=="GopherDamage")][1] 
	
	##Last date of being alive prior to gopher damage
#	Last_alive<-tail(sub$Ordinal_date[which(sub$Status=="Alive" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)
		
	#Last date of being dead prior to gopher damage
#	Last_dead<-tail(sub$Ordinal_date[which(sub$Status=="Dead" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)	
	
	#If the plant has both a Last_alive and Last_dead date
#	if(length(Last_alive)>=1 & length(Last_dead)>=1) {
	
	#If the plant was alive before it was dead and if it died before gopher was first observed - have to separate this if statement from the last because you can't ask if something is less than something else if never occured
#	if(Last_alive<Last_dead & Last_dead<Empty$Date_first_gopher[target_row]) {
		
		#Include plant
#		Empty$Exclude[target_row]<-"Include"
			
#	} } #closes inner if statements
	
	##Lastly, if plant died prior to gopher damage, list it as a natural death
#	if(Empty$Date_first_death_observed[target_row]<Empty$Date_first_gopher[target_row]) {
		
#		Empty$Exclude[target_row]<-"Include"
				
#	} #Closes inner if statement
	
#} #Closes outer if statement

	
	

## issues here####	
#Bolted and Flowered, yes=1 no=0
	#to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times

	ifelse((sum(sub$Reproductive=="1", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Bolted[target_row]<-1, Empty$Bolted[target_row]<-0)
	
	#First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
	ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered[target_row]<-1, Empty$Flowered[target_row]<-0)


#Subset reproductive data #also doesnot work
Repro<-subset(sub, Reproductive=="1" & (Flower_number>=1 | Silique_number>=1))


#Date_flowering (really this is the date that a plant is first observed to be reproductive)
	if(nrow(Repro)>=1) {
	Empty$Date_flowering[target_row]<-Repro$Ordinal_date[1] }
	

#Flower_color - note the nested ifelse statements - this is needed to assign a 1 for plants with pigmented flowers, a 0 to plants with white flowers, and an NA to plants that never flowered (or for which we never observed flowers) 
#	if(any(Repro$Flower_color=="Pigmented")) {
#	Empty$Flower_color[target_row]<-1 }
#	
#	if(all(Repro$Flower_color!="Pigmented") & sum(Repro$Flower_number, na.rm=T)>=1) {
#	Empty$Flower_color[target_row]<-0 }


#Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
#	if(nrow(Repro)>=1) {
#	Empty$Stem_number_flowering[target_row]<-Repro$Stem_number[1] }

#	if(nrow(Repro)>=1) {
#	Empty$Flower_number_flowering[target_row]<-Repro$Flower_number[1] }

#	if(nrow(Repro)>=1) {
#	Empty$Silique_number_flowering[target_row]<-Repro$Silique_number[1] }
	
#	if(nrow(Repro)>=1) {
#	Empty$Silique_length_flowering[target_row]<-Repro$Silique_length[1] }
	

#Height1 through Height3 at first flower
#	if(nrow(Repro)>=1) {
#	Empty$Height1_flowering[target_row]<-Repro$Height1[1] }
			
#	if(nrow(Repro)>=1) {
#	Empty$Height2_flowering[target_row]<-Repro$Height2[1] }

#	if(nrow(Repro)>=1) {
#	Empty$Height3_flowering[target_row]<-Repro$Height3[1] }
	
#		if(nrow(Repro)>=1) {
#	Empty$Height4_flowering[target_row]<-Repro$Height4[1] }
			
#	if(nrow(Repro)>=1) {
#	Empty$Height5_flowering[target_row]<-Repro$Height5[1] }

#	if(nrow(Repro)>=1) {
#	Empty$Height6_flowering[target_row]<-Repro$Height6[1] }


########Note - this if statement contains summary commands for all peak-flowering related data

#if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below

#Date_peak_flowering		
#	max_fl<-max(Repro$Flower_number, na.rm=T)
	
#	Empty$Date_peak_flowering[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
			
			
#Flower_number_peak and Silique_number_peak and Stem_number_peak
#	Empty$Flower_number_peak[target_row]<-max_fl
	
#	Empty$Silique_number_peak[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
	
#	Empty$Stem_number_peak[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]


#Heiestesst1_peak throuestess Heiestesst5_peak
#	Empty$Heiestesst1_peak[target_row]<-Repro$Heiestesst1[which(Repro$Flower_number==max_fl)][1]
		
#	Empty$Heiestesst2_peak[target_row]<-Repro$Heiestesst2[which(Repro$Flower_number==max_fl)][1]

#	Empty$Heiestesst3_peak[target_row]<-Repro$Heiestesst3[which(Repro$Flower_number==max_fl)][1]

#	Empty$Heiestesst4_peak[target_row]<-Repro$Heiestesst4[which(Repro$Flower_number==max_fl)][1]
		
#	Empty$Heiestesst5_peak[target_row]<-Repro$Heiestesst5[which(Repro$Flower_number==max_fl)][1]

#	Empty$Heiestesst6_peak[target_row]<-Repro$Heiestesst6[which(Repro$Flower_number==max_fl)][1]
	

		
#} #Closes if statement




########Data for date of first SO - Note that this is the date of first SO observed after the last flower was observed

#if(sum(Repro$Silique_number, na.rm=T)>=1) { #outer if statement

	#Extract the rowname for the last day a flower was observed, and then the rownumber of that rowname within sub (note that the rownames in sub will come from its name within the master, long-format data file)
#	last_Fl<-rownames(tail(Repro[which(Repro$Flower_number>=1),], n=1))
	
	##If there is a date of last flower, adjust the repro subset to reflect only data after the date of last flower
#	if(length(last_Fl)>0) { #inner if statement

	##Row number of date of last flower in sub
#	last_Fl<-which(grepl(last_Fl, rownames(Repro)))

	##Subset data to exclude censuses on or before the last flower was observed
	
#	Repro<-Repro[-c(1:last_Fl),]
		
#} #Closes inner if statement


#Subset for censuses where a silique was observed (we have to do this second subset to exlude plants that became b or b, tco after flowering)
#	Repro<-subset(Repro, Silique_number>=1)

	
##Now estimate the date of first silique
#	Empty$Date_silique[target_row]<-Repro$Ordinal_date[1]


#Stem_number_silique, Silique_number_silique, and Silique_length_silique
#	Empty$Stem_number_silique[target_row]<-Repro$Stem_number[1] 
	
#	Empty$Silique_number_silique[target_row]<-Repro$Silique_number[1]
	
#	Empty$Silique_length_silique[target_row]<-Repro$Silique_length[1]


#Height1_silique
#	Empty$Height1_silique[target_row]<-Repro$Height1[1]
	
#	Empty$Height2_silique[target_row]<-Repro$Height2[1]

#	Empty$Height3_silique[target_row]<-Repro$Height3[1]

#	Empty$Height4_silique[target_row]<-Repro$Height4[1]
	
#	Empty$Height5_silique[target_row]<-Repro$Height5[1]

#	Empty$Height6_silique[target_row]<-Repro$Height6[1]

	
 #Closing bracket from outer if statement



#Mature_silique_number
	Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)


#Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
	if (sum(sub$Silique_number, na.rm=T)>=1) { 
		
	Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)

	##If the number of failed siliques is negative, set to zero
	if(Empty$Failed_silique_number[target_row]<0) {
		
		Empty$Failed_silique_number[target_row]<-0
		
	} #closes inner if statement
}} #closes outer if statement


#Mature_length_siliques - note that the column numbers have to be updated
#	Empty$Mature_length_siliques[target_row]<-sum(sub[,c(34,36,38,40,42,44,46,48,50,52,54,56,58,60,62)], na.rm=T)


#Total_flower_number
#	Empty$Total_flower_number[target_row]<-sum(sub$Flower_number, na.rm=T)
	
#Mature_entire_silique_length (when %PD is measured) - first vector has column numbers of silique lengths, second vector has column numbers of PD% - columns must be updated before calculating
#	Empty$Mature_PD[target_row]<-(sum(sub[,c(35,37,39,41,43,45,47,49,51,53,55,57,59,61,63)], na.rm=T))/sum(sub$Number_collected_siliques, na.rm=T)
	
	
#} ##Closes entire loop


G15001 <-subset(Empty, PlantID=="G15001")
G15001
G15246 <-subset(Empty, PlantID=="G15246")
G15246
library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "2022_gothic21_summary.xlsx", sheetName="Summary2022", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 









#### gothic 2021 ####

##Data
gothic <- read.csv("Field2021_gothic.csv")

sapply(gothic,class)

head(gothic,20)
tail(gothic)

gothic$Master_notes<-as.character(gothic$Master_notes)
gothic$Number_collected_siliques<-as.numeric(gothic$Number_collected_siliques)
gothic$Mature_length<-as.numeric(gothic$Number_collected_siliques)

reproductive<-subset(gothic, Phenology=="reproductive")
head(reproductive,20)




##Empty data set that we will populate with summary data
Empty <-read.csv("Field2021_gothic2021_empty.csv", header=T)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)
Empty $survival<-as.numeric(Empty $survival)
Empty $bolted<-as.numeric(Empty $bolted)
Empty $flowered<-as.numeric(Empty $flowered)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(gothic$PlantID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(gothic $PlantID)[i]
  sub<-subset(gothic, PlantID==target_ID)
  
  #        ##If the subset has less than 18 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==18)==F)
  { print(sub$PlantID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(gothic$PlantID))) {
  
  ##Subset the data for plantID i
  target_ID<-unique(gothic$PlantID)[i]
  sub<-subset(gothic, PlantID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$PlantID==target_ID),]))
  
  
  #exp_survival - Census 12 was the last full census, alive=1 dead=0
  ifelse((sub$Status[12]=="alive"), Empty$survival[target_row]<-1, Empty$survival[target_row]<-0)
  
  ## SKIP ####
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census or if they had leaves counted during the damage census (note that we can only ask if there are >= 1 leaf if the element is not NA)
  #	if(Empty$Season_survival[target_row]==1 | (Empty$Leaf_Number[target_row]>=1 & is.na(Empty$Leaf_Number[target_row])==F)) {
  #		Empty$Overwinter_survival[target_row]<-1
  #	}
  
  ##I measured leaf damage most censuses, so we will put them all here	
  
  Empty$LAR_1[target_row]<-sub$LAR[1] 
  Empty$LAR_2[target_row]<-sub$LAR[2] 
  Empty$LAR_3[target_row]<-sub$LAR[3] 
  Empty$LAR_4[target_row]<-sub$LAR[4] 
  Empty$LAR_5[target_row]<-sub$LAR[5] 
  Empty$LAR_6[target_row]<-sub$LAR[6] 
  Empty$LAR_7[target_row]<-sub$LAR[7] 
  Empty$LAR_8[target_row]<-sub$LAR[8] 
  Empty$LAR_9[target_row]<-sub$LAR[9] 
  Empty$LAR_10[target_row]<-sub$LAR[10] 
  Empty$LAR_11[target_row]<-sub$LAR[11] 
  Empty$LAR_12[target_row]<-sub$LAR[12]
  
  
  
  ## Master notes
  Empty$Master_notes[target_row]<-sub$Master_notes[18] 
  
  ## vegetative cover measured on census 7
  Empty$vegetative_cover[target_row]<-sub$X..Vegetative.cover[7] 
  
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  #if(Empty$survival[target_row]==0) { #Outer if statement
  
  #Subset data where Status==Dead	
  #First_death<-subset(sub, Status=="dead,check again")
  
  #Extract row name of last census where status==Alive
  #Last_alive<-rownames(tail(sub[which(sub$Status=="alive"),], n=1))
  
  ##Adjust First_death data set to include only data 	after the last Alive status was observed
  #if(length(Last_alive)>=1) { #Inner if statement
  
  ##Adjust subset of First_death
  #First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
  
  #} #Closes inner if statement
  
  ##Assign the date of first death
  #	Empty$Date_first_death_observed[target_row]<-First_death$Date[1] 
  
  #} #Closes outer if statement
  
  
  ###SKIP###	
  ##Gopher damage or Censor - was coded as GopherDamage or Censor in the Exclude column - sometimes a plant marked as being exposed to gopher damage doesn't die - for that reason, this has to meet conditions that it was censused as gopher damaged at least once and that we didn't observe it as alive in the last census
  #	if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive")) {
  #	Empty$Exclude[target_row]<-"GopherDamage"
  #	}
  
  #	if(any(sub$Exclude=="Censor") & (sub$Status[8]!="Alive")) {
  #	Empty$Exclude[target_row]<-"Censor"
  #	}
  
  #	if(is.na(Empty$Exclude[target_row])==T) {
  #	Empty$Exclude[target_row]<-"Include"
  #	}
  
  
  ##Date_first_gopher, again only noted if the plant was not marked as alive in the last census
  ##This code is a little crazy because we are also calculating the date a plant was last alive prior to gopher damage and whether the living plant died prior to gopher damage
  #if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive" )) {
  
  #	Empty$Date_first_gopher[target_row]<-sub$Ordinal_date[which(sub$Exclude=="GopherDamage")][1] 
  
  ##Last date of being alive prior to gopher damage
  #	Last_alive<-tail(sub$Ordinal_date[which(sub$Status=="Alive" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)
  
  #Last date of being dead prior to gopher damage
  #	Last_dead<-tail(sub$Ordinal_date[which(sub$Status=="Dead" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)	
  
  #If the plant has both a Last_alive and Last_dead date
  #	if(length(Last_alive)>=1 & length(Last_dead)>=1) {
  
  #If the plant was alive before it was dead and if it died before gopher was first observed - have to separate this if statement from the last because you can't ask if something is less than something else if never occured
  #	if(Last_alive<Last_dead & Last_dead<Empty$Date_first_gopher[target_row]) {
  
  #Include plant
  #		Empty$Exclude[target_row]<-"Include"
  
  #	} } #closes inner if statements
  
  ##Lastly, if plant died prior to gopher damage, list it as a natural death
  #	if(Empty$Date_first_death_observed[target_row]<Empty$Date_first_gopher[target_row]) {
  
  #		Empty$Exclude[target_row]<-"Include"
  
  #	} #Closes inner if statement
  
  #} #Closes outer if statement
  
  
  
  
  ## issues here####	
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  
  ifelse((sum(sub$Phenology=="reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$bolted[target_row]<-1, Empty$bolted[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered[target_row]<-1, Empty$Flowered[target_row]<-0)
  
  
  #Subset reproductive data #also doesnot work
  Repro<-subset(sub, Phenology=="reproductive" & (Flower_number>=1 | Silique_number>=1))
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering[target_row]<-Repro$Ordinal_date[1] }
  
  
  #Flower_color - note the nested ifelse statements - this is needed to assign a 1 for plants with pigmented flowers, a 0 to plants with white flowers, and an NA to plants that never flowered (or for which we never observed flowers) 
  #	if(any(Repro$Flower_color=="Pigmented")) {
  #	Empty$Flower_color[target_row]<-1 }
  #	
  #	if(all(Repro$Flower_color!="Pigmented") & sum(Repro$Flower_number, na.rm=T)>=1) {
  #	Empty$Flower_color[target_row]<-0 }
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  #	if(nrow(Repro)>=1) {
  #	Empty$Stem_number_flowering[target_row]<-Repro$Stem_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Flower_number_flowering[target_row]<-Repro$Flower_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Silique_number_flowering[target_row]<-Repro$Silique_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Silique_length_flowering[target_row]<-Repro$Silique_length[1] }
  
  
  #Height1 through Height3 at first flower
  #	if(nrow(Repro)>=1) {
  #	Empty$Height1_flowering[target_row]<-Repro$Height1[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height2_flowering[target_row]<-Repro$Height2[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height3_flowering[target_row]<-Repro$Height3[1] }
  
  #		if(nrow(Repro)>=1) {
  #	Empty$Height4_flowering[target_row]<-Repro$Height4[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height5_flowering[target_row]<-Repro$Height5[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height6_flowering[target_row]<-Repro$Height6[1] }
  
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  #if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
  
  #Date_peak_flowering		
  #	max_fl<-max(Repro$Flower_number, na.rm=T)
  
  #	Empty$Date_peak_flowering[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
  
  
  #Flower_number_peak and Silique_number_peak and Stem_number_peak
  #	Empty$Flower_number_peak[target_row]<-max_fl
  
  #	Empty$Silique_number_peak[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Stem_number_peak[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
  
  
  #Heigothict1_peak througothic Heigothict5_peak
  #	Empty$Heigothict1_peak[target_row]<-Repro$Heigothict1[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict2_peak[target_row]<-Repro$Heigothict2[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict3_peak[target_row]<-Repro$Heigothict3[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict4_peak[target_row]<-Repro$Heigothict4[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict5_peak[target_row]<-Repro$Heigothict5[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict6_peak[target_row]<-Repro$Heigothict6[which(Repro$Flower_number==max_fl)][1]
  
  
  
  #} #Closes if statement
  
  
  
  
  ########Data for date of first SO - Note that this is the date of first SO observed after the last flower was observed
  
  #if(sum(Repro$Silique_number, na.rm=T)>=1) { #outer if statement
  
  #Extract the rowname for the last day a flower was observed, and then the rownumber of that rowname within sub (note that the rownames in sub will come from its name within the master, long-format data file)
  #	last_Fl<-rownames(tail(Repro[which(Repro$Flower_number>=1),], n=1))
  
  ##If there is a date of last flower, adjust the repro subset to reflect only data after the date of last flower
  #	if(length(last_Fl)>0) { #inner if statement
  
  ##Row number of date of last flower in sub
  #	last_Fl<-which(grepl(last_Fl, rownames(Repro)))
  
  ##Subset data to exclude censuses on or before the last flower was observed
  
  #	Repro<-Repro[-c(1:last_Fl),]
  
  #} #Closes inner if statement
  
  
  #Subset for censuses where a silique was observed (we have to do this second subset to exlude plants that became b or b, tco after flowering)
  #	Repro<-subset(Repro, Silique_number>=1)
  
  
  ##Now estimate the date of first silique
  #	Empty$Date_silique[target_row]<-Repro$Ordinal_date[1]
  
  
  #Stem_number_silique, Silique_number_silique, and Silique_length_silique
  #	Empty$Stem_number_silique[target_row]<-Repro$Stem_number[1] 
  
  #	Empty$Silique_number_silique[target_row]<-Repro$Silique_number[1]
  
  #	Empty$Silique_length_silique[target_row]<-Repro$Silique_length[1]
  
  
  #Height1_silique
  #	Empty$Height1_silique[target_row]<-Repro$Height1[1]
  
  #	Empty$Height2_silique[target_row]<-Repro$Height2[1]
  
  #	Empty$Height3_silique[target_row]<-Repro$Height3[1]
  
  #	Empty$Height4_silique[target_row]<-Repro$Height4[1]
  
  #	Empty$Height5_silique[target_row]<-Repro$Height5[1]
  
  #	Empty$Height6_silique[target_row]<-Repro$Height6[1]
  
  
  #Closing bracket from outer if statement
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  if (sum(sub$Silique_number, na.rm=T)>=1) { 
    
    Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
    
    ##If the number of failed siliques is negative, set to zero
    if(Empty$Failed_silique_number[target_row]<0) {
      
      Empty$Failed_silique_number[target_row]<-0
      
    } #closes inner if statement
 # }} #closes outer if statement


#Mature_length_siliques - note that the column numbers have to be updated
	Empty$Mature_length_siliques[target_row]<-sum(sub[,c(51,53,55,57,59,61,63,65,67,69,71)], na.rm=T)


#Total_flower_number
#	Empty$Total_flower_number[target_row]<-sum(sub$Flower_number, na.rm=T)

#Mature_entire_silique_length (when %PD is measured) - first vector has column numbers of silique lengths, second vector has column numbers of PD% - columns must be updated before calculating
#	Empty$Mature_PD[target_row]<-(sum(sub[,c(35,37,39,41,43,45,47,49,51,53,55,57,59,61,63)], na.rm=T))/sum(sub$Number_collected_siliques, na.rm=T)

	}} #closes outer if statement
#} ##Closes entire loop


G15000 <-subset(Empty, PlantID=="G15000")
View(G15000)

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "2021_2021gothic_summary.xlsx", sheetName="Summary2021", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 




#### gothic 2020 ####

#MAKE SURE YOU CLEAR THE GLOBAL ENV. THIS CALLS THE SAME VARIABLE NAMES AS IN GOTHIC 2021 COHORT
##Data
gothic <- read.csv("Field2021_gothic2020_merge.csv") #have to use this one bc i had to merge the genotype info with the data

sapply(gothic,class)

head(gothic,20)
tail(gothic)

gothic$Master_notes.x<-as.character(gothic$Master_notes.x)
gothic$Number_collected_siliques<-as.numeric(gothic$Number_collected_siliques)


reproductive<-subset(gothic, Phenology=="reproductive")
head(reproductive,20)




##Empty data set that we will populate with summary data
Empty <-read.csv("Field2021_gothic2020_empty.csv", header=T)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)
Empty $survival<-as.numeric(Empty $survival)
Empty $bolted<-as.numeric(Empty $bolted)
Empty $flowered<-as.numeric(Empty $flowered)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(gothic$PlantID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(gothic $PlantID)[i]
  sub<-subset(gothic, PlantID==target_ID)
  
  #        ##If the subset has less than 13 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==13)==F)
  { print(sub$PlantID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(gothic$PlantID))) {
  
  ##Subset the data for plantID i
  target_ID<-unique(gothic$PlantID)[i]
  sub<-subset(gothic, PlantID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$PlantID==target_ID),]))
  
  
  #exp_survival - Census 12 was the last full census, alive=1 dead=0
  ifelse((sub$Status[13]=="alive"), Empty$survival[target_row]<-1, Empty$survival[target_row]<-0)
  
  ## SKIP ####
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census or if they had leaves counted during the damage census (note that we can only ask if there are >= 1 leaf if the element is not NA)
  #	if(Empty$Season_survival[target_row]==1 | (Empty$Leaf_Number[target_row]>=1 & is.na(Empty$Leaf_Number[target_row])==F)) {
  #		Empty$Overwinter_survival[target_row]<-1
  #	}
  
  ##I measured leaf damage for only a few censuses, so we will put them all here	
  
  Empty$LAR_2[target_row]<-sub$LAR[2] 
  Empty$LAR_3[target_row]<-sub$LAR[3] 
  Empty$LAR_5[target_row]<-sub$LAR[5] #jill did this one
  Empty$LAR_7[target_row]<-sub$LAR[7] 
  
  
  
  
  ## Master notes
  Empty$Master_notes[target_row]<-sub$Master_notes.x[18] 
  
  ## vegetative cover not measured
  #Empty$vegetative_cover[target_row]<-sub$X..Vegetative.cover[7] 
  
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  #if(Empty$survival[target_row]==0) { #Outer if statement
  
  #Subset data where Status==Dead	
  #First_death<-subset(sub, Status=="dead,check again")
  
  #Extract row name of last census where status==Alive
  #Last_alive<-rownames(tail(sub[which(sub$Status=="alive"),], n=1))
  
  ##Adjust First_death data set to include only data 	after the last Alive status was observed
  #if(length(Last_alive)>=1) { #Inner if statement
  
  ##Adjust subset of First_death
  #First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
  
  #} #Closes inner if statement
  
  ##Assign the date of first death
  #	Empty$Date_first_death_observed[target_row]<-First_death$Date[1] 
  
  #} #Closes outer if statement
  
  
  ###SKIP###	
  ##Gopher damage or Censor - was coded as GopherDamage or Censor in the Exclude column - sometimes a plant marked as being exposed to gopher damage doesn't die - for that reason, this has to meet conditions that it was censused as gopher damaged at least once and that we didn't observe it as alive in the last census
  #	if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive")) {
  #	Empty$Exclude[target_row]<-"GopherDamage"
  #	}
  
  #	if(any(sub$Exclude=="Censor") & (sub$Status[8]!="Alive")) {
  #	Empty$Exclude[target_row]<-"Censor"
  #	}
  
  #	if(is.na(Empty$Exclude[target_row])==T) {
  #	Empty$Exclude[target_row]<-"Include"
  #	}
  
  
  ##Date_first_gopher, again only noted if the plant was not marked as alive in the last census
  ##This code is a little crazy because we are also calculating the date a plant was last alive prior to gopher damage and whether the living plant died prior to gopher damage
  #if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive" )) {
  
  #	Empty$Date_first_gopher[target_row]<-sub$Ordinal_date[which(sub$Exclude=="GopherDamage")][1] 
  
  ##Last date of being alive prior to gopher damage
  #	Last_alive<-tail(sub$Ordinal_date[which(sub$Status=="Alive" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)
  
  #Last date of being dead prior to gopher damage
  #	Last_dead<-tail(sub$Ordinal_date[which(sub$Status=="Dead" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)	
  
  #If the plant has both a Last_alive and Last_dead date
  #	if(length(Last_alive)>=1 & length(Last_dead)>=1) {
  
  #If the plant was alive before it was dead and if it died before gopher was first observed - have to separate this if statement from the last because you can't ask if something is less than something else if never occured
  #	if(Last_alive<Last_dead & Last_dead<Empty$Date_first_gopher[target_row]) {
  
  #Include plant
  #		Empty$Exclude[target_row]<-"Include"
  
  #	} } #closes inner if statements
  
  ##Lastly, if plant died prior to gopher damage, list it as a natural death
  #	if(Empty$Date_first_death_observed[target_row]<Empty$Date_first_gopher[target_row]) {
  
  #		Empty$Exclude[target_row]<-"Include"
  
  #	} #Closes inner if statement
  
  #} #Closes outer if statement
  
  
  
  
  ## issues here####	
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  
  ifelse((sum(sub$Phenology=="reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$bolted[target_row]<-1, Empty$bolted[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$flowered[target_row]<-1, Empty$flowered[target_row]<-0)
  
  
  #Subset reproductive data #also doesnot work
  Repro<-subset(sub, Phenology=="reproductive" & (Flower_number>=1 | Silique_number>=1))
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering[target_row]<-Repro$Ordinal_date[1] }
  
  
  #Flower_color - note the nested ifelse statements - this is needed to assign a 1 for plants with pigmented flowers, a 0 to plants with white flowers, and an NA to plants that never flowered (or for which we never observed flowers) 
  #	if(any(Repro$Flower_color=="Pigmented")) {
  #	Empty$Flower_color[target_row]<-1 }
  #	
  #	if(all(Repro$Flower_color!="Pigmented") & sum(Repro$Flower_number, na.rm=T)>=1) {
  #	Empty$Flower_color[target_row]<-0 }
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  #	if(nrow(Repro)>=1) {
  #	Empty$Stem_number_flowering[target_row]<-Repro$Stem_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Flower_number_flowering[target_row]<-Repro$Flower_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Silique_number_flowering[target_row]<-Repro$Silique_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Silique_length_flowering[target_row]<-Repro$Silique_length[1] }
  
  
  #Height1 through Height3 at first flower
  #	if(nrow(Repro)>=1) {
  #	Empty$Height1_flowering[target_row]<-Repro$Height1[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height2_flowering[target_row]<-Repro$Height2[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height3_flowering[target_row]<-Repro$Height3[1] }
  
  #		if(nrow(Repro)>=1) {
  #	Empty$Height4_flowering[target_row]<-Repro$Height4[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height5_flowering[target_row]<-Repro$Height5[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height6_flowering[target_row]<-Repro$Height6[1] }
  
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  #if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
  
  #Date_peak_flowering		
  #	max_fl<-max(Repro$Flower_number, na.rm=T)
  
  #	Empty$Date_peak_flowering[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
  
  
  #Flower_number_peak and Silique_number_peak and Stem_number_peak
  #	Empty$Flower_number_peak[target_row]<-max_fl
  
  #	Empty$Silique_number_peak[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Stem_number_peak[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
  
  
  #Heigothict1_peak througothic Heigothict5_peak
  #	Empty$Heigothict1_peak[target_row]<-Repro$Heigothict1[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict2_peak[target_row]<-Repro$Heigothict2[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict3_peak[target_row]<-Repro$Heigothict3[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict4_peak[target_row]<-Repro$Heigothict4[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict5_peak[target_row]<-Repro$Heigothict5[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heigothict6_peak[target_row]<-Repro$Heigothict6[which(Repro$Flower_number==max_fl)][1]
  
  
  
  #} #Closes if statement
  
  
  
  
  ########Data for date of first SO - Note that this is the date of first SO observed after the last flower was observed
  
  #if(sum(Repro$Silique_number, na.rm=T)>=1) { #outer if statement
  
  #Extract the rowname for the last day a flower was observed, and then the rownumber of that rowname within sub (note that the rownames in sub will come from its name within the master, long-format data file)
  #	last_Fl<-rownames(tail(Repro[which(Repro$Flower_number>=1),], n=1))
  
  ##If there is a date of last flower, adjust the repro subset to reflect only data after the date of last flower
  #	if(length(last_Fl)>0) { #inner if statement
  
  ##Row number of date of last flower in sub
  #	last_Fl<-which(grepl(last_Fl, rownames(Repro)))
  
  ##Subset data to exclude censuses on or before the last flower was observed
  
  #	Repro<-Repro[-c(1:last_Fl),]
  
  #} #Closes inner if statement
  
  
  #Subset for censuses where a silique was observed (we have to do this second subset to exlude plants that became b or b, tco after flowering)
  #	Repro<-subset(Repro, Silique_number>=1)
  
  
  ##Now estimate the date of first silique
  #	Empty$Date_silique[target_row]<-Repro$Ordinal_date[1]
  
  
  #Stem_number_silique, Silique_number_silique, and Silique_length_silique
  #	Empty$Stem_number_silique[target_row]<-Repro$Stem_number[1] 
  
  #	Empty$Silique_number_silique[target_row]<-Repro$Silique_number[1]
  
  #	Empty$Silique_length_silique[target_row]<-Repro$Silique_length[1]
  
  
  #Height1_silique
  #	Empty$Height1_silique[target_row]<-Repro$Height1[1]
  
  #	Empty$Height2_silique[target_row]<-Repro$Height2[1]
  
  #	Empty$Height3_silique[target_row]<-Repro$Height3[1]
  
  #	Empty$Height4_silique[target_row]<-Repro$Height4[1]
  
  #	Empty$Height5_silique[target_row]<-Repro$Height5[1]
  
  #	Empty$Height6_silique[target_row]<-Repro$Height6[1]
  
  
  #Closing bracket from outer if statement
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  if (sum(sub$Silique_number, na.rm=T)>=1) { 
    
    Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
    
    ##If the number of failed siliques is negative, set to zero
    if(Empty$Failed_silique_number[target_row]<0) {
      
      Empty$Failed_silique_number[target_row]<-0
      
    } #closes inner if statement
    # }} #closes outer if statement
    
    
    #Mature_length_siliques - note that the column numbers have to be updated
    #	Empty$Mature_length_siliques[target_row]<-sum(sub[,c(51,53,55,57,59,61,63,65,67,69,71)], na.rm=T)
    
    
    #Total_flower_number
    #	Empty$Total_flower_number[target_row]<-sum(sub$Flower_number, na.rm=T)
    
    #Mature_entire_silique_length (when %PD is measured) - first vector has column numbers of silique lengths, second vector has column numbers of PD% - columns must be updated before calculating
    #	Empty$Mature_PD[target_row]<-(sum(sub[,c(35,37,39,41,43,45,47,49,51,53,55,57,59,61,63)], na.rm=T))/sum(sub$Number_collected_siliques, na.rm=T)
    
  }} #closes outer if statement
#} ##Closes entire loop


G13740 <-subset(Empty, PlantID=="G13740")
View(G13740)

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "2021_2020gothic_summary.xlsx", sheetName="Summary2021", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 


#### schofield 2020 ####

#MAKE SURE YOU CLEAR THE GLOBAL ENV.T
##Data
schofield <- read.csv("Field2021_schofield_merge.csv") #have to use this one bc i had to merge the genotype info with the data

sapply(schofield,class)

head(schofield,20)
tail(schofield)

schofield$Master_notes<-as.character(schofield$Master_notes)
schofield$Number_collected_siliques<-as.numeric(schofield$Number_collected_siliques)


reproductive<-subset(schofield, Phenology=="reproductive")
head(reproductive,20)




##Empty data set that we will populate with summary data
Empty <-read.csv("Field2021_schofield_empty.csv", header=T)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)
Empty $survival<-as.numeric(Empty $survival)
Empty $bolted<-as.numeric(Empty $bolted)
Empty $flowered<-as.numeric(Empty $flowered)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(schofield$PlantID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(schofield $PlantID)[i]
  sub<-subset(schofield, PlantID==target_ID)
  
  #        ##If the subset has less than 8 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==8)==F)
  { print(sub$PlantID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(schofield$PlantID))) {
  
  ##Subset the data for plantID i
  target_ID<-unique(schofield$PlantID)[i]
  sub<-subset(schofield, PlantID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$PlantID==target_ID),]))
  
  
  #exp_survival - Census 8 was the last full census, alive=1 dead=0
  ifelse((sub$Status[8]=="alive"), Empty$survival[target_row]<-1, Empty$survival[target_row]<-0)
  
  ## SKIP ####
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census or if they had leaves counted during the damage census (note that we can only ask if there are >= 1 leaf if the element is not NA)
  #	if(Empty$Season_survival[target_row]==1 | (Empty$Leaf_Number[target_row]>=1 & is.na(Empty$Leaf_Number[target_row])==F)) {
  #		Empty$Overwinter_survival[target_row]<-1
  #	}
  
  ##the file only has a census for the 8 that jill did, I need to come back and reimport the other values will put them all here	
  
  Empty$LAR_8[target_row]<-sub$LAR[8] #jill did this one

  
  
  
  
  ## Master notes
  Empty$Master_notes[target_row]<-sub$Master_notes[8] 
  
  ## vegetative cover 
  Empty$vegetative_cover[target_row]<-sub$X..cover[6] 
  
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  #if(Empty$survival[target_row]==0) { #Outer if statement
  
  #Subset data where Status==Dead	
  #First_death<-subset(sub, Status=="dead,check again")
  
  #Extract row name of last census where status==Alive
  #Last_alive<-rownames(tail(sub[which(sub$Status=="alive"),], n=1))
  
  ##Adjust First_death data set to include only data 	after the last Alive status was observed
  #if(length(Last_alive)>=1) { #Inner if statement
  
  ##Adjust subset of First_death
  #First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
  
  #} #Closes inner if statement
  
  ##Assign the date of first death
  #	Empty$Date_first_death_observed[target_row]<-First_death$Date[1] 
  
  #} #Closes outer if statement
  
  
  ###SKIP###	
  ##Gopher damage or Censor - was coded as GopherDamage or Censor in the Exclude column - sometimes a plant marked as being exposed to gopher damage doesn't die - for that reason, this has to meet conditions that it was censused as gopher damaged at least once and that we didn't observe it as alive in the last census
  #	if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive")) {
  #	Empty$Exclude[target_row]<-"GopherDamage"
  #	}
  
  #	if(any(sub$Exclude=="Censor") & (sub$Status[8]!="Alive")) {
  #	Empty$Exclude[target_row]<-"Censor"
  #	}
  
  #	if(is.na(Empty$Exclude[target_row])==T) {
  #	Empty$Exclude[target_row]<-"Include"
  #	}
  
  
  ##Date_first_gopher, again only noted if the plant was not marked as alive in the last census
  ##This code is a little crazy because we are also calculating the date a plant was last alive prior to gopher damage and whether the living plant died prior to gopher damage
  #if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive" )) {
  
  #	Empty$Date_first_gopher[target_row]<-sub$Ordinal_date[which(sub$Exclude=="GopherDamage")][1] 
  
  ##Last date of being alive prior to gopher damage
  #	Last_alive<-tail(sub$Ordinal_date[which(sub$Status=="Alive" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)
  
  #Last date of being dead prior to gopher damage
  #	Last_dead<-tail(sub$Ordinal_date[which(sub$Status=="Dead" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)	
  
  #If the plant has both a Last_alive and Last_dead date
  #	if(length(Last_alive)>=1 & length(Last_dead)>=1) {
  
  #If the plant was alive before it was dead and if it died before gopher was first observed - have to separate this if statement from the last because you can't ask if something is less than something else if never occured
  #	if(Last_alive<Last_dead & Last_dead<Empty$Date_first_gopher[target_row]) {
  
  #Include plant
  #		Empty$Exclude[target_row]<-"Include"
  
  #	} } #closes inner if statements
  
  ##Lastly, if plant died prior to gopher damage, list it as a natural death
  #	if(Empty$Date_first_death_observed[target_row]<Empty$Date_first_gopher[target_row]) {
  
  #		Empty$Exclude[target_row]<-"Include"
  
  #	} #Closes inner if statement
  
  #} #Closes outer if statement
  
  
  
  
  ## issues here####	
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  
  ifelse((sum(sub$Phenology=="reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$bolted[target_row]<-1, Empty$bolted[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$flowered[target_row]<-1, Empty$flowered[target_row]<-0)
  
  
  #Subset reproductive data #also doesnot work
  Repro<-subset(sub, Phenology=="reproductive" & (Flower_number>=1 | Silique_number>=1))
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering[target_row]<-Repro$Ordinal_date[1] }
  
  
  #Flower_color - note the nested ifelse statements - this is needed to assign a 1 for plants with pigmented flowers, a 0 to plants with white flowers, and an NA to plants that never flowered (or for which we never observed flowers) 
  #	if(any(Repro$Flower_color=="Pigmented")) {
  #	Empty$Flower_color[target_row]<-1 }
  #	
  #	if(all(Repro$Flower_color!="Pigmented") & sum(Repro$Flower_number, na.rm=T)>=1) {
  #	Empty$Flower_color[target_row]<-0 }
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  #	if(nrow(Repro)>=1) {
  #	Empty$Stem_number_flowering[target_row]<-Repro$Stem_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Flower_number_flowering[target_row]<-Repro$Flower_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Silique_number_flowering[target_row]<-Repro$Silique_number[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Silique_length_flowering[target_row]<-Repro$Silique_length[1] }
  
  
  #Height1 through Height3 at first flower
  #	if(nrow(Repro)>=1) {
  #	Empty$Height1_flowering[target_row]<-Repro$Height1[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height2_flowering[target_row]<-Repro$Height2[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height3_flowering[target_row]<-Repro$Height3[1] }
  
  #		if(nrow(Repro)>=1) {
  #	Empty$Height4_flowering[target_row]<-Repro$Height4[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height5_flowering[target_row]<-Repro$Height5[1] }
  
  #	if(nrow(Repro)>=1) {
  #	Empty$Height6_flowering[target_row]<-Repro$Height6[1] }
  
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  #if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
  
  #Date_peak_flowering		
  #	max_fl<-max(Repro$Flower_number, na.rm=T)
  
  #	Empty$Date_peak_flowering[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
  
  
  #Flower_number_peak and Silique_number_peak and Stem_number_peak
  #	Empty$Flower_number_peak[target_row]<-max_fl
  
  #	Empty$Silique_number_peak[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Stem_number_peak[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
  
  
  #Heischofieldt1_peak throuschofield Heischofieldt5_peak
  #	Empty$Heischofieldt1_peak[target_row]<-Repro$Heischofieldt1[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heischofieldt2_peak[target_row]<-Repro$Heischofieldt2[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heischofieldt3_peak[target_row]<-Repro$Heischofieldt3[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heischofieldt4_peak[target_row]<-Repro$Heischofieldt4[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heischofieldt5_peak[target_row]<-Repro$Heischofieldt5[which(Repro$Flower_number==max_fl)][1]
  
  #	Empty$Heischofieldt6_peak[target_row]<-Repro$Heischofieldt6[which(Repro$Flower_number==max_fl)][1]
  
  
  
  #} #Closes if statement
  
  
  
  
  ########Data for date of first SO - Note that this is the date of first SO observed after the last flower was observed
  
  #if(sum(Repro$Silique_number, na.rm=T)>=1) { #outer if statement
  
  #Extract the rowname for the last day a flower was observed, and then the rownumber of that rowname within sub (note that the rownames in sub will come from its name within the master, long-format data file)
  #	last_Fl<-rownames(tail(Repro[which(Repro$Flower_number>=1),], n=1))
  
  ##If there is a date of last flower, adjust the repro subset to reflect only data after the date of last flower
  #	if(length(last_Fl)>0) { #inner if statement
  
  ##Row number of date of last flower in sub
  #	last_Fl<-which(grepl(last_Fl, rownames(Repro)))
  
  ##Subset data to exclude censuses on or before the last flower was observed
  
  #	Repro<-Repro[-c(1:last_Fl),]
  
  #} #Closes inner if statement
  
  
  #Subset for censuses where a silique was observed (we have to do this second subset to exlude plants that became b or b, tco after flowering)
  #	Repro<-subset(Repro, Silique_number>=1)
  
  
  ##Now estimate the date of first silique
  #	Empty$Date_silique[target_row]<-Repro$Ordinal_date[1]
  
  
  #Stem_number_silique, Silique_number_silique, and Silique_length_silique
  #	Empty$Stem_number_silique[target_row]<-Repro$Stem_number[1] 
  
  #	Empty$Silique_number_silique[target_row]<-Repro$Silique_number[1]
  
  #	Empty$Silique_length_silique[target_row]<-Repro$Silique_length[1]
  
  
  #Height1_silique
  #	Empty$Height1_silique[target_row]<-Repro$Height1[1]
  
  #	Empty$Height2_silique[target_row]<-Repro$Height2[1]
  
  #	Empty$Height3_silique[target_row]<-Repro$Height3[1]
  
  #	Empty$Height4_silique[target_row]<-Repro$Height4[1]
  
  #	Empty$Height5_silique[target_row]<-Repro$Height5[1]
  
  #	Empty$Height6_silique[target_row]<-Repro$Height6[1]
  
  
  #Closing bracket from outer if statement
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  if (sum(sub$Silique_number, na.rm=T)>=1) { 
    
    Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
    
    ##If the number of failed siliques is negative, set to zero
    if(Empty$Failed_silique_number[target_row]<0) {
      
      Empty$Failed_silique_number[target_row]<-0
      
    } #closes inner if statement
    # }} #closes outer if statement
    
    
    #Mature_length_siliques - note that the column numbers have to be updated
    #	Empty$Mature_length_siliques[target_row]<-sum(sub[,c(51,53,55,57,59,61,63,65,67,69,71)], na.rm=T)
    
    
    #Total_flower_number
    #	Empty$Total_flower_number[target_row]<-sum(sub$Flower_number, na.rm=T)
    
    #Mature_entire_silique_length (when %PD is measured) - first vector has column numbers of silique lengths, second vector has column numbers of PD% - columns must be updated before calculating
    #	Empty$Mature_PD[target_row]<-(sum(sub[,c(35,37,39,41,43,45,47,49,51,53,55,57,59,61,63)], na.rm=T))/sum(sub$Number_collected_siliques, na.rm=T)
    
  }} #closes outer if statement
#} ##Closes entire loop


S10599 <-subset(Empty, PlantID=="S10599")
View(S10599)

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "2021_schofield_summary.xlsx", sheetName="Summary2021", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 
