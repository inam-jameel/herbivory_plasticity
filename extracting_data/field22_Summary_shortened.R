setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/field/2022_files/census_files")

getwd()

library(writexl)


###### gothic 2021 cohort ######


##Data
gothic<-read.csv("2021_ginam_03oct22.csv", header=T)
sapply(gothic,class)

head(gothic,20)
tail(gothic)

gothic$Master_notes<-as.character(gothic$Master_notes)
gothic$Flower_number <as.integer(gothic$Flower_number)

reproductive<-subset(gothic, State=="Reproductive")
head(reproductive,20)
tail(reproductive,20)




##Empty data set that we will populate with summary data


Empty_input<-read.csv("inam_gothic2022_empty21.csv", header=T)
Empty <- Empty_input[order(Empty_input$PlantID),]
row.names(Empty) <- 1:nrow(Empty)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(gothic $PlantID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(gothic $PlantID)[i]
  sub<-subset(gothic, PlantID==target_ID)
  
  #        ##If the subset has less than 21 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==10)==F)
  { print(sub$PlantID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(gothic$PlantID))){ 
  
  ##Subset the data for plantID i
  target_ID<-unique(gothic$PlantID)[i]
  sub<-subset(gothic, PlantID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$PlantID==target_ID),]))
  
  #overwinter_survival - Census 1 was the first full census, alive=1 dead=0
  ifelse((sub$Status[1]=="Alive"), Empty$Overwinter_survival_2022[target_row]<-1, Empty$Overwinter_survival_2022[target_row]<-0)
  
  
  #season_survival - Census 11 was the last full census, alive=1 dead=0
  ifelse((sub$Status[10]=="Alive"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  
  #season_Include - Census 10 was the last full census, alive=1 dead=0
  #ifelse((sub$Include[10]=="Alive"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  Empty$Exclude_2022[target_row]<-sub$Include[8] 
  
  #error here
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census
  
  if(Empty$Season_survival_2022[target_row]==1) {
    Empty$Overwinter_survival_2022[target_row]<-1
  }
  
  
  # leaf damage not in file
  
  #Empty$num_total3_2022[target_row]<-sub$Leaf_number[3] 
  #Empty$num_damaged3_2022[target_row]<-sub$Damaged_leaves[3] 
  #Empty$total_damage3_2022[target_row]<-sub$total_damage[3]
  #Empty$avg_damage3_2022[target_row]<-sub$average_damage[3]
  #Empty$LAR_3_2022[target_row]<-sub$LAR[3] 
  
  #Empty$num_total7_2022[target_row]<-sub$Leaf_number[7] 
  #Empty$num_damaged7_2022[target_row]<-sub$Damaged_leaves[7] 
  #Empty$total_damage7_2022[target_row]<-sub$total_damage[7]
  #Empty$avg_damage7_2022[target_row]<-sub$average_damage[7]
  #Empty$LAR_7_2022[target_row]<-sub$LAR[7] 
  
  
  ## Master notes
  Empty$Master_notes_2022[target_row]<-sub$Master_notes[2] 
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  if(Empty$Season_survival_2022[target_row]==0) { #Outer if statement
    
    #Subset data where Status==Dead	
    First_death<-subset(sub, Status=="Dead")
    
    #Extract row name of last census where status==Alive
    Last_alive<-rownames(tail(sub[which(sub$Status=="Alive"),], n=1))
    
    ##Adjust First_death data set to include only data 	after the last Alive status was observed
    if(length(Last_alive)>=1) { #Inner if statement
      
      ##Adjust subset of First_death
      First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
      
    } #Closes inner if statement
    
    ##Assign the date of first death
    Empty$Date_first_death_observed_2022[target_row]<-First_death$Ordinal_date[1] 
    
  } #Closes outer if statement
  
  ##Gopher damage or Censor - was coded as GopherDamage or Censor in the Exclude column - sometimes a plant marked as being exposed to gopher damage doesn't die - for that reason, this has to meet conditions that it was censused as gopher damaged at least once and that we didn't observe it as alive in the last census
  #if(any(sub$Exclude[10]=="GopherDamage") & (sub$Status[10]!="Alive")) {
  #  Empty$Exclude_2022[target_row]<-"GopherDamage"
  #}
  
  #if(any(sub$Exclude[10]=="Censor") & (sub$Status[10]!="Alive")) {
  #  Empty$Exclude_2022[target_row]<-"Censor"
  #}
  
  #if(is.na(Empty$Exclude_2022[target_row])==T) {
  #  Empty$Exclude_2022[target_row]<-"Include"
  #}
  
  
  
  ##Date_first_gopher, again only noted if the plant was not marked as alive in the last census
  ##This code is a little crazy because we are also calculating the date a plant was last alive prior to gopher damage and whether the living plant died prior to gopher damage
  #if(any(sub$Exclude[10]=="GopherDamage") & (sub$Status[10]!="Alive" )) {
  
  #  Empty$Date_first_gopher_2022[target_row]<-sub$Ordinal_date[which(sub$Exclude=="GopherDamage")][1] 
  
  ##Last date of being alive prior to gopher damage
  #  Last_alive<-tail(sub$Ordinal_date[which(sub$Status=="Alive" & sub$Ordinal_date<Empty$Date_first_gopher_2022[target_row])], n=1)
  
  #Last date of being dead prior to gopher damage
  #  Last_dead<-tail(sub$Ordinal_date[which(sub$Status=="Dead" & sub$Ordinal_date<Empty$Date_first_gopher_2022[target_row])], n=1)	
  
  #If the plant has both a Last_alive and Last_dead date
  #  if(length(Last_alive)>=1 & length(Last_dead)>=1) {
  
  #If the plant was alive before it was dead and if it died before gopher was first observed - have to separate this if statement from the last because you can't ask if something is less than something else if never occured
  #    if(Last_alive<Last_dead & Last_dead<Empty$Date_first_gopher_2022[target_row]) {
  
  #Include plant
  #      Empty$Exclude_2022[target_row]<-"Include"
  
  #    } } #closes inner if statements
  
  ##Lastly, if plant died prior to gopher damage, list it as a natural death
  #  if(Empty$Date_first_death_observed_2022[target_row]<Empty$Date_first_gopher_2022[target_row]) {
  
  #    Empty$Exclude_2022[target_row]<-"Include"
  
  #  } #Closes inner if statement
  
  #} #Closes outer if statement
  
  
  
  
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  ifelse((sum(sub$State=="Reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Bolted_2022[target_row]<-1, Empty$Bolted_2022[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered_2022[target_row]<-1, Empty$Flowered_2022[target_row]<-0)
  
  ####Subset reproductive data
  Repro<-subset(sub, State=="Reproductive" & (Flower_number>=1 | Silique_number>=1)) #this filters out the rows where the plant was reproductive, but is missing flower or silique number
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering_2022[target_row]<-Repro$Ordinal_date[1] }
  
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  if(nrow(Repro)>=1) {
    Empty$Stem_number_flowering_2022[target_row]<-Repro$Stem_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Flower_number_flowering_2022[target_row]<-Repro$Flower_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_number_flowering_2022[target_row]<-Repro$Silique_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_length_flowering_2022[target_row]<-Repro$Silique_length[1] 
    #<-Repro$Silique_length[1] 
  }
  
  
  
  
  #Height1 through Height3 at first flower
  if(nrow(Repro)>=1) {
    Empty$Height1_flowering_2022[target_row]<-Repro$Height1[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height2_flowering_2022[target_row]<-Repro$Height2[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height3_flowering_2022[target_row]<-Repro$Height3[1] }
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
    
    #Date_peak_flowering		
    max_fl<-max(Repro$Flower_number, na.rm=T)
    
    Empty$Date_peak_flowering_2022[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
    
    
    #Flower_number_peak and Silique_number_peak and Stem_number_peak
    Empty$Flower_number_peak_2022[target_row]<-max_fl
    
    Empty$Silique_number_peak_2022[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Stem_number_peak_2022[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
    
    
    #Height1_peak through Height3_peak
    Empty$Height1_peak_2022[target_row]<-Repro$Height1[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height2_peak_2022[target_row]<-Repro$Height2[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height3_peak_2022[target_row]<-Repro$Height3[which(Repro$Flower_number==max_fl)][1]
    
    
  } #Closes if statement
  
  
  
  
  ########Data for date of first SO - Note that this is the date of first SO observed after the last flower was observed
  
  if(sum(Repro$Silique_number, na.rm=T)>=1) { #outer if statement
    
    #Extract the rowname for the last day a flower was observed, and then the rownumber of that rowname within sub (note that the rownames in sub will come from its name within the master, long-format data file)
    last_Fl<-rownames(tail(Repro[which(Repro$Flower_number>=1),], n=1))
    
    ##If there is a date of last flower, adjust the repro subset to reflect only data after the date of last flower
    if(length(last_Fl)>0) { #inner if statement
      
      ##Row number of date of last flower in sub
      last_Fl<-which(grepl(last_Fl, rownames(Repro)))
      
      ##Subset data to exclude censuses on or before the last flower was observed
      
      Repro<-Repro[-c(1:last_Fl),]
      
    } #Closes inner if statement
    
    
    #Subset for censuses where a silique was observed (we have to do this second subset to exlude plants that became b or b, tco after flowering)
    Repro<-subset(Repro, Silique_number>=1)
    
    
    ##Now estimate the date of first silique
    Empty$Date_silique_2022[target_row]<-Repro$Ordinal_date[1]
    
    
    #Stem_number_silique, Silique_number_silique, and Silique_length_silique
    Empty$Stem_number_silique_2022[target_row]<-Repro$Stem_number[1] 
    
    Empty$Silique_number_silique_2022[target_row]<-Repro$Silique_number[1]
    
    Empty$Silique_length_silique_2022[target_row]<-Repro$Silique_length[1]
    
    
    #Height3_silique
    Empty$Height1_silique_2022[target_row]<-Repro$Height1[1]
    
    Empty$Height2_silique_2022[target_row]<-Repro$Height2[1]
    
    Empty$Height3_silique_2022[target_row]<-Repro$Height3[1]
    
    
    
  } #Closing bracket from outer if statement
  
  
  
  
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number_2022[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number 
  Empty$Failed_silique_number_2022[target_row]<-sum(sub$Failed_silique_number, na.rm=T)
  
  ##Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  #if (sum(sub$Silique_number, na.rm=T)>=1) { 
  
  # Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
  
  ##If the number of failed siliques is negative, set to zero
  #if(Empty$Failed_silique_number[target_row]<0) {
  
  # Empty$Failed_silique_number[target_row]<-0
  
  #Mature_length_siliques - note that the column numbers have to be updated
  Empty$Mature_length_siliques_2022[target_row]<-sum(sub[,c(62,
                                                            64,
                                                            66,
                                                            68,
                                                            70,
                                                            72,
                                                            74,
                                                            76,
                                                            78,
                                                            80,
                                                            82,
                                                            84,
                                                            86,
                                                            88)], na.rm=T)
  
} #closes inner if statement



check1 <-subset(Empty, PlantID=="G16375")
check2 <-subset(Empty, PlantID=="G16316")

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "gothic2021_2022_summary.xlsx", sheetName="Summary2022", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 

###### gothic 2022 cohort ######


##Data
gothic<-read.csv("2022_ginam_census22.csv", header=T)
sapply(gothic,class)

head(gothic,20)
tail(gothic)

gothic$Master_notes<-as.character(gothic$Master_notes)
gothic$Flower_number <as.integer(gothic$Flower_number)

reproductive<-subset(gothic, State=="Reproductive")
head(reproductive,20)
tail(reproductive,20)




##Empty data set that we will populate with summary data


Empty_input<-read.csv("2022_ginam_empty22.csv", header=T)
Empty <- Empty_input[order(Empty_input$PlantID),]
row.names(Empty) <- 1:nrow(Empty)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(gothic $PlantID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(gothic $PlantID)[i]
  sub<-subset(gothic, PlantID==target_ID)
  
  #        ##If the subset has less than 21 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==9)==F)
  { print(sub$PlantID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(gothic$PlantID))){ 
  
  ##Subset the data for plantID i
  target_ID<-unique(gothic$PlantID)[i]
  sub<-subset(gothic, PlantID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$PlantID==target_ID),]))
  
  #overwinter_survival - Census 1 was the first full census, alive=1 dead=0
  ifelse((sub$Status[1]=="Alive"), Empty$Overwinter_survival_2022[target_row]<-1, Empty$Overwinter_survival_2022[target_row]<-0)
  
  
  #season_survival - Census 11 was the last full census, alive=1 dead=0
  ifelse((sub$Status[8]=="Alive"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  
  #season_Include - Census 1 has the plants that died due to trasnplant shock , alive=1 dead=0
  #ifelse((sub$Include[10]=="Alive"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  Empty$Exclude_2022[target_row]<-sub$Include[1] 
  
  #error here
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census
  
  if(Empty$Season_survival_2022[target_row]==1) {
    Empty$Overwinter_survival_2022[target_row]<-1
  }
  
  
  # leaf damage was assessed on  censuses 1,3,5,7,8
  
  Empty$LAR_1_2022[target_row]<-sub$LAR[1] 
  Empty$LAR_2_2022[target_row]<-sub$LAR[3] 
  Empty$LAR_3_2022[target_row]<-sub$LAR[5] 
  Empty$LAR_4_2022[target_row]<-sub$LAR[7] 
  Empty$LAR_5_2022[target_row]<-sub$LAR[8] 
  
  ## Master notes
  Empty$Master_notes_2022[target_row]<-sub$Master_notes[2] 
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  if(Empty$Season_survival_2022[target_row]==0) { #Outer if statement
    
    #Subset data where Status==Dead	
    First_death<-subset(sub, Status=="Dead")
    
    #Extract row name of last census where status==Alive
    Last_alive<-rownames(tail(sub[which(sub$Status=="Alive"),], n=1))
    
    ##Adjust First_death data set to include only data 	after the last Alive status was observed
    if(length(Last_alive)>=1) { #Inner if statement
      
      ##Adjust subset of First_death
      First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
      
    } #Closes inner if statement
    
    ##Assign the date of first death
    Empty$Date_first_death_observed_2022[target_row]<-First_death$Ordinal_date[1] 
    
  } #Closes outer if statement
  
  ##Gopher damage or Censor - was coded as GopherDamage or Censor in the Exclude column - sometimes a plant marked as being exposed to gopher damage doesn't die - for that reason, this has to meet conditions that it was censused as gopher damaged at least once and that we didn't observe it as alive in the last census
  #if(any(sub$Exclude[10]=="GopherDamage") & (sub$Status[10]!="Alive")) {
  #  Empty$Exclude_2022[target_row]<-"GopherDamage"
  #}
  
  #if(any(sub$Exclude[10]=="Censor") & (sub$Status[10]!="Alive")) {
  #  Empty$Exclude_2022[target_row]<-"Censor"
  #}
  
  #if(is.na(Empty$Exclude_2022[target_row])==T) {
  #  Empty$Exclude_2022[target_row]<-"Include"
  #}
  
  
  
  ##Date_first_gopher, again only noted if the plant was not marked as alive in the last census
  ##This code is a little crazy because we are also calculating the date a plant was last alive prior to gopher damage and whether the living plant died prior to gopher damage
  #if(any(sub$Exclude[10]=="GopherDamage") & (sub$Status[10]!="Alive" )) {
  
  #  Empty$Date_first_gopher_2022[target_row]<-sub$Ordinal_date[which(sub$Exclude=="GopherDamage")][1] 
  
  ##Last date of being alive prior to gopher damage
  #  Last_alive<-tail(sub$Ordinal_date[which(sub$Status=="Alive" & sub$Ordinal_date<Empty$Date_first_gopher_2022[target_row])], n=1)
  
  #Last date of being dead prior to gopher damage
  #  Last_dead<-tail(sub$Ordinal_date[which(sub$Status=="Dead" & sub$Ordinal_date<Empty$Date_first_gopher_2022[target_row])], n=1)	
  
  #If the plant has both a Last_alive and Last_dead date
  #  if(length(Last_alive)>=1 & length(Last_dead)>=1) {
  
  #If the plant was alive before it was dead and if it died before gopher was first observed - have to separate this if statement from the last because you can't ask if something is less than something else if never occured
  #    if(Last_alive<Last_dead & Last_dead<Empty$Date_first_gopher_2022[target_row]) {
  
  #Include plant
  #      Empty$Exclude_2022[target_row]<-"Include"
  
  #    } } #closes inner if statements
  
  ##Lastly, if plant died prior to gopher damage, list it as a natural death
  #  if(Empty$Date_first_death_observed_2022[target_row]<Empty$Date_first_gopher_2022[target_row]) {
  
  #    Empty$Exclude_2022[target_row]<-"Include"
  
  #  } #Closes inner if statement
  
  #} #Closes outer if statement
  
  
  
  
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  ifelse((sum(sub$State=="Reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Bolted_2022[target_row]<-1, Empty$Bolted_2022[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered_2022[target_row]<-1, Empty$Flowered_2022[target_row]<-0)
  
  ####Subset reproductive data
  Repro<-subset(sub, State=="Reproductive" & (Flower_number>=1 | Silique_number>=1)) #this filters out the rows where the plant was reproductive, but is missing flower or silique number
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering_2022[target_row]<-Repro$Ordinal_date[1] }
  
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  if(nrow(Repro)>=1) {
    Empty$Stem_number_flowering_2022[target_row]<-Repro$Stem_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Flower_number_flowering_2022[target_row]<-Repro$Flower_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_number_flowering_2022[target_row]<-Repro$Silique_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_length_flowering_2022[target_row]<-Repro$Silique_length[1] 
    #<-Repro$Silique_length[1] 
  }
  
  
  
  
  #Height1 through Height3 at first flower
  if(nrow(Repro)>=1) {
    Empty$Height1_flowering_2022[target_row]<-Repro$Height1[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height2_flowering_2022[target_row]<-Repro$Height2[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height3_flowering_2022[target_row]<-Repro$Height3[1] }
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
    
    #Date_peak_flowering		
    max_fl<-max(Repro$Flower_number, na.rm=T)
    
    Empty$Date_peak_flowering_2022[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
    
    
    #Flower_number_peak and Silique_number_peak and Stem_number_peak
    Empty$Flower_number_peak_2022[target_row]<-max_fl
    
    Empty$Silique_number_peak_2022[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Stem_number_peak_2022[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
    
    
    #Height1_peak through Height3_peak
    Empty$Height1_peak_2022[target_row]<-Repro$Height1[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height2_peak_2022[target_row]<-Repro$Height2[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height3_peak_2022[target_row]<-Repro$Height3[which(Repro$Flower_number==max_fl)][1]
    
    
  } #Closes if statement
  
  
  
  
  ########Data for date of first SO - Note that this is the date of first SO observed after the last flower was observed
  
  if(sum(Repro$Silique_number, na.rm=T)>=1) { #outer if statement
    
    #Extract the rowname for the last day a flower was observed, and then the rownumber of that rowname within sub (note that the rownames in sub will come from its name within the master, long-format data file)
    last_Fl<-rownames(tail(Repro[which(Repro$Flower_number>=1),], n=1))
    
    ##If there is a date of last flower, adjust the repro subset to reflect only data after the date of last flower
    if(length(last_Fl)>0) { #inner if statement
      
      ##Row number of date of last flower in sub
      last_Fl<-which(grepl(last_Fl, rownames(Repro)))
      
      ##Subset data to exclude censuses on or before the last flower was observed
      
      Repro<-Repro[-c(1:last_Fl),]
      
    } #Closes inner if statement
    
    
    #Subset for censuses where a silique was observed (we have to do this second subset to exlude plants that became b or b, tco after flowering)
    Repro<-subset(Repro, Silique_number>=1)
    
    
    ##Now estimate the date of first silique
    Empty$Date_silique_2022[target_row]<-Repro$Ordinal_date[1]
    
    
    #Stem_number_silique, Silique_number_silique, and Silique_length_silique
    Empty$Stem_number_silique_2022[target_row]<-Repro$Stem_number[1] 
    
    Empty$Silique_number_silique_2022[target_row]<-Repro$Silique_number[1]
    
    Empty$Silique_length_silique_2022[target_row]<-Repro$Silique_length[1]
    
    
    #Height3_silique
    Empty$Height1_silique_2022[target_row]<-Repro$Height1[1]
    
    Empty$Height2_silique_2022[target_row]<-Repro$Height2[1]
    
    Empty$Height3_silique_2022[target_row]<-Repro$Height3[1]
    
    
    
  } #Closing bracket from outer if statement
  
  
  
  
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number_2022[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number 
  Empty$Failed_silique_number_2022[target_row]<-sum(sub$Failed_silique_number, na.rm=T)
  
  ##Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  #if (sum(sub$Silique_number, na.rm=T)>=1) { 
  
  # Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
  
  ##If the number of failed siliques is negative, set to zero
  #if(Empty$Failed_silique_number[target_row]<0) {
  
  # Empty$Failed_silique_number[target_row]<-0
  
  #Mature_length_siliques - note that the column numbers have to be updated
  Empty$Mature_length_siliques_2022[target_row]<-sum(sub[,c(40,
                                                            42,
                                                            44,
                                                            46,
                                                            48,
                                                            50,
                                                            52,
                                                            54,
                                                            56)], na.rm=T)
  
} #closes inner if statement



check1 <-subset(Empty, PlantID=="G16375")
check2 <-subset(Empty, PlantID=="G16316")

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "gothic2022_2022_summary.xlsx", sheetName="Summary2022", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 

###### schofield 2022 cohort ######


##Data
scho<-read.csv("2022_sinam_census22.csv", header=T)
sapply(scho,class)

head(scho,20)
tail(scho)

scho$Master_notes<-as.character(scho$Master_notes)
scho$Flower_number <as.integer(scho$Flower_number)

reproductive<-subset(scho, State=="Reproductive")
head(reproductive,20)
tail(reproductive,20)




##Empty data set that we will populate with summary data


Empty_input<-read.csv("2022_sinam_empty22.csv", header=T)
Empty <- Empty_input[order(Empty_input$PlantID),]
row.names(Empty) <- 1:nrow(Empty)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(scho $PlantID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(scho $PlantID)[i]
  sub<-subset(scho, PlantID==target_ID)
  
  #        ##If the subset has less than 21 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==3)==F)
  { print(sub$PlantID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(scho$PlantID))){ 
  
  ##Subset the data for plantID i
  target_ID<-unique(scho$PlantID)[i]
  sub<-subset(scho, PlantID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$PlantID==target_ID),]))
  
  #overwinter_survival - Census 1 was the first full census, alive=1 dead=0
  ifelse((sub$Status[1]=="Alive"), Empty$Overwinter_survival_2022[target_row]<-1, Empty$Overwinter_survival_2022[target_row]<-0)
  
  
  #season_survival - Census 11 was the last full census, alive=1 dead=0
  ifelse((sub$Status[3]=="Alive"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  
  #season_Include - Census 1 has the plants that died due to trasnplant shock , alive=1 dead=0
  #ifelse((sub$Include[10]=="Alive"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  Empty$Exclude_2022[target_row]<-sub$Include[3] 
  
  #error here
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census
  
  if(Empty$Season_survival_2022[target_row]==1) {
    Empty$Overwinter_survival_2022[target_row]<-1
  }
  
  
  # leaf damage was assessed on  censuses 1,2
  
  Empty$LAR_1_2022[target_row]<-sub$LAR[1] 
  Empty$LAR_2_2022[target_row]<-sub$LAR[2] 

  
  ## Master notes
  Empty$Master_notes_2022[target_row]<-sub$Master_notes[3] 
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  if(Empty$Season_survival_2022[target_row]==0) { #Outer if statement
    
    #Subset data where Status==Dead	
    First_death<-subset(sub, Status=="Dead")
    
    #Extract row name of last census where status==Alive
    Last_alive<-rownames(tail(sub[which(sub$Status=="Alive"),], n=1))
    
    ##Adjust First_death data set to include only data 	after the last Alive status was observed
    if(length(Last_alive)>=1) { #Inner if statement
      
      ##Adjust subset of First_death
      First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
      
    } #Closes inner if statement
    
    ##Assign the date of first death
    Empty$Date_first_death_observed_2022[target_row]<-First_death$Ordinal_date[1] 
    
  } #Closes outer if statement
  
  ##Gopher damage or Censor - was coded as GopherDamage or Censor in the Exclude column - sometimes a plant marked as being exposed to gopher damage doesn't die - for that reason, this has to meet conditions that it was censused as gopher damaged at least once and that we didn't observe it as alive in the last census
  #if(any(sub$Exclude[10]=="GopherDamage") & (sub$Status[10]!="Alive")) {
  #  Empty$Exclude_2022[target_row]<-"GopherDamage"
  #}
  
  #if(any(sub$Exclude[10]=="Censor") & (sub$Status[10]!="Alive")) {
  #  Empty$Exclude_2022[target_row]<-"Censor"
  #}
  
  #if(is.na(Empty$Exclude_2022[target_row])==T) {
  #  Empty$Exclude_2022[target_row]<-"Include"
  #}
  
  
  
  ##Date_first_gopher, again only noted if the plant was not marked as alive in the last census
  ##This code is a little crazy because we are also calculating the date a plant was last alive prior to gopher damage and whether the living plant died prior to gopher damage
  #if(any(sub$Exclude[10]=="GopherDamage") & (sub$Status[10]!="Alive" )) {
  
  #  Empty$Date_first_gopher_2022[target_row]<-sub$Ordinal_date[which(sub$Exclude=="GopherDamage")][1] 
  
  ##Last date of being alive prior to gopher damage
  #  Last_alive<-tail(sub$Ordinal_date[which(sub$Status=="Alive" & sub$Ordinal_date<Empty$Date_first_gopher_2022[target_row])], n=1)
  
  #Last date of being dead prior to gopher damage
  #  Last_dead<-tail(sub$Ordinal_date[which(sub$Status=="Dead" & sub$Ordinal_date<Empty$Date_first_gopher_2022[target_row])], n=1)	
  
  #If the plant has both a Last_alive and Last_dead date
  #  if(length(Last_alive)>=1 & length(Last_dead)>=1) {
  
  #If the plant was alive before it was dead and if it died before gopher was first observed - have to separate this if statement from the last because you can't ask if something is less than something else if never occured
  #    if(Last_alive<Last_dead & Last_dead<Empty$Date_first_gopher_2022[target_row]) {
  
  #Include plant
  #      Empty$Exclude_2022[target_row]<-"Include"
  
  #    } } #closes inner if statements
  
  ##Lastly, if plant died prior to gopher damage, list it as a natural death
  #  if(Empty$Date_first_death_observed_2022[target_row]<Empty$Date_first_gopher_2022[target_row]) {
  
  #    Empty$Exclude_2022[target_row]<-"Include"
  
  #  } #Closes inner if statement
  
  #} #Closes outer if statement
  
  
  
  
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  ifelse((sum(sub$State=="Reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Bolted_2022[target_row]<-1, Empty$Bolted_2022[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered_2022[target_row]<-1, Empty$Flowered_2022[target_row]<-0)
  
  ####Subset reproductive data
  Repro<-subset(sub, State=="Reproductive" & (Flower_number>=1 | Silique_number>=1)) #this filters out the rows where the plant was reproductive, but is missing flower or silique number
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering_2022[target_row]<-Repro$Ordinal_date[1] }
  
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  if(nrow(Repro)>=1) {
    Empty$Stem_number_flowering_2022[target_row]<-Repro$Stem_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Flower_number_flowering_2022[target_row]<-Repro$Flower_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_number_flowering_2022[target_row]<-Repro$Silique_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_length_flowering_2022[target_row]<-Repro$Silique_length[1] 
    #<-Repro$Silique_length[1] 
  }
  
  
  
  
  #Height1 through Height3 at first flower
  if(nrow(Repro)>=1) {
    Empty$Height1_flowering_2022[target_row]<-Repro$Height1[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height2_flowering_2022[target_row]<-Repro$Height2[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height3_flowering_2022[target_row]<-Repro$Height3[1] }
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
    
    #Date_peak_flowering		
    max_fl<-max(Repro$Flower_number, na.rm=T)
    
    Empty$Date_peak_flowering_2022[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
    
    
    #Flower_number_peak and Silique_number_peak and Stem_number_peak
    Empty$Flower_number_peak_2022[target_row]<-max_fl
    
    Empty$Silique_number_peak_2022[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Stem_number_peak_2022[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
    
    
    #Height1_peak through Height3_peak
    Empty$Height1_peak_2022[target_row]<-Repro$Height1[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height2_peak_2022[target_row]<-Repro$Height2[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height3_peak_2022[target_row]<-Repro$Height3[which(Repro$Flower_number==max_fl)][1]
    
    
  } #Closes if statement
  
  
  
  
  ########Data for date of first SO - Note that this is the date of first SO observed after the last flower was observed
  
  if(sum(Repro$Silique_number, na.rm=T)>=1) { #outer if statement
    
    #Extract the rowname for the last day a flower was observed, and then the rownumber of that rowname within sub (note that the rownames in sub will come from its name within the master, long-format data file)
    last_Fl<-rownames(tail(Repro[which(Repro$Flower_number>=1),], n=1))
    
    ##If there is a date of last flower, adjust the repro subset to reflect only data after the date of last flower
    if(length(last_Fl)>0) { #inner if statement
      
      ##Row number of date of last flower in sub
      last_Fl<-which(grepl(last_Fl, rownames(Repro)))
      
      ##Subset data to exclude censuses on or before the last flower was observed
      
      Repro<-Repro[-c(1:last_Fl),]
      
    } #Closes inner if statement
    
    
    #Subset for censuses where a silique was observed (we have to do this second subset to exlude plants that became b or b, tco after flowering)
    Repro<-subset(Repro, Silique_number>=1)
    
    
    ##Now estimate the date of first silique
    Empty$Date_silique_2022[target_row]<-Repro$Ordinal_date[1]
    
    
    #Stem_number_silique, Silique_number_silique, and Silique_length_silique
    Empty$Stem_number_silique_2022[target_row]<-Repro$Stem_number[1] 
    
    Empty$Silique_number_silique_2022[target_row]<-Repro$Silique_number[1]
    
    Empty$Silique_length_silique_2022[target_row]<-Repro$Silique_length[1]
    
    
    #Height3_silique
    Empty$Height1_silique_2022[target_row]<-Repro$Height1[1]
    
    Empty$Height2_silique_2022[target_row]<-Repro$Height2[1]
    
    Empty$Height3_silique_2022[target_row]<-Repro$Height3[1]
    
    
    
  } #Closing bracket from outer if statement
  
  
  
  
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number_2022[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number 
  Empty$Failed_silique_number_2022[target_row]<-sum(sub$Failed_silique_number, na.rm=T)
  
  ##Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  #if (sum(sub$Silique_number, na.rm=T)>=1) { 
  
  # Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
  
  ##If the number of failed siliques is negative, set to zero
  #if(Empty$Failed_silique_number[target_row]<0) {
  
  # Empty$Failed_silique_number[target_row]<-0
  
  #Mature_length_siliques - note that the column numbers have to be updated
  Empty$Mature_length_siliques_2022[target_row]<-sum(sub[,c(44,46,48)], na.rm=T)
  
} #closes inner if statement



check1 <-subset(Empty, PlantID=="G16375")
check2 <-subset(Empty, PlantID=="G16316")

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "schofield2022_2022_summary.xlsx", sheetName="Summary2022", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 

