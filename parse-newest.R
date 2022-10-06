library(readxl)
library(stringr)

exclude.ol = T

# read data from Excel sheets and covert to dataframes
# between-generation data
gen.tib = read_excel("Data/ALL_mother-progeny_2022-08-11.xlsx")
gen.df = data.frame(gen.tib)
gen.df = gen.df[!is.na(gen.df$progenyID),]
# within-child data
tissue.tib = read_excel("Data/ALL_within plant_2022-08-11.xlsx")
tissue.df = data.frame(tissue.tib)
# decide whether to exclude OL data
if(exclude.ol == T) {
  tissue.df = tissue.df[tissue.df$tissue_time != "OL",]
}
tissue.df = tissue.df[!is.na(tissue.df$individual),]

gen.df$sibling = sapply(str_split(gen.df$progenyID, "[^B-C0-9]"), tail, 1)
tissue.df$sibling = tissue.df$individual

# loop through organelle-background sets
organelles = c("mito", "plastid")
backgrounds = c("MSH1", "WILD")

for(organelle in organelles) {
  for(background in backgrounds) {
    message(paste(c("Parsing ", organelle, " ", background), sep=""))
    if(organelle == "plastid" & background == "WILD") { break }
    
   # subset out this case
  gen.sub = gen.df[gen.df$organelle==organelle & gen.df$background==background & !is.na(gen.df$time),]
  tissue.sub = tissue.df[tissue.df$organelle==organelle & tissue.df$background==background & !is.na(tissue.df$time),]
  
  # assign individual IDs
  gen.sub$individual = gen.sub$sibling
  
  if(nrow(tissue.sub) > 0) {
    tissue.sub$devtime = 0
    for(i in 1:nrow(tissue.sub)) {
      tissue.sub$devtime[i] = ifelse(tissue.sub$tissue_time[i] == "OL" & tissue.sub$time.point[i] == 2, 1, tissue.sub$time.point[i])
    }
  }
  
  # pull useful data into dataframes
  gen.final = data.frame(family = tolower(gen.sub$unique.familiy.ID), individual = gen.sub$individual, time = gen.sub$time, h0 = gen.sub$mother_.altSNV, h = gen.sub$individual_.altSNV, src = gen.sub$notes)
  tissue.final = data.frame(family = tolower(tissue.sub$unique.familiy.ID), individual = tissue.sub$individual, time = tissue.sub$time.point, h0 = tissue.sub$Mom_.altSNV, h = tissue.sub$Tissue_.alt, src = tissue.sub$Notes)
  
  # pull together and remove initial mother observations
  all.final = rbind(gen.final, tissue.final)
  all.final = all.final[all.final$time != 0,]
  
  # process these observations separately
  allfamilies = unique(all.final$family)
  for(i in 1:length(allfamilies)) {
    all.final = rbind(all.final, data.frame(family = allfamilies[i], individual = 0, time = 0, h0 = -1, h = all.final$h0[which(all.final$family == allfamilies[i])[1]], src = all.final$src[which(all.final$family == allfamilies[i])[1]]))
  }
  all.final = all.final[order(all.final$family),]
  
  # assign a unique string reference
  all.final$uniqueID = ""
  all.final$uniqueref = all.final$familyref = 0
  for(i in 1:nrow(all.final)) {
    all.final$uniqueID[i] = paste(c(all.final$family[i], "-", all.final$individual[i]), collapse="")
  }
  allIDs = unique(all.final$uniqueID)
  
  for(i in 1:nrow(all.final)) {
    all.final$uniqueref[i] = which(allIDs == all.final$uniqueID[i])
    all.final$familyref[i] = which(allfamilies == all.final$family[i])
  }
  
  # order by reference and bound h measurements
  all.final = all.final[order(all.final$uniqueref),]
  all.final$h[all.final$h > 100] = 100
  all.final$h[all.final$h < 0] = 0
  
  # write to file
  filename = paste(c("newest-data-all-", organelle, "-", background, ".csv"), collapse="")
  write.csv(data.frame(id = all.final$uniqueref, family = all.final$familyref, time = all.final$time, h = all.final$h), filename, row.names = FALSE)
  
  old.final = all.final[all.final$src == "pub",]
  old.final$uniqueref = old.final$uniqueref - min(old.final$uniqueref) + 1
  old.final$familyref = old.final$familyref - min(old.final$familyref) + 1
  filename = paste(c("newest-data-old-", organelle, "-", background, ".csv"), collapse="")
  write.csv(data.frame(id = old.final$uniqueref, family = old.final$familyref, time = old.final$time, h = old.final$h), filename, row.names = FALSE)
  
  new.final = all.final[all.final$src == "new",]
  new.final$uniqueref = new.final$uniqueref - min(new.final$uniqueref) + 1
  new.final$familyref = new.final$familyref - min(new.final$familyref) + 1
  filename = paste(c("newest-data-new-", organelle, "-", background, ".csv"), collapse="")
  write.csv(data.frame(id = new.final$uniqueref, family = new.final$familyref, time = new.final$time, h = new.final$h), filename, row.names = FALSE)
  }
}
  
  