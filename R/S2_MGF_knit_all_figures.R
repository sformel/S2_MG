#S2_MGF

#Knit all figures into separate html folder

#last updated 11/11/2021
#by Steve Formel

FL <- list.files(path = "R", pattern = ".Rmd")

#trim extension
FL <- str_remove(FL, ".Rmd")

#knit
lapply(FL, function(x){
  rmarkdown::render(paste0('R/', x, '.Rmd'), 
  output_file = paste0('html_output/', x, '.html'))
})
