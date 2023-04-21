

# ipack function: install and load multiple R packages. 
# check to see if packages are installed. Install them if they are not, then load them into the R session. 

ipack <- function(pkg){ 
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])] 
  if (length(new.pkg))  
    install.packages(new.pkg, dependencies = TRUE) 
  sapply(pkg, require, character.only = TRUE) 
  
}


packages=c("data.table",  # function :: rbindlist
           "plyr",        # function :: ddply/join
           "ggplot2",
           "rootSolve",       #multiroot
           "ismev",
           "Lmoments",
           "numDeriv",
           "homtest",
           "lmomco",
           "lmom",
           "evd",
           "psych",
           "dplyr",
           "plyr",
           "foreach",
           "writexl",
           "zipfR",
           "actuaryr",
           "Rsolnp",
           'latex2exp')       


ipack(packages)