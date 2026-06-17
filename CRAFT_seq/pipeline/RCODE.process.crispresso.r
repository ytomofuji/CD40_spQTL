## ---------------------------
##
## Script name: RCODE.process.crispresso.r
##
## Purpose of script: Process the CRISPRESSO output from a plate into a single file (supposed to be used in the wrapper script [craft_seq_process_gDNA.sh])
##
##
## Date Created: 2025-05-06
##
## Email:
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------



source("misc_funcs.R")
library(stringr)
library(tidyverse)

argv=commandArgs(T)

file_path=argv[1]
name=argv[2]

if(FALSE){
file_path="/path/to/gDNA/CRISPRESSO/gDNA1"
name="Plate1"
}

AllelesFunction <- 
function (filepath, number_of_alleles = 10) #Filepath is the folder with all the allele tables. 
{
    filepath <- sub("/$", "", filepath)
    files <- list.files(path = filepath, pattern = "Alleles_frequency_table.txt$", full.names = T, 
        recursive = T) 
        
    matrix <- rbindlist(lapply(files, function(f) {
    # Extract Plate_Well ID from the file path
    well_id <- str_extract(f, "(?<=CRISPResso_on_)[^/]+")
    
        # Read file and add Plate_Well column
        mat_file <- fread(f, nrows = number_of_alleles, select = c(1, 5, 6, 7, 8, 9))
        mat_file[, Well_ID := well_id]  # Assign the extracted Plate_Well ID

        return(mat_file)  # Return the modified table
    }), use.names = TRUE, fill = TRUE)  # Bind all tables together

    matrix = cbind(matrix[, (ncol(matrix)):(ncol(matrix))], matrix[, 1:(ncol(matrix)-1)])
    
    return(matrix)
}

#need to load conda R env
alleles_DNA <- AllelesFunction(file_path)
alleles_DNA <- alleles_DNA %>% mutate(DNA_Barcode=name)%>%as_tibble()

write_tsv(alleles_DNA, str_c(name,".alleles_DNA.tsv"))

