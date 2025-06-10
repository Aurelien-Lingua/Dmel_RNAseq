# Download Drosophila genome 
# Code edited from 
# "Download Saccharomyces Cerevisiae genome
# copyright (c) 2022 - Danny Arends"
# to enable download of all genome fasta files for latest drosophila melanogaster release
# Modifier - Aurelien Lingua
#
# Last modified - 6.10.2025
#   Changes: changed uri and base strings to 
#   allow seamless download from most current accession 
# Note: as Kallisto is used for pseudoalignment, using primary assembly or soft-masked (or toplevel)
# doesn't matter much since the mapping is done to the transcriptome. Change to "_sm" version of 
# accession if using whole genome alignment tools (STAR, HiSAT2).

# Install and load the RCurl package
install.packages("RCurl")
library(RCurl)

# Define the FTP URL and credentials
ftp_url <- "https://ftp.ensembl.org/pub/current/fasta/drosophila_melanogaster/dna/"
# userpwd <- "demo:password"  # Username:Password
# Directory to download into
dir <- paste0(getwd(),"/genome/")
# Get the directory listing from the FTP server
ftp_listing <- getURL(ftp_url, ftp.use.epsv = FALSE, .opts = "dirlistonly")

# Split the listing into individual file names
files <- strsplit(ftp_listing, "\r*\n")[[1]]

# extracting only the fasta files filenames that are unmasked
fasta_names <- regmatches(files, 
                          regexpr("Drosophila_melanogaster.*?\\.dna\\..*?\\.fa\\.gz", files))
cat("Files available on FTP server:\n", fasta_names, "\n")

# remove toplevel file
fasta_names<- fasta_names[-grep("toplevel",fasta_names)]

# Loop through the list of files and download each one
for (file_name in fasta_names) {
  print(file_name)
  # Construct full FTP URL for each file
  file_url <- paste0(ftp_url, file_name)

  # Define the local path where the file will be saved
  local_file <- paste0(dir,file_name)

  # Download the file from the FTP server
  binary_data <- getBinaryURL(file_url)
  writeBin(binary_data, local_file)

  cat("Downloaded:", file_name, "\n")
}

cat("All files downloaded successfully.")

chrs <- c("nonchromosomal", "primary_assembly.2L","primary_assembly.2R","primary_assembly.3L","primary_assembly.3R","primary_assembly.4","primary_assembly.X","primary_assembly.Y",
          "primary_assembly.mitochondrion_genome")


# Create an empty file for the whole assembly
cat("", file = "genome/Drosophila_melanogaster.BDGP6.54.dna.primary_assembly.fa")
# loopover files and add them to whole assembly
for(file_name in fasta_names){
  fname <- paste0("genome/", file_name)
  # Extract and merge into a fast file
  cmd <- paste0("wsl zcat ", fname, " >> genome/Drosophila_melanogaster.BDGP6.54.dna.primary_assembly.fa")
  cat(cmd, "\n")
  system(cmd)
}

# Compress the fasta file using bgzip (keep original with -k) --> ned to do in cmdline atm
cmd <- paste0("wsl bgzip genome/Drosophila_melanogaster.BDGP6.54.dna.primary_assembly.fa -k")
cat(cmd, "\n")

