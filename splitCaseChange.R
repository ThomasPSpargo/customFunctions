#---
#Date: 19/10/2022
#Author: Thomas Spargo <thomas.spargo@kcl.ac.uk>
#---

#Quick function for splitting a string when case changes.
#Original purpose is to facilitate identification of how protein residues are distributed across exons of a gene.

## INPUT
# x     = A character string to be split at points where adjacent characters change between upper/lower case.
# delim = Defaults to " ", indicates the delimiter to use when splitting the string;
#         included as an argument in case the default " " delimiter is already used within the string and should not be used for splitting.

## OUTPUT
# List:
# $x_elements = The split string sequence, returned as a vector.
# $pos_range  = List of start and end positions for each $x_elements. These are returned cumulatively with respect to the previous element of $x_elements.
#               (e.g. If x is an amino acid sequence to be split by exon, then this list is start and end positions of residues associated with each exon.

splitCaseChange <- function(x,delim=" "){

#gsub to identify string separation points at case change, insert delimiter.
string<- gsub('([[:upper:]])([[:lower:]])', paste0('\\1',delim,'\\2'), x)
string<- gsub('([[:lower:]])([[:upper:]])', paste0('\\1',delim,'\\2'), string)

#Split the string by delimiter
splitstring <- strsplit(string,split=delim)[[1]]

#Identify length of each element
ct_length <- nchar(splitstring)

#Create list with generic names to store cumulative start and endpoints of each element
st_end <- vector(mode="list",length(ct_length))
names(st_end) <- paste0("Element_",1:length(st_end))
st_end[[1]] <- c(1,ct_length[1]) #store first element

#Cumulatively sum the remaining elements, and add to list
for(i in 2:length(ct_length)){
  ct_length[i]<- ct_length[i]+ct_length[i-1]
  
  st_end[[i]] <- c(ct_length[i-1]+1,ct_length[i])
}

#Return split string elements and cumulative positions in a list
return(list(x_elements=splitstring,pos_range=st_end))
}

#Functionality test with SOD1 AA sequence, obtained from ENSEMBL, with exons alternating by case
#https://useast.ensembl.org/Homo_sapiens/Transcript/Sequence_Protein?db=core;g=ENSG00000142168;r=21:31659666-31668931;t=ENST00000270142


# seq <- "MATKAVCVLKGDGPVQGIINFEQKesngpvkvwgsikglteglhgfhvhefgdntaGCTSAGPHFNPLSRKHGGPKDEERhvgdlgnvtadkdgvadvsiedsvislsgdhciigrtlvVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ"
# 
# seq
# splitCaseChange(seq)
