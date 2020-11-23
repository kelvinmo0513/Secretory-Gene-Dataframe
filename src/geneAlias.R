### To install package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
                     
### Reference and tutorials
# https://urldefense.com/v3/__https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf__;!!Mih3wA!Snp5TxGWavvZYJFUl2rpDQhZTD2SCpSsFLqvqSgFUuCUCh73EHFbYZfeOVMu$ 


##Load packages
library(org.Hs.eg.db)
library(readxl)

## Upload curations (note you will need to upload the rest, either manually or put all in a folder and loop through the files in the folder)
C1 <- read_excel('/Users/Kelvi/Desktop/Curation_Angel.xlsx')
C2 <- read_excel('/Users/Kelvi/Desktop/Curation_AM.xlsx')
C3 <- read_excel('/Users/Kelvi/Desktop/Curation_Zerong.xlsx')
C4 <- read_excel('/Users/Kelvi/Desktop/Curation_Moji.xlsx')
C5 <- read_excel('/Users/Kelvi/Desktop/Curation_HM.xlsx')
C6 <- read_excel('/Users/Kelvi/Desktop/Curation_Nat.xlsx')
C7 <- read_excel('/Users/Kelvi/Desktop/Curation_Isaacs.xlsx')
C8 <- read_excel('/Users/Kelvi/Desktop/Curation_HB.xlsx')
C9 <- read_excel('/Users/Kelvi/Desktop/Curation_Curtis.xlsx')
C10 <- read_excel('/Users/Kelvi/Desktop/Curation_Caressa.xlsx')

## Extract all the unique genes (again you will need to include all the others, either manually or with a loop)
uGenes <- unique(c(C1$`Generic Gene ID`, C2$`Generic Gene ID`, C3$`Generic Gene ID`, C4$`Generic Gene ID`, C5$`Generic Gene ID`, C6$`Generic Gene ID`, C7$`Generic Gene ID`, C8$`Generic Gene ID`, C9$`Generic Gene ID`, C10$`Generic Gene ID`))

### ____________ Modify gene list to ensure they are either Alias or Symbol keys in org db object

## Some useful functions to explore the db object
columns(org.Hs.eg.db) #see what kinds of data rea retrievable via select
keytypes(org.Hs.eg.db) # list what kinds of fields can be used as keys to query database (similar to columns)
help("ALIAS") #use help to know more about them
keys(org.Hs.eg.db, keytype="SYMBOL") #use 'keys' extract keys of a particular type, in this case 'SYMBOL'

## Find which genes are symbol keys
symKEY <- intersect(uGenes, keys(org.Hs.eg.db, keytype="SYMBOL"))
length(symKEY) 
# note: not all genes are found in the symbol keys
noSYM <- setdiff(uGenes, keys(org.Hs.eg.db, keytype = "SYMBOL"))
length(noSYM)

## For remaining genes, which are alias keys
aliasKEY <- intersect(noSYM, keys(org.Hs.eg.db, keytype="ALIAS"))
length(aliasKEY) #note: still some genes in neither alias nor symbol keys

## For remaining keys, investigate manually
setdiff(uGenes, c(symKEY, aliasKEY)) # to display all the genes that were not able to map

# For example: ST6Gal1 and ST6Gal2 need to be capitalized and then you can find them in SYMBOL, edit uGenes list
intersect(c('ST6GAL1', 'ST6GAL2'), keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<ST6Gal1\\>","ST6GAL1",uGenes)
uGenes_mod <- gsub("\\<ST6Gal2\\>","ST6GAL2",uGenes_mod)

# Another example: ALG02 is found in SYMBOL as ALG2 (probably similary for ALG09)
intersect('ALG2', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<ALG02\\>","ALG2",uGenes_mod)

# Kelvin: ALG09 is found in SYMBOL as ALG9
intersect('ALG9', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<ALG09\\>","ALG9",uGenes_mod)

# Kelvin: GALNT19 is found in SYMBOL as GALNT9
intersect('GALNT9', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<GALNT19\\>","GALNT9",uGenes_mod)

# Kelvin: MGAT5 mannosyl -glycoprotein beta-1,6-N-acetyl-glucosaminyltransferase is found in SYMBOL as MGAT5
intersect('MGAT5', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<MGAT5 mannosyl -glycoprotein beta-1,6-N-acetyl-glucosaminyltransferase\\>","MGAT5",uGenes_mod)

# Kelvin: Rer1p is in SYMBOL as RER1
intersect('RER1', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<Rer1p\\>","RER1",uGenes_mod)

# Kelvin: CK2 is found in SYMBOL as CSNK2A1
intersect('CSNK2A1', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<CK2\\>","CSNK2A1",uGenes_mod)

# Kelvin: HSP70 is found in SYMBOL as HSPA1A
intersect('HSPA1A', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<HSP70\\>","HSPA1A",uGenes_mod)

# Kelvin: HSP90 is found in SYMBOL as HSP90AA1
intersect('HSP90AA1', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<HSP90\\>","HSP90AA1",uGenes_mod)

# Kelvin: UBQLN5 is found in SYMBOL as UBQLNL
intersect('UBQLNL', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<UBQLN5\\>","UBQLNL",uGenes_mod)

# Another example: Google 'UGGT gene' and open genecards (one of my fav gene websites).  UGGT seems to be an alias for UGGT1
intersect('UGGT1', keys(org.Hs.eg.db, keytype='SYMBOL')) #check to see if UGGT1 is actually in the db object, yes it is
uGenes_mod <- gsub("\\<UGGT\\>","UGGT1",uGenes_mod)

# Another example: DPM 1 needs to have the space removed
intersect('DPM1', keys(org.Hs.eg.db, keytype='SYMBOL'))
uGenes_mod <- gsub("\\<DPM 1\\>","DPM1",uGenes_mod)

## For weird stuff like NA and 43530.0, look at curations to see if those rows actually contain any information
#NA
C1[is.na(C1$`Generic Gene ID`),]
C2[is.na(C2$`Generic Gene ID`),]
C3[is.na(C3$`Generic Gene ID`),]
C4[is.na(C4$`Generic Gene ID`),]
C5[is.na(C5$`Generic Gene ID`),]
C6[is.na(C6$`Generic Gene ID`),]
C7[is.na(C7$`Generic Gene ID`),]
C8[is.na(C8$`Generic Gene ID`),]
C9[is.na(C9$`Generic Gene ID`),]
C10[is.na(C10$`Generic Gene ID`),]
# You can see that NA does not contain any information so can remove from gene list
uGenes_mod <- uGenes_mod[!is.na(uGenes_mod)]
#43530.0
C1[C1$`Generic Gene ID`=='43530.0',]
C2[C2$`Generic Gene ID`=='43530.0',]
C3[C3$`Generic Gene ID`=='43530.0',]
C4[C4$`Generic Gene ID`=='43530.0',]
C5[C5$`Generic Gene ID`=='43530.0',]
C6[C6$`Generic Gene ID`=='43530.0',]
C7[C7$`Generic Gene ID`=='43530.0',]
C8[C8$`Generic Gene ID`=='43530.0',]
C9[C9$`Generic Gene ID`=='43530.0',]
C10[C10$`Generic Gene ID`=='43530.0',]
#This one is tricky, you can see upon investiation of Curation_Helen that excel has turned gene MARCH6 into a date, which then got even more screwed up upon loading into R
intersect('MARCH6', keys(org.Hs.eg.db, keytype='ALIAS')) #tried symbol first, wasn't there, tried alias next... MARCH6 is there
uGenes_mod <- gsub("\\<43530.0\\>","MARCH6",uGenes_mod)

#43530
C1[C1$`Generic Gene ID`=='43530',]
C2[C2$`Generic Gene ID`=='43530',]
C3[C3$`Generic Gene ID`=='43530',]
C4[C4$`Generic Gene ID`=='43530',]
C5[C5$`Generic Gene ID`=='43530',]
C6[C6$`Generic Gene ID`=='43530',]
C7[C7$`Generic Gene ID`=='43530',]
C8[C8$`Generic Gene ID`=='43530',]
C9[C9$`Generic Gene ID`=='43530',]
C10[C10$`Generic Gene ID`=='43530',]
uGenes_mod <- gsub("\\<43530\\>","MARCH6",uGenes_mod)

#43530
C1[C1$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C2[C2$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C3[C3$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C4[C4$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C5[C5$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C6[C6$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C7[C7$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C8[C8$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C9[C9$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]
C10[C10$`Generic Gene ID`=='HMGCR, NPLOC4, UFD1, IGF1R, RNF19A',]

##We have to be careful when we manually change gene names.  Things to take into account: 
  #1) we may introduce dupicates so at the end we have to take the unique values'
  #2) we have to be sure to manyally include the origical gene name found in the curation as an alias so it can be mapped back later during the actual mering process
uGenes_mod <- unique(uGenes_mod)

### ____________ Retrieve gene alias for each gene in list
## Use either the select method or mapIds method to extract alternate geneIDs (refer to link above)
## At this point all of the genes in uGenes_mod should belong to either the symbol or alias keys
symKEY = intersect(uGenes_mod, keys(org.Hs.eg.db, keytype="SYMBOL"))
aliasKEY = intersect(setdiff(uGenes_mod, symKEY), keys(org.Hs.eg.db, keytype='ALIAS'))
setdiff(uGenes_mod, c(symKEY, aliasKEY)) #This should theoretically be empty

#dfGenes <- data.frame(matrix("", ncol = 2, nrow = length(c(symKEY, aliasKEY))))
#names(dfGenes) <- c("GenericGeneSymbol", "Alias")
#dfGenes$GenericGeneSymbol[1:length(symKEY)] <- symKEY
#dfGenes$Alias[(length(symKEY)+1):nrow(dfGenes)] <- aliasKEY

## Find aliases for gene symbols
gAlias <- select(org.Hs.eg.db, keys=symKEY, columns=c("ALIAS"), keytype="SYMBOL")
gAlias_mod <- aggregate(.~SYMBOL , gAlias, paste, collapse = ",")

## Find gene sybols for aliases
gSymbol <- select(org.Hs.eg.db, keys=aliasKEY, columns=c('SYMBOL'), keytype='ALIAS') #this maps the alias to the original gene symbol
gSymbol2 <- select(org.Hs.eg.db, keys=gSymbol$SYMBOL, columns=('ALIAS'), keytype='SYMBOL') #extract ALL aliases for the gene symbol
gSymbol_mod <- aggregate(.~SYMBOL , gSymbol2, paste, collapse = ",")

## Combine
dfGenes <- rbind(gSymbol_mod, gAlias_mod)

## Manually add aliases that you used in the previous manal conversions
# For example: ST6Gal1 and ST6Gal2
dfGenes[dfGenes$SYMBOL=='ST6GAL1',2] <- paste(dfGenes[dfGenes$SYMBOL=='ST6GAL1',2], 'ST6Gal1', sep=',')
dfGenes[dfGenes$SYMBOL=='ST6GAL2',2] <- paste(dfGenes[dfGenes$SYMBOL=='ST6GAL2',2], 'ST6Gal2', sep=',')
# Another example: ALG02
dfGenes[dfGenes$SYMBOL=='ALG2',2] <- paste(dfGenes[dfGenes$SYMBOL=='ALG2',2], 'ALG02', sep=',')
# ... keep going to add aliases for all the genes you had to manipulate manyally before
dfGenes[dfGenes$SYMBOL=='ALG9',2] <- paste(dfGenes[dfGenes$SYMBOL=='ALG9',2], 'ALG09', sep=',')
dfGenes[dfGenes$SYMBOL=='DPM1',2] <- paste(dfGenes[dfGenes$SYMBOL=='DPM1',2], 'DPM 1', sep=',')
dfGenes[dfGenes$SYMBOL=='GALNT9',2] <- paste(dfGenes[dfGenes$SYMBOL=='GALNT9',2], 'GALNT19', sep=',')
dfGenes[dfGenes$SYMBOL=='MGAT5',2] <- paste(dfGenes[dfGenes$SYMBOL=='MGAT5',2], 'MGAT5 mannosyl -glycoprotein beta-1,6-N-acetyl-glucosaminyltransferase', sep=',')
dfGenes[dfGenes$SYMBOL=='RER1',2] <- paste(dfGenes[dfGenes$SYMBOL=='RER1',2], 'Rer1p', sep=',')
dfGenes[dfGenes$SYMBOL=='CSNK2A1',2] <- paste(dfGenes[dfGenes$SYMBOL=='CSNK2A1',2], 'CK2', sep=',')
dfGenes[dfGenes$SYMBOL=='HSPA1A',2] <- paste(dfGenes[dfGenes$SYMBOL=='HSPA1A',2], 'HSP70', sep=',')
dfGenes[dfGenes$SYMBOL=='HSP90AA1',2] <- paste(dfGenes[dfGenes$SYMBOL=='HSP90AA1',2], 'HSP90', sep=',')
dfGenes[dfGenes$SYMBOL=='UBQLNL',2] <- paste(dfGenes[dfGenes$SYMBOL=='UBQLNL',2], 'UBQLN5', sep=',')
dfGenes[dfGenes$SYMBOL=='UGGT1',2] <- paste(dfGenes[dfGenes$SYMBOL=='UGGT1',2], 'UGGT', sep=',')

# Manually add in APF4 gene since it is not int the database
LastSet <- c('APF4', 'APF4')
dfGenes <- rbind(dfGenes, LastSet)

# Export the database into excel file
install.packages("writexl")
library("writexl")
write_xlsx(dfGenes, "/Users/Kelvi/Desktop/dfGenes.xlsx")

#######
# The end result should be a two column data frame that contains gene SYMBOLs and their corresponding ALIASes
# Each of the genes from the curation template should map to either a gene name in SYMBOL, or a gene name in ALIAS
# The rows in this dataframe are the final rows that you'll use for the merge

