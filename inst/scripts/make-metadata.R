
meta <- data.frame(
    Title = "Mapping table of metabolite IDs and common names from different sources",
    Description = paste0("Large mapping table including metabolite common names ",
    		         "and 9 distinct metabolite ID formats that orignated ",
			 "from 4 different databases: HMDB, Comptox Dashboard, ",
			 "ChEBi, and graphite R package."),
    BiocVersion = "3.12",
    Genome = NA,
    SourceType = "TSV",
    SourceUrl = "https://github.com/yigbt/metaboliteIDmapping/blob/master/inst/scripts/make-data.R",
    SourceVersion = "Jun 15 2020",
    Species = NA,
    TaxonomyId = NA,
    Coordinate_1_based = NA,
    DataProvider = "HMDB, EMBL-EBI, EPA",
    Maintainer = "Sebastian Canzler <sebastian.canzler@ufz.de>",
    RDataClass = "Tibble",
    DispatchClass = "Rda",
    RDataPath = "metaboliteIDmapping/v2/metabolitesMapping.rda",
    Tags = "metabolites:mapping:HMDB:KEGG:ChEBI:Pubchem:Comptox:CAS:Drugbank"
)

write.csv(meta, file="metadata.csv", row.names=FALSE)
