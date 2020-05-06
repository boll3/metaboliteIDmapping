
meta <- data.frame(
    Title = "Metabolite ID mapping from different sources",
    Description = paste0("Large mapping table including 9 distinct metabolite ",
                         "ID formats that orignated from 4 different databases: ",
                         "HMDB, Comptox Dashboard, ChEBi, and graphite R package."),
    BiocVersion = "3.12",
    Genome = NA,
    SourceType = NA,
    SourceUrl = "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip,
                ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv,
                https://comptox.epa.gov/dashboard/downloads",
    SourceVersion = "Apr 30 2020",
    Species = NA,
    TaxonomyId = NA,
    Coordinate_1_based = NA,
    DataProvider = "HMDB, EMBL-EBI, EPA",
    Maintainer = "Sebastian Canzler <sebastian.canzler@ufz.de>",
    RDataClass = "Tibble",
    DispatchClass = "Rda",
    RDataPath = "metaboliteIDmapping/v1/metabolitesMapping.rda",
    Tags = "metabolites:mapping:HMDB:KEGG:ChEBI:Pubchem:Comptox:CAS:Drugbank"
)

write.csv(meta, file="metadata.csv", row.names=FALSE)
