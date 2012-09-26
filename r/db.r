library(RSQLite)
conn <- dbConnect("SQLite", "data/bug_metadata.db")
dbWriteTable(conn, "ibi", ibiv4)
dbWriteTable(conn, "taxonomy", taxonomy_v5)
dbWriteTable(conn, "maxmin", maxmin)
dbWriteTable(conn, "otu_crosswalk", otu_crosswalk)
dbWriteTable(conn, "example_pred", pred)
dbWriteTable(conn, "example_bugs", bugs)
dbWriteTable(conn, "predcal", predcal)
dbWriteTable(conn, "bugcal_pa", bugcal.pa)
grps <- data.frame(grps.final)
dbWriteTable(conn, "grps_final", grps)

test1 <- dbGetQuery(conn, "SELECT * FROM grps_final")$grps_final

forestsdb <- local({load("data/FinalForests.Rdata"); rf.mod <- rf.mod; environment()})
tools:::makeLazyLoadDB(forestsdb, "forestsdb")


dbGetQuery(conn, paste("Select * FROM ibi WHERE FinalID is 'Trichoptera'"))

dbGetQuery(conn, paste("UPDATE ibi SET SAFIT1 = 'Missing' WHERE SAFIT1 is null"))
dbGetQuery(conn, paste("UPDATE ibi SET SAFIT1 = FinalID WHERE SAFIT1 == 'Missing'"))