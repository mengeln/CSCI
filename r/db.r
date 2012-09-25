library(RSQLite)
conn <- dbConnect("SQLite", "data/bug_metadata.db")
dbWriteTable(conn, "ibi", ibiv4)
dbWriteTable(conn, "taxonomy", taxonomy_v5)
dbWriteTable(conn, "maxmin", maxmin)
dbWriteTable(conn, "otu_crosswalk", otu_crosswalk)
dbWriteTable(conn, "example_pred", pred)
dbWriteTable(conn, "example_bugs", bugs)
dbWriteTable(conn, "predcal", predcal)
grps <- data.frame(grps.final)
dbWriteTable(conn, "grps_final", grps)

test1 <- dbGetQuery(conn, "SELECT * FROM grps_final")$grps_final

forestsdb = local({load("data/FinalForests.Rdata"); environment()})
tools:::makeLazyLoadDB(forestsdb, "forestsdb")
