sch29 <- read.csv("http://faculty.missouri.edu/huangf/data/MLM/diagmath.csv")
str(sch29)
sn = names(table(sch29$sch_id))
d2 <- data.frame(sch_id = sn, schid = 1:29)
d2
sch29 <- merge(d2, sch29, by = 'sch_id')
sch29$sch_id <- NULL
sch29 <- dplyr::arrange(sch29, schid)
table(sch29$byhomewk)

str(sch29)
getwd()
save(sch29, file = 'data/sch29.rda')
