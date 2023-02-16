###A function to summarise coefficients produced by glmmTMB into a latex xtable

library(tools) ### Contains function to title case

glmmTMB.latex <- function(glmm.model) {
 
glmm.table  <- as.data.frame(round(summary(glmm.model)$coefficients$cond, digits=3))
glmm.table  <- glmm.table[order(abs(glmm.table$Estimate), decreasing = TRUE), ]
rownames(glmm.table)  <-  gsub("scale", "", rownames(glmm.table))
rownames(glmm.table)  <-  gsub("\\(", "", rownames(glmm.table))
rownames(glmm.table)  <-  gsub("\\)", "", rownames(glmm.table))
rownames(glmm.table)  <-  gsub("dem.", "", rownames(glmm.table))
rownames(glmm.table)  <-  toTitleCase(rownames(glmm.table))
rownames(glmm.table)  <-gsub("^(x)", "Easting (m)", rownames(glmm.table)) ###first character in string only
rownames(glmm.table)  <-gsub("^(y)", "Northing (m)", rownames(glmm.table)) 
rownames(glmm.table)  <-  gsub(":", " x ", rownames(glmm.table))
names(glmm.table)  <- c("Estimate", "Estimate Std Error", "Z value", "P value")
glmm.table[glmm.table[ , "P value"] == 0, "P value" ]  <- "0.001"     

##return(xtable(glmm.table))
return(glmm.table)

}
