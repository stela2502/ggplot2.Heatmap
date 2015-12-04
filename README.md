# ggplot2.Heatmap
A R package to create heatmaps using ggplot2 trying to implement the functionality of the gplots Heatmap.3 function.
With main focus on the color groups. Dendrograms are not supported in any case and are not on the to do list.

This package will not work as standard R package at the moment.
You instead have to include the R script file R/ggplot2.Heatmap.R

source('R/ggplot2.Heatmap.R')
source('load_libs.R')
obj <- create( data= <data frame with e.g. expression data genes in rows, samples in columns>, samples=<data frame describing the samples same order as in data>, annotation=<data frame describing the genes same order data rows> )

Afterwards you are able to define row or column level annotations to plot.
These annotations have to be included in the sample or annotation tables.

 obj <- set.SampleGroup ( obj, 'GroupName' ) ## where GroupName is a column in the samples table
 obj <-  set.GeneGroup ( obj, 'clusters', color= gey.colors( length(unique(as.vector(obj$exprSet$annotation$cluster))) ## where clusters is a column in the annotation table.
 
If the variables in the table are factors, the order in the factors should be used to order the groups (untested at the moment)
If no color argument is given a rainbow scale will be created.

A final step is to
 
plot(obj) 


Which will create a ggplot2 plot of the data.

Hope you can use the tool!

