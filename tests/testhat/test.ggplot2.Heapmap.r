source('../../load_libs.R') ## this has to be changed to the real next line once documentation is finished
#library('ggplot2.Heatmap')
source ('../../R/ggplot2.Heatmap.R')
obj <- list (data = matrix ( sample( 0:12, 10000, replace=T), nrow=100 ) )
print ( dim(obj$data))
obj$samples <- data.frame( SampleName = paste( rep('Sample', 100), 1:100, sep='_') )
obj$samples$GroupName<- paste( 'Group',kmeans( t(obj$data), centers=10)$cluster)
obj$samples$GroupName<- factor( obj$samples$GroupName, levels= paste( 'Group',1:10))
obj$samples$RGroup <- factor(rep( LETTERS[1:10], 10 ), levels=LETTERS[1:10] )
colnames(obj$data) <- obj$samples$SampleName
rownames(obj$data) <-  paste('Gene', 1:100,sep='')
obj$annotation <- data.frame( 
		ProbeName = rownames(obj$data), 
		clusters = factor( paste( 'GGroup',kmeans( obj$data, centers=10)$cluster ), labels=paste('GGroup', 1:10) ),
		RFgroup = factor( paste( 'RFGroup',kmeans( obj$data, centers=14)$cluster ), labels=paste('RFGroup', 1:14) )
)

t <- create( ExpressionSet= obj )
t <- set.SampleGroup ( t, 'GroupName', color=rainbow(length(unique(obj$samples$GroupName))))
te <- create( ExpressionSet= obj )
te <- set.SampleGroup ( te, 'GroupName' )

all.equal(te$sampleGroups, t$sampleGroups )

te <- set.SampleGroup ( te, 'GroupName' )
all.equal(te$sampleGroups, t$sampleGroups )


t <- set.GeneGroup ( t, 'clusters', color=rainbow(length(unique(obj$samples$GroupName))))
te <- set.GeneGroup ( te, 'clusters')

all.equal(t$geneGroups, te$geneGroups)

t <- melt(t, probeNames='ProbeName')
all.equal(dim(t$melted),c(10000,5) )
## now I need a second SampleGroup option
t <- set.SampleGroup ( t, 'RGroup' )
t <- melt(t, probeNames='ProbeName')
all.equal(dim(t$melted),c(10100,5) )

## now lets add another gene group
t <- set.GeneGroup ( t, 'RFgroup', color=grey.colors(14) )
t <- melt(t, probeNames='ProbeName')
all.equal(dim(t$melted),c(10200,5) )

plot(t)
## great working!!!!



