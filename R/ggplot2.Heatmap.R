create <-function ( ExpressionSet=NULL, data=NULL, samples=NULL, annotation=NULL ){
	obj <- list(sampleGroups = list(), geneGroups = list()  )
	if ( !is.null(ExpressionSet) ) {
		obj$exprSet <- ExpressionSet
	}
	else {
		err= 'Expression set is not defined'
		if ( is.null(data) ) {
			err <- paste( err, ' and data is missing')
		} 
		if ( is.null(samples) ) {
			err <- paste( err, ' and samples is missing')
		}
		if ( is.null(annotation) ) {
			err <- paste( err, ' and annotation is missing')
		}
		if ( nchar(err) > 30 ) {
			stop( err)
		}
		obj$exprSet <- list( data = data, samples=samples, annotation =annotation )
	}
	class(obj) <- 'ggplot2.Heatmap'
	obj
}

## the first SampleGroup does not necessarily need a color as it can be covered by the facet_grid blocks
set.SampleGroup <- function( x, groupName='GroupName', color=NULL ) {
	UseMethod('add.SampleGroup', x)
}
set.SampleGroup.ggplot2.Heatmap <- function (x, groupName='GroupName', color=NULL ) {
	if ( length( grep ( groupName, names( x$sampleGroups) )) == 0 ){
		if ( length(grep( groupName,colnames(x$exprSet$samples))) == 0 ){
			stop ( paste("Sorry I do not know the sample table column ", groupName )) 
		}
		l <- length(unique(x$exprSet$samples[,groupName]))
		if ( is.null(color) ) {
			color = rainbow( l )
		}
		if ( length(color) != l ) {
			color = rainbow( l )
		}
		id <- length( names( x$sampleGroups) ) +1
		x$sampleGroups[[id]] <- color
		names(x$sampleGroups)[id] <- groupName
	}
	x
}

set.GeneGroup <- function( x, groupName='GroupName', color=NULL ) {
	UseMethod('add.GeneGroup', x)
}
set.GeneGroup.ggplot2.Heatmap <- function (x, groupName='GroupName', color=NULL ) {
	if ( length( grep ( groupName, names( x$geneGroups) )) == 0 ){
		if ( length(grep( groupName,colnames(x$exprSet$annotation))) == 0 ){
			stop ( paste("Sorry I do not know the annotation table column ", groupName )) 
		}
		l <- length(unique(x$exprSet$annotation[,groupName]))
		if ( is.null(color) ) {
			color = rainbow( l )
		}
		if ( length(color) != l ) {
			color = rainbow( l )
		}
		id <- length( names( x$geneGroups) ) +1
		x$geneGroups[[id]] <- color
		names(x$geneGroups)[id] <- groupName
	}
	x
}


melt4heatmap <- function( x, sampleGroup='GroupName', sampleColGroup='GroupName', sampleColor=NULL, 
		geneNames="Gene.Symbol", geneGroup=NULL, geneColGroup=NULL, geneColor=NULL ){
	UseMethod('melt4heatmap', x)
}

melt4heatmap.ggplot2.Heatmap <- function (x, sampleGroups=c('GroupName'), 
		sampleColor=NULL, geneNames="Gene.Symbol", geneGroup=NULL, geneColGroup=NULL, geneColor=NULL ){
	
	melted <- melt( x, groupcol='GroupName', colCol='GroupName', probeNames="Gene.Symbol" )
	if ( is.null(sampleColor) ){
		sampleColor <- rainbow( length(unique(melted$ColorGroup)))
	}
	
	melted_new <- NULL
	melted_new <- as.matrix(melted[1:nrow(obj$data),])
	x <- as.vector( t(melted[nrow(obj$data),]))
	melted_new <- rbind( melted_new, matrix(c('SampleGroup', as.character(x[2]), as.character(x[5]), as.character(x[4]), as.character(x[5]), as.character(x[6]) ), nrow=1) )
	
	for ( i in 2:ncol(obj$data)) {
		melted_new <- rbind( melted_new, as.matrix(melted[(1+nrow(obj$data)*(i-1)):(nrow(obj$data)*i),]) )
		x <- as.vector( t(melted[nrow(obj$data)*i,]))
		melted_new <- rbind( melted_new, matrix(c('SampleGroup', as.character(x[2]), as.character(x[5]), as.character(x[4]), as.character(x[5]), as.character(x[6]) ), nrow=1) )
	}
	
}

addSampleColGroup <- function ( x, melted, name='GroupName', color=NULL ) {
	if ( is.null(color) ){
		color <- rainbow( length(unique(x$samples[,name])))
	}
	## have I added any other sample groups before this group?
	
	
}