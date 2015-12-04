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


## brks has to be the breaks in wich to cut the expression values will be
## default = c( -20.1, quantile(melted$Expression[which(melted$Expression != -20)],seq(0,1,by=0.1)) )
plot.ggplot2.Heatmap <- function( x, brks=NULL, col=NULL ) {
	x <- melt(x, brks=brks)
	## built the color vector
	## + finalize the Expression by creating the required factor levels
	l <- x$datalevels ## created during the melt call
	if ( is.null(col) | length(col) != length(l) ){
		col = c( 'gray', bluered(length(x$datalevels) -1) )
	}
	n <- names(x$sampleGroups)
	if ( length(n) > 1 ){
		for ( id in 2:length(names(x$sampleGroups))) {
			col <- c( col, x$sampleGroups[[id]] )
			if ( is.factor(x$exprSet$samples[,n[id] ])){
				l <- c(l, levels(x$exprSet$samples[,n[id] ]) )
			}else{
				l <- c(l, unique(x$exprSet$samples[,n[id] ]) )
			}
			
		}
	}
	n <- names(x$geneGroups)
	if ( length(n) > 1 ){
		for ( id in 2:length(names(x$geneGroups))) {
			col <- c( col, x$geneGroups[[id]] )
			if ( is.factor(x$exprSet$annotation[,n[id] ])){
				l <- c(l, levels(x$exprSet$annotation[,n[id] ]) )
			}else{
				l <- c(l, unique(x$exprSet$annotation[,n[id] ]) )
			}
		}
	}
	
	x$melted$Expression <- factor( as.vector(x$melted$Expression), levels= l )
	p = ( ggplot(x$melted, aes(x=SampleName,y=ProbeName))
				+ geom_tile(aes(fill=Expression)) 
				+ scale_fill_manual( values = col ) 
				+ theme(
						axis.text.x=element_blank(),
#axis.ticks.x=element_line(color=ss$colrss),
						axis.ticks.length=unit(0.00,"cm")
				)+ labs( y='') )
	if ( ncol(x$melted) == 5 ){
		p <- p + facet_grid( GeneGroup ~ Group,scales="free", space='free')
	}else if ( ncol(x$melted) == 4 ) {
		p <- p + facet_grid( . ~ Group,scales="free", space='free')
	}
	p
}

## the first SampleGroup does not necessarily need a color as it can be covered by the facet_grid blocks
set.SampleGroup <- function( x, groupName='GroupName', color=NULL ) {
	UseMethod('set.SampleGroup', x)
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
	UseMethod('set.GeneGroup', x)
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


## brks has to be the breaks in wich to cut the expression values will be
## default = c( -20.1, quantile(melted$Expression[which(melted$Expression != -20)],seq(0,1,by=0.1)) )
melt.ggplot2.Heatmap <- function (x, probeNames=NULL, brks=NULL){
#	if ( ! exists ('melted', where =x) ) {
		ma  <- x$exprSet$data
		if ( ! is.null(probeNames) ) {
			rownames(ma) <- forceAbsoluteUniqueSample(x$exprSet$annotation[, probeNames] )
		}
		x$melted <- melt( cbind(rownames(ma),ma) )
		
		if ( length( which ( x$melted[,2] == '') ) > 0 ){
			x$melted <- x$melted[ - which ( x$melted[,2] == ''),]
		}
		## bin the expression values
		if ( is.null(brks)){
			x$melted[,3] <- as.numeric(as.vector(x$melted[,3]))
			brks= c( -20.1, quantile(x$melted[which(x$melted[,3] != -20),3],seq(0,1,by=0.1)) )
			brks[length(brks)] = brks[length(brks)] + 0.1
		}
		x$melted[,3] <- cut( x$melted[,3], breaks= brks)
		x$datalevels <- levels(x$melted[,3])
		
		##now we need to add the sample and gene colors
		x <- addSampleColGroup ( x )
		x <- addGeneColGroup( x )
		colnames(x$melted) <- c('ProbeName', 'SampleName', 'Expression', 'Group', 'GeneGroup')[1:ncol(x$melted)]
#	}
	x
}

addSampleColGroup <- function ( x ) {
	if ( length(names(x$sampleGroups)) > 0 ) {
		n <-  names(x$sampleGroups)
		x$x_facet <- n[1]
		grps <- NULL
		datarows = nrow(x$exprSet$data)
		for ( i in as.vector(obj$samples[, n[1]]) ){
			grps <- c( grps, rep( i, datarows ))
		}
		if ( is.factor(x$exprSet$samples[,x$x_facet])){
			x$melted$Group <- factor( grps, levels=levels(x$exprSet$samples[,x$x_facet]))
		}else{
			x$melted$Group <- factor( grps, levels=unique(x$exprSet$samples[,x$x_facet]))
		}
		if ( length(n) > 1 ){
			## color groups!!
			for ( GNid in 2:length(n)){
				le <- datarows + GNid -2
				melted_new <- NULL
				for (sid in 1:datarows) {
					melted_new <- rbind( melted_new, as.matrix(x$melted[(1+le*(sid-1)):(le*sid),]) )
					line <- as.vector( t(x$melted[le*sid,]))
					melted_new <- rbind(melted_new,  matrix(c('SampleGroup', as.character(line[2]), as.character(x$exprSet$samples[sid,n[GNid]]), as.character(line[4]) ), nrow=1) )
				}
				x$melted <- data.frame(melted_new,row.names= 1:nrow(melted_new))
			}
		}
	}
	else {
		if ( exisits('x_facet',where=x) ) {
			rm(x$x_facet)
		}
	}
	x
}

addGeneColGroup <- function ( x ) {
	if ( length(names(x$geneGroups)) > 0 ) {
		n <-  names(x$geneGroups)
		x$y_facet <- n[1]
		grps <- NULL
		datarows = nrow(x$exprSet$data)
		sampleNULL <- c()
		if ( length(names(x$sampleGroups)) > 1 ) {
			sampleNULL <- rep (NA, length(names(x$sampleGroups)) -1 )
		}
		for ( i in as.vector(obj$annotation[, n[1]]) ){
			grps <- c( grps, as.vector(obj$annotation[, n[1]]), sampleNULL)
		}
		x$melted$GeneGroup <- grps
		if ( length(n) > 1 ){
			## color groups!! here we need to add a matrix with the colors after all other data blocks
			## keep, keep, <Gene group name>, 'Gene Group', NA
			for ( GNid in 2:length(n)){
				add  <- as.matrix(x$melted[1:nrow(x$exprSet$data),])
				add[,3] <- as.vector(x$exprSet$annotation[,n[GNid]])
				add[,4] <-'GeneGroup'
				x$melted <- rbind( x$melted, add )
			}
		}
	}
	else {
		if ( exisits('y_facet',where=x) ) {
			rm(x$y_facet)
		}
	}
	x
}

## coped from ExpressionSet to not depend on the ExpressionSet lib (which is not working at the moment)
forceAbsoluteUniqueSample <- function ( x ,separator='_') {
	last = ''
	ret <- vector(length=length(x))
	for ( i in 1:length(x) ){
		if ( is.null(ret) ){
			last = x[i]
			ret[i] <- last
		}
		else{
			last = x[i]
			if ( ! is.na(match( last, ret )) ){
				last <- paste(last,separator,sum( ! is.na(match( x[1:i], last )))-1, sep = '')
			}
			ret[i] <- last
		}
	}
	ret
}
