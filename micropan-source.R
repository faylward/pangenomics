#' @name binomixEstimate
#' @aliases binomixEstimate
#' 
#' @title Binomial mixture model estimates
#' 
#' @description Fits binomial mixture models to the data given as a pan-matrix. From the fitted models
#' both estimates of pan-genome size and core-genome size are available.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' @param K.range The range of model complexities to explore. The vector of integers specify the number
#' of binomial densities to combine in the mixture models.
#' @param core.detect.prob The detection probability of core genes. This should almost always be 1.0,
#' since a core gene is by definition always present in all genomes, but can be set fractionally smaller.
#' @param verbose Logical indicating if textual output should be given to monitor the progress of the
#' computations.
#' 
#' @details  A binomial mixture model can be used to describe the distribution of gene clusters across
#' genomes in a pan-genome. The idea and the details of the computations are given in Hogg et al (2007),
#' Snipen et al (2009) and Snipen & Ussery (2012).
#' 
#' Central to the concept is the idea that every gene has a detection probability, i.e. a probability of
#' being present in a genome. Genes who are always present in all genomes are called core genes, and these
#' should have a detection probability of 1.0. Other genes are only present in a subset of the genomes, and
#' these have smaller detection probabilities. Some genes are only present in one single genome, denoted
#' ORFan genes, and an unknown number of genes have yet to be observed. If the number of genomes investigated
#' is large these latter must have a very small detection probability. 
#' 
#' A binomial mixture model with \samp{K} components estimates \samp{K} detection probabilities from the
#' data. The more components you choose, the better you can fit the (present) data, at the cost of less
#' precision in the estimates due to less degrees of freedom. \code{\link{binomixEstimate}} allows you to
#' fit several models, and the input \samp{K.range} specifies which values of \samp{K} to try out. There no
#' real point using \samp{K} less than 3, and the default is \samp{K.range=3:5}. In general, the more genomes
#' you have the larger you can choose \samp{K} without overfitting.  Computations will be slower for larger
#' values of \samp{K}. In order to choose the optimal value for \samp{K}, \code{\link{binomixEstimate}}
#' computes the BIC-criterion, see below.
#' 
#' As the number of genomes grow, we tend to observe an increasing number of gene clusters. Once a
#' \samp{K}-component binomial mixture has been fitted, we can estimate the number of gene clusters not yet
#' observed, and thereby the pan-genome size. Also, as the number of genomes grows we tend to observe fewer
#' core genes. The fitted binomial mixture model also gives an estimate of the final number of core gene
#' clusters, i.e. those still left after having observed \sQuote{infinite} many genomes.
#' 
#' The detection probability of core genes should be 1.0, but can at times be set fractionally smaller.
#' This means you accept that even core genes are not always detected in every genome, e.g. they may be
#' there, but your gene prediction has missed them. Notice that setting the \samp{core.detect.prob} to less
#' than 1.0 may affect the core gene size estimate dramatically.
#' 
#' @return \code{\link{binomixEstimate}} returns a \code{Binomix} object, which is a small (S3) extension
#' of a \code{list} with two components. These two components are named \samp{BIC.table} and \samp{Mix.list}.
#' 
#' The \samp{BIC.table} is a matrix listing, in each row, the results for each number of components used,
#' given by the input \samp{K.range}. The column \samp{Core.size} is the estimated number of core gene families,
#' the column \samp{Pan.size} is the estimated pan-genome size. The column \samp{BIC} is the Bayesian
#' Information Criterion (Schwarz, 1978) that should be used to choose the optimal value for \samp{K}.
#' The number of components where \samp{BIC} is minimized is the optimal. If minimum \samp{BIC} is reached
#' for the largest \samp{K} value you should extend the \samp{K.range} and re-fit. The function will issue
#' a \code{warning} to remind you of this.
#' 
#' The \samp{Mix.list} is a list with one element for each number of components tested. The content of each
#' \samp{Mix.list} element is a matrix describing one particular fitted binomial mixture model. A fitted model
#' is characterized by two vectors (rows) denoted \samp{Detect.prob} and \samp{Mixing.prop}. \samp{Detect.prob}
#' are the estimated detection probabilities, sorted in ascending order. The \samp{Mixing.prop} are the
#' corresponding mixing proportions. A mixing proportion is the proportion of the gene clusters having the
#' corresponding detection probability.
#' 
#' The generic functions \code{\link{plot.Binomix}} and \code{\link{summary.Binomix}}
#' are available for \code{Binomix} objects.
#' 
#' @references
#' Hogg, J.S., Hu, F.Z, Janto, B., Boissy, R., Hayes, J., Keefe, R., Post, J.C., Ehrlich, G.D. (2007).
#' Characterization and modeling of the Haemophilus influenzae core- and supra-genomes based on the
#' complete genomic sequences of Rd and 12 clinical nontypeable strains. Genome Biology, 8:R103.
#' 
#' Snipen, L., Almoy, T., Ussery, D.W. (2009). Microbial comparative pan-genomics using binomial
#' mixture models. BMC Genomics, 10:385.
#' 
#' Snipen, L., Ussery, D.W. (2012). A domain sequence approach to pangenomics: Applications to
#' Escherichia coli. F1000 Research, 1:19.
#' 
#' Schwarz, G. (1978). Estimating the Dimension of a Model. The Annals of Statistics, 6(2):461-464.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{chao}}, \code{\link{plot.Binomix}},
#' \code{\link{summary.Binomix}}.
#' 
#' @examples
#' # Loading a Panmat object in the micropan package
#' data(list="Mpneumoniae.blast.panmat",package="micropan")
#' 
#' # Estimating binomial mixture models
#' bino <- binomixEstimate(Mpneumoniae.blast.panmat,K.range=3:8)  # using 3,4,...,8 components
#' print(bino$BIC.table) # minimum BIC at 3 components
#' 
#' # Plotting the optimal model, and printing the summary
#' plot(bino)
#' summary(bino)
#' 
#' # Plotting the 8-component model as well
#' plot(bino,ncomp=8)  # clearly overfitted, we do not need this many sectors
#' 
#' # Plotting the distribution in a single genome
#' plot(bino,type="single")  # completely dominated by core genes
#' 
#' @export binomixEstimate
binomixEstimate <- function( pan.matrix, K.range=3:5, core.detect.prob=1.0, verbose=TRUE ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  y <- table( factor( colSums( pan.matrix ), levels=1:dim( pan.matrix )[1] ) )
  bic.tab <- matrix( NA, nrow=length(K.range), ncol=3 )
  colnames( bic.tab ) <- c( "Core.size", "Pan.size", "BIC" )
  rownames( bic.tab ) <- paste( K.range, "components" )
  mix.list <- vector( "list", length( K.range ) )
  for( i in 1:length( K.range ) ){
    if( verbose ) cat( "binomixEstimate: Fitting", K.range[i], "component model...\n" )
    lst <- binomixMachine( y, K.range[i], core.detect.prob )
    bic.tab[i,] <- lst[[1]]
    mix.list[[i]] <- lst[[2]]
  }
  if( bic.tab[length(K.range),3] == min( bic.tab[,3] ) ) warning( "Minimum BIC at maximum K, increase upper limit of K.range" )
  binomix <- list( BIC.table=bic.tab, Mix.list=mix.list )
  class( binomix ) <- c( "Binomix", "list" )
  return( binomix )
}

#' @importFrom stats constrOptim as.dendrogram dendrapply dist is.leaf optim prcomp sd
binomixMachine <- function( y, K, core.detect.prob=1.0 ){
  n <- sum( y )
  G <- length( y )
  ctr <- list( maxit=300, reltol=1e-6 )
  np <- K - 1
    
  pmix0 <- rep( 1, np )/K            # flat mixture proportions
  pdet0 <- (1:np)/(np+1)              # "all" possible detection probabilities
  p.initial <- c( pmix0, pdet0 )      # initial values for parameters
  # the inequality constraints...
  A <- rbind( c( rep( 1, np ), rep( 0, np ) ), c( rep( -1, np ), rep( 0, np ) ), diag( np+np ), -1*diag( np+np ) )
  b <- c( 0, -1, rep( 0, np+np ), rep( -1, np+np ) )
  
  # The estimation, maximizing the negative truncated log-likelihood function
  est <- constrOptim( theta=p.initial, f=negTruncLogLike, grad=NULL, method="Nelder-Mead", control=ctr, ui=A, ci=b, y=y, core.p=core.detect.prob )
  
  estimates <- numeric( 3 )
  names( estimates ) <- c( "Core.size", "Pan.size", "BIC" )
  estimates[3] <- 2*est$value + log(n)*(np+K)                         # the BIC-criterion
  p.mix <- c( 1 - sum( est$par[1:np] ), est$par[1:np] )               # the mixing proportions
  p.det <- c( core.detect.prob, est$par[(np+1):length( est$par )] )   # the detection probabilities
  ixx <- order( p.det )
  p.det <- p.det[ixx]
  p.mix <- p.mix[ixx]
    
  theta_0 <- choose( G, 0 ) * sum( p.mix * (1-p.det)^G )
  y_0 <- n * theta_0/(1-theta_0)
  estimates[2] <- n + round( y_0 )
  ixx <- which( p.det >= core.detect.prob )
  estimates[1] <- round( estimates[2] * sum( p.mix[ixx] ) )
    
  mixmod <- matrix( c( p.det, p.mix ), nrow=2, byrow=T )
  rownames( mixmod ) <- c( "Detection.prob", "Mixing.prop" )
  colnames( mixmod ) <- paste( "Comp_", 1:K, sep="" )

  return( list( estimates, mixmod ) )
}

negTruncLogLike <- function( p, y, core.p ){
  np <- length( p )/2
  p.det <- c( core.p, p[(np+1):length(p)] )
  p.mix <- c( 1-sum( p[1:np] ), p[1:np] )
  G <- length( y )
  K <- length( p.mix )
  n <- sum( y )
    
  theta_0 <- choose( G, 0 ) * sum( p.mix * (1-p.det)^G )
  L <- -n * log( 1 - theta_0 )
  for( g in 1:G ){
    theta_g <- choose( G, g ) * sum( p.mix * p.det^g * (1-p.det)^(G-g) )
    L <- L + y[g] * log( theta_g )
  }
  return( -L )
}

#' @rdname generic.Binomix
#' @name plot.Binomix
#' 
#' @title Plot and summary of \code{Binomix} objects
#' 
#' @description Generic functions for \code{Binomix} objects.
#' 
#' @param x A \code{Binomix} object, see below.
#' @param object A \code{Binomix} object, see below.
#' @param type Type of plot, default is \samp{type="pan"} which means the pie chart shows distribution
#' over the entire pan-genome. The alternative is \samp{type="single"} which means the pie chart will
#' show the distribution within a single (average) genome.
#' @param cex Plot symbol scaling.
#' @param ncomp Which model to display. You can override the display of the optimal (minimum BIC) model
#' by specifying the number of components here, e.g. \samp{ncomp=5} will always display the model with
#' \samp{5} components regardless of its BIC value.
#' @param show.bar Logical indicating if a colorbar should be displayed next to the pie.
#' @param \dots Optional graphical arguments.
#' 
#' @details A \code{Binomix} object contains a series of fitted binomial mixture models. It is a small
#' (S3) extension to a \code{list}, having two components. These are named \samp{BIC.table} and
#' \samp{Mix.list}, see \code{\link{binomixEstimate}} for more details.
#' 
#' The \code{\link{plot.Binomix}} function will display a \code{Binomix} object as a pie chart. Only
#' the model with the smallest BIC-criterion value is displayed. The BIC-criterion is used to rank the
#' various fitted models, and minimum BIC is an objective criterion for finding the best model complexity.
#' Each sector of the pie chart is a component, the color of the sector indicates its detection probability
#' and the size of the sector its mixing proportion. This pie chart illustrates how gene clusters are
#' distributed within the pan-genome. Sectors of (dark) blue color are highly conserved gene clusters
#' (core genes), sectors of greenish colors are medium conserved clusters (shell genes) and sectors of
#' orange/pink colors are non-conserved clusters (cloud genes).
#' 
#' The \code{\link{summary.Binomix}} function will print the estimated core size and pan-genome size for
#' the optimal component model.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{binomixEstimate}}.
#' 
#' @examples # See examples in the Help-file for binomixEstimate.
#' 
#' @importFrom graphics axis barplot box layout par pie plot points rect text
#' 
#' @export
plot.Binomix <- function( x, type="pan", cex=2, ncomp=NA, show.bar=TRUE, ... ){
  Binomix <- x
  if( is.na(ncomp) ){
    ncomp <- which( Binomix$BIC.table[,3] == min( Binomix$BIC.table[,3] ) )[1]
  } else {
    ncomp <- which( as.numeric( gsub( " components", "", rownames( Binomix$BIC.table ) ) ) == ncomp )
    if( length( ncomp ) != 1 ) stop( "Specified value of ncomp has not been fitted for this model" )
  }
  cpar <- par()$mar
  par( mar=c(2,2,2,2) )
  typ <- grep( type, c( "pan", "single" ) )
  fit <- Binomix$Mix.list[[ncomp]]
  dprob <- as.character( round( fit[1,]*1000 )/1000 )
  if( show.bar ){
    layout( matrix( c(1,1,1,1,1,1,2), nrow=1 ) )
    if( length( typ ) == 0 ){
      stop( "Unknown type specified" )
    } else if( typ == 1 ){
      pie( fit[2,], clockwise=T, col=panColor( fit[1,] ), labels=dprob, radius=1.0, cex=cex, ... )
    } else {
      eg <- fit[2,]*fit[1,]
      pie( eg/sum( eg ), clockwise=T, col=panColor( fit[1,] ), labels=dprob, radius=1.0, cex=cex, ... )
    }
    
    par( mar=c(1,1,1,3) )
    p <- (0:100)
    plot( rep( 0, 101 ), p, cex=0, col="white", xlim=c(0,1), ylim=c(0,100), 
          xaxt="n", yaxt="n", xlab="", ylab="" )
    box( lwd=3, col="white" )
    cols <- panColor( p/100 )
    xl <- rep( 0, 100 )
    xr <- rep( 1, 100 )
    yb <- (0:100)
    yt <- 1:101
    rect( xl, yb, xr, yt, col=cols, border=cols )
    axis( side=4, at=seq(0,100,10), labels=as.character( seq(0,100,10)/100) )
  } else {
    if( length( typ ) == 0 ){
      stop( "Unknown type specified" )
    } else if( typ == 1 ){
      pie( fit[2,], clockwise=T, col=panColor( fit[1,] ), labels=dprob, radius=1.0, cex=cex, ... )
    } else {
      eg <- fit[2,]*fit[1,]
      pie( eg/sum( eg ), clockwise=T, col=panColor( fit[1,] ), labels=dprob, radius=1.0, cex=cex, ... )
    }
  }
  par( mar=cpar )
}
#' @rdname generic.Binomix
#' @export
summary.Binomix <- function( object, ... ){
  ncomp <- which( object$BIC.table[,3] == min( object$BIC.table[,3] ) )[1]
  cat( "Minimum BIC model at", rownames( object$BIC.table )[ncomp], "\nFor this model:\n" )
  cat( "Estimated core size:", object$BIC.table[ncomp,1], "clusters\n" )
  cat( "Estimated pangenome size:", object$BIC.table[ncomp,2], "clusters\n" )
}

plotBar <- function(){
  cpar <- par()$mar
  par( mar=c(1,1,1,3) )
  p <- (0:100)
  plot( rep( 0, 101 ), p, cex=0, col="white", xlim=c(0,1), ylim=c(0,100), 
        xaxt="n", yaxt="n", xlab="", ylab="" )
  box( lwd=3, col="white" )
  cols <- panColor( p/100 )
  xl <- rep( 0, 100 )
  xr <- rep( 1, 100 )
  yb <- (0:100)
  yt <- 1:101
  rect( xl, yb, xr, yt, col=cols, border=cols )
  axis( side=4, at=seq(0,100,10), labels=as.character( seq(0,100,10)/100) )
}

#' @importFrom grDevices colorRampPalette
panColor <- function( p.vec ){
  level <- pretty( c(0,1), 100 )
  nlevel <- length( level )
  crp <- colorRampPalette( c("pink","orange","green","cyan","blue") )(nlevel)
  return( crp[1+round( 100*p.vec )] )
}



#' @name fluidity
#' @title Computing genomic fluidity for a pan-genome
#' 
#' @description Computes the genomic fluidity, which is a measure of population diversity.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' @param n.sim An integer specifying the number of random samples to use in the computations.
#' 
#' @details  The genomic fluidity between two genomes is defined as the number of unique gene
#' families divided by the total number of gene families (Kislyuk et al, 2011). This is averaged
#' over \samp{n.sim} random pairs of genomes to obtain a population estimate.
#' 
#' The genomic fluidity between two genomes describes their degree of overlap with respect to gene
#' cluster content. If the fluidity is 0.0, the two genomes contain identical gene clusters. If it
#' is 1.0 the two genomes are non-overlapping. The difference between a Jaccard distance (see
#' \code{\link{distJaccard}}) and genomic fluidity is small, they both measure overlap between
#' genomes, but fluidity is computed for the population by averaging over many pairs, while Jaccard
#' distances are computed for every pair. Note that only presence/absence of gene clusters are
#' considered, not multiple occurrences.
#' 
#' The input \samp{pan.matrix} is typically constructed by \code{\link{panMatrix}}.
#' 
#' @return A list with two elements, the mean fluidity and its sample standard deviation over
#' the \samp{n.sim} computed values.
#' 
#' @references Kislyuk, A.O., Haegeman, B., Bergman, N.H., Weitz, J.S. (2011). Genomic fluidity:
#' an integrative view of gene diversity within microbial populations. BMC Genomics, 12:32.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{distJaccard}}.
#' 
#' @examples 
#' # Loading two Panmat objects in the micropan package
#' data(list=c("Mpneumoniae.blast.panmat","Mpneumoniae.domain.panmat"),package="micropan")
#' 
#' # Fluidity based on a BLAST clustering Panmat object
#' fluid.blast <- fluidity(Mpneumoniae.blast.panmat)
#' 
#' # Fluidity based on domain sequence clustering Panmat object
#' fluid.domains <- fluidity(Mpneumoniae.domain.panmat)
#' 
#' @importFrom stats sd
#' 
#' @export fluidity
#' 
fluidity <- function( pan.matrix, n.sim=10 ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  flu <- rep( 0, n.sim )
  for( i in 1:n.sim ){
    ii <- sample( ng, 2 )
    g1 <- pan.matrix[ii[1],]
    g2 <- pan.matrix[ii[2],]
    flu[i] <- (sum( g1>0 & g2==0 )+sum( g1==0 & g2>0 ))/(sum(g1)+sum(g2))
  }
  flu.list <- list( Mean=mean( flu ), Std=sd( flu ) )
  return( flu.list )
}


#' @name distJaccard
#' @title Computing Jaccard distances between genomes
#' 
#' @description Computes the Jaccard distances between all pairs of genomes.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' 
#' @details The Jaccard index between two sets is defined as the size of the interesection of
#' the sets divided by the size of the union. The Jaccard distance is simply 1 minus the Jaccard index.
#' 
#' The Jaccard distance between two genomes describes their degree of overlap with respect to gene
#' cluster content. If the Jaccard distance is 0.0, the two genomes contain identical gene clusters.
#' If it is 1.0 the two genomes are non-overlapping. The difference between a genomic fluidity (see
#' \code{\link{fluidity}}) and a Jaccard distance is small, they both measure overlap between genomes,
#' but fluidity is computed for the population by averaging over many pairs, while Jaccard distances are
#' computed for every pair. Note that only presence/absence of gene clusters are considered, not multiple
#' occurrences.
#' 
#' The input \samp{pan.matrix} is typically constructed by \code{\link{panMatrix}}.
#' 
#' @return A \code{dist} object (see \code{\link{dist}}) containing all pairwise Jaccard distances
#' between genomes.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{fluidity}}, \code{\link{dist}}.
#' 
#' @examples
#' # Loading two Panmat objects in the micropan package
#' data(list=c("Mpneumoniae.blast.panmat","Mpneumoniae.domain.panmat"),package="micropan")
#' 
#' # Jaccard distances based on a BLAST clustering Panmat object
#' Jdist.blast <- distJaccard(Mpneumoniae.blast.panmat)
#' 
#' # Jaccard distances based on domain sequence clustering Panmat object
#' Jdist.domains <- distJaccard(Mpneumoniae.domain.panmat) 
#' 
#' @export distJaccard
#' 
distJaccard <- function( pan.matrix ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  dtab <- matrix( 0, nrow=ng, ncol=ng )
  rownames( dtab ) <- rownames( pan.matrix )
  colnames( dtab ) <- rownames( pan.matrix )
  for( i in 1:(ng-1)){
    g1 <- pan.matrix[i,]
    for( j in (i+1):ng ){
      cs <- g1+pan.matrix[j,]
      dtab[j,i] <- dtab[i,j] <- 1 - sum( cs>1 )/sum( cs>0 )
    }
  }
  return( as.dist( dtab ) )
}


#' @name distManhattan
#' @title Computing Manhattan distances between genomes
#' 
#' @description Computes the (weighted) Manhattan distances beween all pairs of genomes.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' @param scale An optional scale to control how copy numbers should affect the distances.
#' @param weights Vector of optional weights of gene clusters.
#' 
#' @details The Manhattan distance is defined as the sum of absolute elementwise differences between
#' two vectors. Each genome is represented as a vector (row) of integers in \samp{pan.matrix}. The
#' Manhattan distance between two genomes is the sum of absolute difference between these rows. If
#' two rows (genomes) of the \samp{pan.matrix} are identical, the corresponding Manhattan distance
#' is \samp{0.0}.
#' 
#' The \samp{scale} can be used to control how copy number differences play a role in the distances
#' computed. Usually we assume that going from 0 to 1 copy of a gene is the big change of the genome,
#' and going from 1 to 2 (or more) copies is less. Prior to computing the Manhattan distance, the
#' \samp{pan.matrix} is transformed according to the following affine mapping: If the original value in
#' \samp{pan.matrix} is \samp{x}, and \samp{x} is not 0, then the transformed value is \samp{1 + (x-1)*scale}.
#' Note that with \samp{scale=0.0} (default) this will result in 1 regardless of how large \samp{x} was.
#' In this case the Manhattan distance only distinguish between presence and absence of gene clusters.
#' If \samp{scale=1.0} the value \samp{x} is left untransformed. In this case the difference between 1
#' copy and 2 copies is just as big as between 1 copy and 0 copies. For any \samp{scale} between 0.0 and
#' 1.0 the transformed value is shrunk towards 1, but a certain effect of larger copy numbers is still
#' present. In this way you can decide if the distances between genomes should be affected, and to what
#' degree, by differences in copy numbers beyond 1. Notice that as long as \samp{scale=0.0} (and no
#' weighting) the Manhattan distance has a nice interpretation, namely the number of gene clusters that
#' differ in present/absent status between two genomes.
#' 
#' When summing the difference across gene clusters we can also up- or downweight some clusters compared
#' to others. The vector \samp{weights} must contain one value for each column in \samp{pan.matrix}. The
#' default is to use flat weights, i.e. all clusters count equal. See \code{\link{geneWeights}} for
#' alternative weighting strategies.
#' 
#' @return A \code{dist} object (see \code{\link{dist}}) containing all pairwise Manhattan distances
#' between genomes.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{distJaccard}}, \code{\link{geneWeights}},
#' \code{\link{panTree}}.
#' 
#' @examples 
#' # Loading two Panmat objects in the micropan package
#' data(list=c("Mpneumoniae.blast.panmat","Mpneumoniae.domain.panmat"),package="micropan")
#' 
#' # Manhattan distances based on a BLAST clustering Panmat object
#' Mdist.blast <- distManhattan(Mpneumoniae.blast.panmat)
#' 
#' # Manhattan distances based on domain sequence clustering Panmat object
#' Mdist.domains <- distManhattan(Mpneumoniae.domain.panmat,scale=0.5)
#' 
#' @importFrom stats dist
#' 
#' @export distManhattan
#' 
distManhattan <- function( pan.matrix, scale=0.0, weights=rep( 1, dim( pan.matrix )[2] ) ){
  if( (scale>1)|(scale<0) ){
    warning( "scale should be between 0.0 and 1.0, using scale=0.0" )
    scale <- 0.0
  }
  idx <- which( pan.matrix > 0, arr.ind=T )
  pan.matrix[idx] <- 1 + (pan.matrix[idx]-1)*scale
  
  pan.matrix <- pan.matrix * matrix( weights, nrow=dim( pan.matrix )[1], ncol=dim( pan.matrix )[2], byrow=T )
  dtab <- dist( pan.matrix, method="manhattan" )
  return( dtab )
}


#' @name geneWeights
#' @title Gene cluster weighting
#' 
#' @description This function computes weights for gene cluster according to their distribution in a pan-genome.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' @param type A text indicating the weighting strategy.
#' 
#' @details When computing distances between genomes or a PCA, it is possible to give weights to the
#' different gene clusters, emphasizing certain aspects.
#' 
#' As proposed by Snipen & Ussery (2010), we have implemented two types of weighting: The default
#' \samp{"shell"} type means gene families occuring frequently in the genomes, denoted shell-genes, are
#' given large weight (close to 1) while those occurring rarely are given small weight (close to 0).
#' The opposite is the \samp{"cloud"} type of weighting. Genes observed in a minority of the genomes are
#' referred to as cloud-genes. Presumeably, the \samp{"shell"} weighting will give distances/PCA reflecting
#' a more long-term evolution, since emphasis is put on genes who have just barely diverged away from the
#' core. The \samp{"cloud"} weighting emphasizes those gene clusters seen rarely. Genomes with similar
#' patterns among these genes may have common recent history. A \samp{"cloud"} weighting typically gives
#' a more erratic or \sQuote{noisy} picture than the \samp{"shell"} weighting.
#' 
#' @return A vector of weights, one for each column in \code{pan.matrix}.
#' 
#' @references Snipen, L., Ussery, D.W. (2010). Standard operating procedure for computing pangenome
#' trees. Standards in Genomic Sciences, 2:135-141.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{distManhattan}}.
#' 
#' @examples 
#' # Loading a Panmat object in the micropan package
#' data(list="Mpneumoniae.blast.panmat",package="micropan")
#' 
#' # Weighted Manhattan distances based on a BLAST clustering Panmat object
#' w <- geneWeights(Mpneumoniae.blast.panmat,type="shell")
#' Mdist.blast <- distManhattan(Mpneumoniae.blast.panmat,weights=w)
#' 
#' @export geneWeights
#' 
geneWeights <- function( pan.matrix, type=c("shell","cloud") ){
  ng <- dim( pan.matrix )[1]
  nf <- dim( pan.matrix )[2]
  pan.matrix[which( pan.matrix>0, arr.ind=T )] <- 1
  cs <- colSums( pan.matrix )
  
  midx <- grep( type[1], c( "shell", "cloud" ) )
  if( length( midx ) == 0 ){
    warning( "Unknown weighting:", type, ", using shell weights" )
    midx <- 1
  }
  W <- rep( 1, nf )
  x <- 1:ng
  ww <- 1/(1+exp( ((x-1)-(max(x)-1)/2)/((max(x)-1)/10) ))
  if( midx == 1 ) ww <- 1-ww
  for( i in 1:ng ) W[which( cs == i )] <- ww[i]
  return( W )
}
#' @name panMatrix
#' @title Computing the pan-matrix for a set of gene clusters
#' 
#' @description A pan-matrix has one row for each genome and one column for each gene cluster, and
#' cell \samp{[i,j]} indicates how many members genome \samp{i} has in gene family \samp{j}.
#' 
#' @param clustering A vector of integers indicating the gene cluster for every sequence. Sequences
#' with the same number belong to the same cluster. The name of each element is the tag identifying
#' the sequence.
#' 
#' @details The pan-matrix is a central data structure for pan-genomic analysis. It is a matrix with
#' one row for each genome in the study, and one column for each gene cluster. Cell \samp{[i,j]}
#' contains an integer indicating how many members genome \samp{i} has in cluster \samp{j}.
#' 
#' The input \code{clustering} must be an integer vector with one element for each sequence in the study,
#' typically produced by either \code{\link{bClust}} or \code{\link{dClust}}. The name of each element
#' is a text identifying every sequence. The value of each element indicates the cluster, i.e. those
#' sequences with identical values are in the same cluster. IMPORTANT: The name of each sequence must
#' contain the GID-tag for each genome, i.e. they must of the form \samp{GID111_seq1}, \samp{GID111_seq2},...
#' where the \samp{GIDxxx} part indicates which genome the sequence belongs to. See \code{\link{panPrep}}
#' for details.
#' 
#' The rows of the pan-matrix is named by the GID-tag for every genome. The columns are just named
#' \samp{Cluster_x} where \samp{x} is an integer copied from \samp{clustering}.
#' 
#' @return The returned object belongs to the class \code{Panmat}, which is a small (S3) extension to a
#' matrix. It can be treated as a matrix, but the generic functions \code{\link{plot.Panmat}} and
#' \code{\link{summary.Panmat}} are defined for a \code{Panmat} object.
#' The input vector \samp{clustering} is attached as the attribute \samp{clustering} to the \code{Panmat}
#' object.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{bClust}}, \code{\link{dClust}}, \code{\link{distManhattan}},
#' \code{\link{distJaccard}}, \code{\link{fluidity}}, \code{\link{chao}},
#' \code{\link{binomixEstimate}}, \code{\link{heaps}}, \code{\link{rarefaction}}.
#' 
#' @examples 
#' # Loading clustering data in the micropan package
#' data(list=c("Mpneumoniae.blast.clustering","Mpneumoniae.domain.clustering"),package="micropan")
#' 
#' # Pan-matrix based on BLAST clustering
#' panmat.blast <- panMatrix(Mpneumoniae.blast.clustering)
#' 
#' # Pan-matrix based on domain sequence clustering
#' panmat.domains <- panMatrix(Mpneumoniae.domain.clustering)
#' 
#' # Plotting the first pan-matrix, and then printing its summary
#' plot(panmat.blast)
#' summary(panmat.blast)
#' 
#' @export panMatrix
#' 
panMatrix <- function( clustering ){
  gids <- sapply( microseq::gregexpr( "GID[0-9]+", names( clustering ), extract=T ), function(x){x[1]} )
  ugids <- sort( unique( gids ) )
  ngids <- length( ugids )
  uclst <- sort( unique( clustering ) )
  nclst <- length( uclst )
  pan.matrix <- matrix( 0, nrow=ngids, ncol=nclst )
  rownames( pan.matrix ) <- ugids
  colnames( pan.matrix ) <- paste( "Cluster", uclst, sep="_" )

  for( i in 1:ngids ){
    idx <- which( gids == ugids[i] )
    clst <- clustering[idx]
    tab <- table( clst )
    idd <- as.numeric( names( tab ) )
    ixx <- which( uclst %in% idd )
    pan.matrix[i,ixx] <- tab
  }
  attr( pan.matrix, "clustering" ) <- clustering
  class( pan.matrix ) <- c( "Panmat", "matrix" )
  return( pan.matrix )
}


#' @rdname generic.Panmat
#' @name plot.Panmat
#' @title Plot and summary of \code{Panmat} objects
#' 
#' @description Generic functions for plotting and printing the content of a \code{Panmat} object.
#' 
#' @param x A \code{Panmat} object, see below.
#' @param object A \code{Panmat} object, see below.
#' @param col The color, default is \samp{"black"}, of interior and borders of the bars in the barplot.
#' @param xlab The label of the X axis.
#' @param ylab The label of the Y axis.
#' @param \dots Optional (graphical) arguments.
#' 
#' @details A \code{Panmat} object contains a pan-matrix, which is the fundamental data structure
#' for pan-genome analyses. It is a small (S3) extension to a \code{matrix}. It has one row for each
#' genome in the study, and one column for each gene cluster. The number in cell \samp{[i,j]} is the
#' number of sequences in genome \samp{i} that belongs to cluster \samp{j}. A \code{Panmat} object is
#' typically created by the function \code{\link{panMatrix}}.
#' 
#' The \code{\link{plot.Panmat}} function will display the content of the \code{Panmat} object as a bar
#' chart showing the number of clusters found in 1,2,...,G genomes, where G is the total number of genomes
#' in the study (rows in \samp{Panmat}).
#' 
#' The \code{\link{summary.Panmat}} function will display a text giving the same information as
#' \code{\link{plot.Panmat}}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}.
#' 
#' @examples # See examples in the Help-file for panMatrix.
#' 
#' @export
plot.Panmat <- function( x, col="black", xlab="Number of genomes", ylab="Number of clusters", ... ){
  # x is a Panmat
  x[which( x > 0, arr.ind=T )] <- 1
  levs <- 1:dim( x )[1]
  y <- table( factor( colSums( x ), levels=levs ) )
  barplot( y, col=col, border=col, names.arg=levs, xlab=xlab, ylab=ylab, ... )
}
#' @rdname generic.Panmat
#' @export
summary.Panmat <- function( object, ... ){
  # object is a Panmat
  object[which( object > 0, arr.ind=T )] <- 1
  levs <- 1:nrow( object )
  y <- table( factor( colSums( object ), levels=levs ) )
  for( i in 1:length( levs ) ){
    cat( y[i], "clusters found in", levs[i], "genomes\n" )
  }
}


#' @name panpca
#' @title Principal component analysis of a pan-matrix
#' 
#' @description Computes a principal component decomposition of a pan-matrix, with possible
#' scaling and weightings.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' @param scale An optional scale to control how copy numbers should affect the distances.
#' @param weights Vector of optional weights of gene clusters.
#' 
#' @details A principal component analysis (PCA) can be computed for any matrix, also a pan-matrix.
#' The principal components will in this case be linear combinations of the gene clusters. One major
#' idea behind PCA is to truncate the space, e.g. instead of considering the genomes as points in a
#' high-dimensional space spanned by all gene clusters, we look for a few \sQuote{smart} combinations
#' of the gene clusters, and visualize the genomes in a low-dimensional space spanned by these directions.
#' 
#' The \samp{scale} can be used to control how copy number differences play a role in the PCA. Usually
#' we assume that going from 0 to 1 copy of a gene is the big change of the genome, and going from 1 to
#' 2 (or more) copies is less. Prior to computing the PCA, the \samp{pan.matrix} is transformed according
#' to the following affine mapping: If the original value in \samp{pan.matrix} is \samp{x}, and \samp{x}
#' is not 0, then the transformed value is \samp{1 + (x-1)*scale}. Note that with \samp{scale=0.0}
#' (default) this will result in 1 regardless of how large \samp{x} was. In this case the PCA only
#' distinguish between presence and absence of gene clusters. If \samp{scale=1.0} the value \samp{x} is
#' left untransformed. In this case the difference between 1 copy and 2 copies is just as big as between
#' 1 copy and 0 copies. For any \samp{scale} between 0.0 and 1.0 the transformed value is shrunk towards
#' 1, but a certain effect of larger copy numbers is still present. In this way you can decide if the PCA
#' should be affected, and to what degree, by differences in copy numbers beyond 1.
#' 
#' The PCA can also up- or downweight some clusters compared to others. The vector \samp{weights} must
#' contain one value for each column in \samp{pan.matrix}. The default is to use flat weights, i.e. all
#' clusters count equal. See \code{\link{geneWeights}} for alternative weighting strategies.
#' 
#' The functions \code{\link{plotScores}} and \code{\link{plotLoadings}} can be used to visualize the
#' results of \code{\link{panpca}}.
#' 
#' @return A \code{Panpca} object is returned from this function. This is a small (S3) extension of
#' a \code{list} with elements \samp{Evar}, \samp{Scores}, \samp{Loadings}, \samp{Scale} and \samp{Weights}. 
#' 
#' \samp{Evar} is a vector with one number for each principal component. It contains the relative
#' explained variance for each component, and it always sums to 1.0. This value indicates the importance of
#' each component, and it is always in descending order, the first component being the most important.
#' The \samp{Evar} is typically the first result you look at after a PCA has been computed, as it indicates
#' how many components (directions) you need to capture the bulk of the total variation in the data.
#' 
#' \samp{Scores} is a matrix with one column for each principal component and one row for each genome. The
#' columns are ordered corresponding to the elements in \samp{Evar}. The scores are the coordinates of
#' each genome in the principal component space. See \code{\link{plotScores}} for how to visualize genomes
#' in the score-space.
#' 
#' \samp{Loadings} is a matrix with one column for each principal component and one row for each gene
#' cluster. The columns are ordered corresponding to the elements in \samp{Evar}. The loadings are the
#' contribution from each original gene cluster to the principal component directions. NOTE: Only gene
#' clusters having a non-zero variance is used in a PCA. Gene clusters with the same value for every
#' genome have no impact and are discarded from the \samp{Loadings}. See \code{\link{plotLoadings}} for
#' how to visualize gene clusters in the loading space.
#' 
#' \samp{Scale} and \samp{Weights} are copies of the corresponding input arguments.
#' 
#' The generic functions \code{\link{plot.Panpca}} and \code{\link{summary.Panpca}}
#' are available for \code{Panpca} objects.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{plotScores}}, \code{\link{plotLoadings}}, \code{\link{panTree}},
#' \code{\link{distManhattan}}, \code{\link{geneWeights}}.
#' 
#' @examples 
#' # Loading two Panmat objects in the micropan package
#' data(list=c("Mpneumoniae.blast.panmat","Mpneumoniae.domain.panmat"),package="micropan")
#' 
#' # Panpca based on a BLAST clustering Panmat object
#' ppca.blast <- panpca(Mpneumoniae.blast.panmat)
#' plot(ppca.blast) # The generic plot function
#' plotScores(ppca.blast) # A score-plot
#' 
#' # Panpca based on domain sequence clustering Panmat object
#' w <- geneWeights(Mpneumoniae.domain.panmat,type="shell")
#' ppca.domains <- panpca(Mpneumoniae.domain.panmat,scale=0.5,weights=w)
#' summary(ppca.domains)
#' plotLoadings(ppca.domains)
#' 
#' @export panpca
#' 
panpca <- function( pan.matrix, scale=0.0, weights=rep( 1, dim( pan.matrix )[2] ) ){
  if( (scale>1)|(scale<0) ){
    warning( "scale should be between 0.0 and 1.0, using scale=0.0" )
    scale <- 0.0
  }
  idx <- which( pan.matrix > 0, arr.ind=T )
  pan.matrix[idx] <- 1 + (pan.matrix[idx]-1)*scale
  pan.matrix <- pan.matrix * matrix( weights, nrow=dim( pan.matrix )[1], ncol=dim( pan.matrix )[2], byrow=T )
  idx <- which( apply( pan.matrix, 2, sd ) > 0 )
  X <- pan.matrix[,idx]

  pca <- prcomp( X )
  evar <- pca$sdev^2/sum( pca$sdev^2 )
  scores <- pca$x
  loadings <- pca$rotation

  pan.pca <- list( Evar=evar, Scores=scores, Loadings=loadings, Scale=scale, Weights=weights )
  class( pan.pca ) <- c( "Panpca", "list" )
  return( pan.pca )
}


#' @rdname generic.Panpca
#' @name plot.Panpca
#' @title Plot and summary of \code{Panpca} objects
#' 
#' @description Generic functions for \code{Panpca} objects.
#' 
#' @param x A \code{Panpca} object, see below.
#' @param object A \code{Panpca} object, see below.
#' @param cum Logical, default is \samp{FALSE}, indicating if explained variance should be plotted
#' per component or cumulative.
#' @param col Color, default is \samp{"black"}, of interior and border of bars in the barplot.
#' @param \dots Optional graphical arguments.
#' 
#' @details A \code{Panpca} object contains the results from a principal component analysis (PCA) on
#' a pan-matrix, and is the output from the function \code{\link{panpca}}. It is a small (S3) extension
#' of a \code{list}, and contains the elements \samp{Evar}, \samp{Scores}, \samp{Loadings}, \samp{Scale}
#' and \samp{Weights}.
#' 
#' The basic idea of a PCA is to find alternative directions in the space spanned by the pan-matrix
#' columns, in order to be able to visualize or in other ways extract the most relevant information in
#' a small number of dimensions. The variable \samp{Evar} contains the explained variance for each
#' principal component, scaled such that summed over all components it is 1.0. This quantity indicates
#' the importance of each component, larger values of \samp{Evar} indicates directions (components) with
#' more information.
#' 
#' The \code{\link{plot.Panpca}} function shows the \samp{Evar} values in a barplot. You can either plot
#' the \samp{Evar} value of each component separately (\samp{cum=FALSE}) or the cumulative value
#' (\samp{cum=TRUE}). This is the basic plot to follow any principal component decomposition, since it
#' tells you how many components you need to capture the bulk of the information in the data. If e.g.
#' component 1, 2 and 3 have \samp{Evar} values of 0.4, 0.3 and 0.2, respectively, it means these three
#' direction capture 90\% (0.4+0.3+0.2=0.9) of all the variation in the data. For some pan-matrices almost
#' all variation can be found in the very few first directions, but more often it is scattered between many.
#' See \code{\link{plotScores}} and \code{\link{plotLoadings}} for other informative graphical displays of
#' a \code{Panpca} object.
#' 
#' The \code{\link{summary.Panpca}} function will print the same information as plotted by
#' \code{\link{plot.Panpca}}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panpca}}, \code{\link{plotScores}}, \code{\link{plotLoadings}}.
#' 
#' @examples # See examples in the Help-file for panpca.
#' 
#' @export
plot.Panpca <- function( x, cum=FALSE, col="black", ... ){
  Panpca <- x
  if( cum ){
    y <- cumsum( Panpca$Evar )
  } else {
    y <- Panpca$Evar
  }
  barplot( y, col=col, names.arg=1:length( y ), xlab="Principal components", ylab="Relative explained variance", ... )
}
#' @rdname generic.Panpca
#' @export
summary.Panpca <- function( object, ... ){
  for( i in 1:length( object$Evar ) ){
    cat( "Principal component ", i, " explains ", round( 1000*object$Evar[i] )/10, "% of the variation\n", sep="" )
  }
}


#' @rdname scores.loadings
#' @name plotScores
#' @title Plotting scores and loadings in a \code{Panpca} object
#' 
#' @description Creates informative plots for a principal component analysis of a pan-matrix.
#' 
#' @param pan.pca A \code{Panpca} object, see \code{\link{panpca}} for details.
#' @param x The component to display along the horizontal axis.
#' @param y The component to display along the vertical axis.
#' @param show.labels Logical indicating if labels should be displayed.
#' @param labels Alternative labels to use in the score-plot, see below.
#' @param col Colors for the points/labels, see below.
#' @param pch Marker type, see \code{\link{points}}.
#' @param \dots Additional arguments passed on to \code{points} or \code{text} (if labels are specified).
#' 
#' @details A \code{Panpca} object contains the results of a principal component analysis on a pan-matrix,
#' see \code{\link{panpca}} for details.
#' 
#' The \code{\link{plotScores}} gives a visual overview of how the genomes are positioned relative to
#' each other in the pan-genome space. The score-matrix of a \code{Panpca} has one row for each genome.
#' The original pan-matrix also has one row for each genome. Two genomes can be compared by their
#' corresponding rows in the pan-matrix, but can also be compared by their rows in the score-matrix,
#' and the latter matrix has (much) fewer columns designed to contain maximum of the original data
#' variation. A plot of the scores will give an approximate overview of how the genomes are located
#' relative to each other.
#' 
#' The \code{\link{plotLoadings}} gives a visual overview of how the gene clusters affect the principal
#' components. The loadings is a matrix with one row for each of the original non-core gene clusters
#' (core gene clusters have no variation across genomes). Clusters located close to the origin have
#' little impact. Clusters far from the origin has high impact, indicating they separate groups of genomes.
#' 
#' These two plots together can reveal information about the pan-genome: The score-plot shows if genomes
#' are grouped/separated, and the loading-plot can then tell you which gene clusters have high impact on
#' this grouping/separation.
#' 
#' The arguments \samp{x} and \samp{y} can be used to plot other components than component 1 and 2
#' (which is always the most informative). In some cases more components are needed to establish a
#' good picture, i.e. the explained variance is low for component 1 and 2 (see \code{\link{plot.Panpca}}
#' for more on explained variance). It is quite common to plot component 1 versus 2, then 1 versus 3
#' and finally 2 versus 3.
#' 
#' The argument \samp{show.labels} can be used to turn off the display of labels, only markers (dots)
#' will appear.
#' 
#' In \code{\link{plotScores}} you can specify alternative labels in \samp{labels}. By default, the
#' GID-tag is used for each genome. You can supply a vector of alternative labels. The labels may be
#' in any order, but the vector must be named by the GID-tags, i.e. each element in \samp{labels} must
#' have a name which is a valid GID-tag for some genome. This is necessary to ensure the alternative
#' labels are placed correctly in the score-space.
#' 
#' There is no alternative labelling of loading-plots, since the gene clusters lack a GID-tag-like system.
#' You can, however, change the gene cluster names by editing the column names of the pan-matrix directly
#' before you do the \code{\link{panpca}}.
#' 
#' You may color each label/marker individually. In \code{\link{plotScores}} you can again supply a vector
#' of colors, and name every element with a GID-tag to make certain they are used correctly. In
#' \code{\link{plotLoadings}} you can supply a vector of colors, but you must arrange them in proper
#' order yourself.
#' 
#' Additional arguments are passed on to \code{\link{text}} if \samp{show.labels=TRUE} and to
#' \code{\link{points}} if \samp{show.labels=FALSE}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panpca}}, \code{\link{plot.Panpca}}.
#' 
#' @examples 
#' # Loading a Panmat object in the micropan package
#' data(list=c("Mpneumoniae.blast.panmat","Mpneumoniae.domain.panmat"),package="micropan")
#' ppca.blast <- panpca(Mpneumoniae.blast.panmat)
#' 
#' # Plotting scores and loadings
#' plotScores(ppca.blast) # A score-plot
#' plotLoadings(ppca.blast) # A loading plot
#' 
#' # Plotting score with alternative labels and colors
#' data(list="Mpneumoniae.table",package="micropan")
#' labels <- Mpneumoniae.table$Strain
#' names(labels) <- Mpneumoniae.table$GID.tag
#' cols <- Mpneumoniae.table$Color
#' names(cols) <- Mpneumoniae.table$GID.tag
#' plotScores(ppca.blast,labels=labels,col=cols)
#' 
#' @export
plotScores <- function( pan.pca, x=1, y=2, show.labels=TRUE, labels=NULL, col="black", pch=16, ... ){
  Z <- pan.pca$Scores
  xr <- range( Z[,x] )
  yr <- range( Z[,y] )
  args <- list(...)
  ii <- match( c("xlab","ylab"), names( args ) )
  if( is.na(ii[1]) ){
    xlab=paste( "PC", x, " (", round( 100*pan.pca$Evar[x] ), "%)", sep="" )
  } else {
    xlab <- args[[ii[1]]]
  }
  if( is.na(ii[2]) ){
    ylab=paste( "PC", y, "(", round( 100*pan.pca$Evar[y] ), "%)", sep="" )
  } else {
    ylab <- args[[ii[2]]]
  }
  
  plot( xr, c(0,0), type="l", col="gray", xlim=xr, ylim=yr, xlab=xlab, ylab=ylab )
  points( c(0,0), yr, type="l", col="gray" )
  
  pm.gid <- rownames( Z )
  if( length( col )>1 ){
    if( is.null( names( col ) ) ) stop( "Each element in col must be named by its GID.tag" )
    gid <- names( col )
    idx <- match( pm.gid, gid )
    if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the pan-matrix" )
    cols <- col[idx]
    if( length( cols ) != length( pm.gid ) ) stop( "The number of elements in col does not correspond to the number of genomes in the pan-matrix" )
  } else {
    cols <- rep( col, length.out=length( pm.gid ) )
  }
  if( is.null( labels ) ){
    labs <- rownames( Z )
  } else {
    if( is.null( names( labels ) ) ) stop( "Each element in labels must be named by its GID.tag" )
    gid <- names( labels )
    idx <- match( pm.gid, gid )
    if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the pan-matrix" )
    labs <- as.character( labels[idx] )
    if( length( labs ) != length( pm.gid ) ) stop( "The number of elements in labels does not correspond to the number of genomes in the pan-matrix" )
  }
  if( show.labels ){
    text( Z[,x], Z[,y], labs, col=cols, ... )
  } else {
    points( Z[,x], Z[,y], pch=pch, col=cols, ... )
  }
}
#' @rdname scores.loadings
#' @export
plotLoadings <- function( pan.pca, x=1, y=2, show.labels=TRUE, col="black", pch=16, ... ){
  L <- pan.pca$Loadings
  xr <- range( L[,x] )
  yr <- range( L[,y] )
  args <- list(...)
  ii <- match( c("xlab","ylab"), names( args ) )
  if( is.na(ii[1]) ){
    xlab=paste( "PC", x, " (", round( 100*pan.pca$Evar[x] ), "%)", sep="" )
  } else {
    xlab <- args[[ii[1]]]
  }
  if( is.na(ii[2]) ){
    ylab=paste( "PC", y, "(", round( 100*pan.pca$Evar[y] ), "%)", sep="" )
  } else {
    ylab <- args[[ii[2]]]
  }
  
  plot( xr, c(0,0), type="l", col="gray", xlim=xr, ylim=yr, xlab=xlab, ylab=ylab )
  points( c(0,0), yr, type="l", col="gray" )
  
  if( show.labels ){
    text( L[,x], L[,y], rownames( L ), col=col, ... )
  } else {
    points( L[,x], L[,y], pch=pch, col=col, ... )
  }
}
#' @name chao
#' @title The Chao lower bound estimate of pan-genome size
#' 
#' @description Computes the Chao lower bound estimated number of gene clusters in a pan-genome.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' 
#' @details The size of a pan-genome is the number of gene clusters in it, both those observed and those
#' not yet observed.
#' 
#' The input \samp{pan.matrix} is a \code{Panmat} object, i.e. it is a matrix with one row for each
#' genome and one column for each observed gene cluster in the pan-genome. See \code{\link{panMatrix}}
#' for how to construct such objects.
#' 
#' The number of observed gene clusters is simply the number of columns in \samp{pan.matrix}. The
#' number of gene clusters not yet observed is estimated by the Chao lower bound estimator (Chao, 1987).
#' This is based solely on the number of clusters observed in 1 and 2 genomes. It is a very simple and
#' conservative estimator, i.e. it is more likely to be too small than too large. 
#' 
#' @return The function returns an integer, the estimated pan-genome size. This includes both the number
#' of gene clusters observed so far, as well as the estimated number not yet seen.
#' 
#' @references Chao, A. (1987). Estimating the population size for capture-recapture data with unequal
#' catchability. Biometrics, 43:783-791.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{panMatrix}}, \code{\link{binomixEstimate}}.
#' 
#' @examples 
#' # Loading a Panmat object in the micropan package
#' data(list="Mpneumoniae.blast.panmat",package="micropan")
#' 
#' # Estimating the pan-genome size using the Chao estimator
#' chao.pansize <- chao(Mpneumoniae.blast.panmat)
#' 
#' @export
chao <- function( pan.matrix ){  
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  y <- table( factor( colSums( pan.matrix ), levels=1:dim( pan.matrix )[1] ) )
  if( y[2]==0 ){
    stop( "Cannot compute Chao estimate since there are 0 gene clusters observed in 2 genomes!\n" )
  } else {
    pan.size <- round( sum( y ) + y[1]^2/(2*y[2]) )
    names( pan.size ) <- NULL
    return( pan.size )
  }
}


#' @name heaps
#' @title Heaps law estimate
#' 
#' @description Estimating if a pan-genome is open or closed based on a Heaps law model.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' @param n.perm The number of random permutations of genome ordering.
#' 
#' @details An open pan-genome means there will always be new gene clusters observed as long as new genomes
#' are being sequenced. This may sound controversial, but in a pragmatic view, an open pan-genome indicates
#' that the number of new gene clusters to be observed in future genomes is \sQuote{large} (but not literally
#' infinite). Opposite, a closed pan-genome indicates we are approaching the end of new gene clusters. 
#' 
#' This function is based on a Heaps law approach suggested by Tettelin et al (2008). The Heaps law model
#' is fitted to the number of new gene clusters observed when genomes are ordered in a random way. The model
#' has two parameters, an intercept and a decay parameter called \samp{alpha}. If \samp{alpha>1.0} the
#' pan-genome is closed, if \samp{alpha<1.0} it is open.
#' 
#' The number of permutations, \samp{n.perm}, should be as large as possible, limited by computation time.
#' The default value of 100 is certainly a minimum.
#' 
#' Word of caution: The Heaps law assumes independent sampling. If some of the genomes in the data set
#' form distinct sub-groups in the population, this may affect the results of this analysis severely.
#' 
#' @return A vector of two estimated parameters: The \samp{Intercept} and the decay parameter \samp{alpha}.
#' If \samp{alpha<1.0} the pan-genome is open, if \samp{alpha>1.0} it is closed.
#' 
#' @references Tettelin, H., Riley, D., Cattuto, C., Medini, D. (2008). Comparative genomics: the
#' bacterial pan-genome. Current Opinions in Microbiology, 12:472-477.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{binomixEstimate}}, \code{\link{chao}}, \code{\link{rarefaction}}.
#' 
#' @examples 
#' # Loading a Panmat object in the micropan package 
#' data(list="Mpneumoniae.blast.panmat",package="micropan")
#' 
#' # Estimating population openness
#' h.est <- heaps(Mpneumoniae.blast.panmat,n.perm=500)
#' if(h.est[2]>1){
#'   cat("Population is closed with alpha =",h.est[2], "\n")
#' } else {
#'   cat("Population is open with alpha =",h.est[2], "\n")
#' }
#' 
#' @importFrom stats optim
#' 
#' @export
heaps <- function( pan.matrix, n.perm=100 ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  nmat <- matrix( 0, nrow=(ng-1), ncol=n.perm )
  cat( "permuting:\n" )
  for( i in 1:n.perm ){
    cm <- apply( pan.matrix[sample( ng ),], 2, cumsum )
    nmat[,i] <- rowSums( (cm==1)[2:ng,] & (cm==0)[1:(ng-1),] )
    cat( "." )
  }
  cat( "\n" )
  x <- rep( (2:dim( pan.matrix )[1]), times=n.perm )
  y <- as.numeric( nmat )
  p0 <- c( mean( y[which( x == 2 )] ), 1 )
  fit <- optim( p0, objectFun, gr=NULL, x, y, method="L-BFGS-B", lower=c(0,0), upper=c(10000,2) )
  p.hat <- fit$par
  names( p.hat ) <- c( "Intercept", "alpha" )
  return( p.hat )
}

objectFun <- function( p, x, y ){
  y.hat <- p[1]*x^(-p[2])
  J <- sqrt( sum( (y - y.hat)^2 ) )/length( x )
  return( J )
}
#' @name rarefaction
#' @title Rarefaction curves for a pan-genome
#' 
#' @description Computes rarefaction curves for a number of random permutations of genomes.
#' 
#' @param pan.matrix A \code{Panmat} object, see \code{\link{panMatrix}} for details.
#' @param n.perm The number of random genome orderings to use. If \samp{n.perm=1} the fixed order of
#' the genomes in \samp{pan.matrix} is used.
#' 
#' @details A rarefaction curve is simply the cumulative number of unique gene clusters we observe as
#' more and more genomes are being considered. The shape of this curve will depend on the order of the
#' genomes. This function will typically compute rarefaction curves for a number of (\samp{n.perm})
#' orderings. By using a large number of permutations, and then averaging over the results, the effect
#' of any particular ordering is smoothed away.
#' 
#' The averaged curve illustrates how many new gene clusters we observe for each new genome. If this
#' levels out and becomes flat, it means we expect few, if any, new gene clusters by sequencing more
#' genomes. The function \code{\link{heaps}} can be used to estimate population openness based on this
#' principle.
#' 
#' @return This function returns a \code{Rarefac} object, which is a small extension to a matrix. The
#' generic functions \code{\link{plot.Rarefac}} and \code{\link{summary.Rarefac}}
#' are available for such objects.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{heaps}}, \code{\link{panMatrix}}, \code{\link{plot.Rarefac}},
#' \code{\link{summary.Rarefac}}.
#' 
#' @examples 
#' # Loading two Panmat objects in the micropan package 
#' data(list=c("Mpneumoniae.blast.panmat","Mpneumoniae.domain.panmat"),package="micropan")
#' 
#' # Rarefaction based on a BLAST clustering Panmat object
#' rarefac.blast <- rarefaction(Mpneumoniae.blast.panmat,n.perm=100)
#' plot(rarefac.blast)
#' 
#' # Rarefaction based on domain sequence clustering Panmat object
#' rarefac.domains <- rarefaction(Mpneumoniae.domain.panmat,n.perm=1000)
#' summary(rarefac.domains)
#' 
#' @export
rarefaction <- function( pan.matrix, n.perm=1 ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  nmat <- matrix( 0, nrow=ng, ncol=n.perm )
  cm <- apply( pan.matrix, 2, cumsum )
  nmat[,1] <- rowSums( cm>0 )
  if( n.perm > 1 ){
    cat( "permuting:\n" )
    for( i in 2:n.perm ){
      cm <- apply( pan.matrix[sample( ng ),], 2, cumsum )
      nmat[,i] <- rowSums( cm>0 )
      cat( "." )
      if( (i/100)==round(i/100) ) cat( "\n" )
    }
    cat( "\n" )
  }
  rownames( nmat ) <- paste( 1:ng, "genomes" )
  colnames( nmat ) <- paste( "Permutation", 1:n.perm )
  class( nmat ) <- c( "Rarefac", "matrix" )
  return( nmat )
}



#' @rdname generic.Rarefac
#' @name plot.Rarefac
#' @title Plot and summary of \code{Rarefac} objects
#' 
#' @description Generic functions for \code{Rarefac} object.
#' 
#' @param x A \code{Rarefac} object, see below.
#' @param object A \code{Rarefac} object, see below.
#' @param type Type of plot, default is \samp{"b"}, giving markers with lines between.
#' @param pch Marker type, default is \samp{16}, a filled circle.
#' @param xlab Text for horizontal axis.
#' @param ylab Text for vertical axis.
#' @param \dots Optional graphical arguments.
#' 
#' @details A \code{Rarefac} object is a small (S3) extension to a matrix. The first column contains
#' the cumulative number of unique gene clusters found when considering 1,2,...,G genomes in a pan-matrix.
#' Thus, the \code{Rarefac} object is a matrix with G rows. Any additional columns will hold similar
#' numbers, but for random shufflings of the genome's ordering. A \code{Rarefac} object is typically
#' created by the function \code{\link{rarefaction}}.
#' 
#' The \code{\link{plot.Rarefac}} function will display the content of the \code{Rarefac} object as a plot
#' of the mean value in rows 1,2,...,G, where G is the total number of genomes in the study.
#' 
#' The \code{\link{summary.Rarefac}} function will display a text giving the same information as
#' \code{\link{plot.Rarefac}}.
#' 
#' @author Lars Snipen and Kristian Hovde Liland.
#' 
#' @seealso \code{\link{rarefaction}}, \code{\link{heaps}}.
#' 
#' @examples # See examples in the Help-file for rarefaction.
#' 
#' @export
plot.Rarefac <- function( x, type="b", pch=16, xlab="Genomes", ylab="Number of unique gene clusters", ... ){
  Rarefac <- x
  plot( 1:dim( Rarefac )[1], rowMeans( Rarefac ), type=type, pch=pch, xlab=xlab, ylab=ylab, ... )
}
#' @rdname generic.Rarefac
#' @export
summary.Rarefac <- function( object, ... ){
  cat( "For", 1, "genome we observe on average", round( mean( object[1,] ) ), "unique gene clusters\n" )
  for( i in 2:nrow( object ) ){
    cat( "For", i, "genomes we observe on average", round( mean( object[i,] ) ), "unique gene clusters\n" )
  }
}

# define function as.panmat
as.panmat <- function(x) { # y is the proteinortho matrix, newmat is the new matrix
	y <- x[,4:dim(x)[2]]
	newmat <- matrix(ncol=dim(y)[2], nrow=dim(y)[1])
	for(i in 1:dim(y)[1]) {
		for(j in 1:dim(y)[2]) {
			entry <- as.character(y[i,j])
			if(entry == "*") {newmat[i,j] <- as.numeric(0)}
			else{
			spl <- strsplit(entry, ",")
			l <- length(spl[[1]])
			newmat[i,j] <- as.numeric(l)
			}
		}
	}
	names <- paste("cluster", c(1:dim(y)[1]), sep="_")
	colnames(newmat) <- colnames(y)
	row.names(newmat) <- names
	class(newmat) <- c("Panmat", "matrix")
	return (t(newmat))
}
# end

# define pairwise fluidity function
fluidity.pw <- function(panmat) {
	newmat <- matrix(ncol=dim(panmat)[1], nrow=dim(panmat)[1])
	for(i in 1:dim(panmat)[1]) {
		for(j in 1:dim(panmat)[1]) {
			fl <- fluidity(panmat[c(i,j),])[[1]]
			newmat[i,j] <- fl
		}	
	}
	return(newmat)
}
# end


