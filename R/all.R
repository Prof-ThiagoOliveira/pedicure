#' Check a pedigree file
#'
#' Check a pedigree file for missing or out of order founders.
#'
#' @param pedigree
#'
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively. The row
#' giving the pedigree of an individual must appear before any row
#' where that individual appears as a parent. Founders use 0 (zero) or
#' NA in the parental columns.
#'
#' @param fgen
#'
#' An optional list of length 2 where \code{fgen[[1]]} is a character
#' string naming the column in \code{pedigree} that contains the level
#' of selfing or the level of inbreeding of an individual. In
#' \code{pedigree[,fgen[[1]]]}, 0 indicates a simple cross, 1
#' indicates selfed once, 2 indicates selfed twice, etc. A value
#' between 0 and 1 for a base individual is taken as its inbreeding
#' value. If the pedigree has implicit individuals (they appear as
#' parents but not as individuals), they will be assumed base
#' non-inbred individuals unless their inbreeding level is set with
#' \code{fgen[[2]]} where \code{0<fgen[[2]]<1} is the inbreeding level
#' of such individuals.
#'
#' @param gender
#'
#' An optional character string naming the column of \code{pedigree}
#' that codes for the gender of an individual. \code{pedigree[,
#' gender]} is coerced to a factor and must only have two (arbitrary)
#' levels, the first of which is taken to mean "male". An inverse
#' relationship matrix is formed for the X chromosome as described by
#' \cite{Fernando and Grossman (1990)} for species where the male is
#' XY and the female is XX.
#'
#' @param mv
#'
#' Missing value indicator; elements of \code{pedigree} that exactly
#' match any element of \code{mv} are treated as missing.
#'
#' @param verbose
#' 
#' If \code{TRUE} details of pedigree insertions and relocations are printed;
#' the default is \code{FALSE} to supress printing.
#' 
#' @return
#'
#' A data frame containing the expanded or re-ordered pedigree.
#'
#' @references
#'
#' % bibentry: Fernando:1990
#'
checkPedigree <- function(pedigree,fgen=list(character(0),0.01),
                          gender=character(0), mv=c("0","*"," "), verbose=FALSE)
{
  is.missing <- function(pcol,mv) {
    isna <- is.na(pcol) | is.element(as.character(pcol),mv) | trimws(as.character(pcol)) == ""

    return(isna)
  }
  if(any(dup <- duplicated(as.character(pedigree[,1])))) {
    stop(cat("Duplicated individuals","\n",as.character(pedigree[,1])[dup],"\n"))
  }
  which <- seq(1,nrow(pedigree))[as.character(pedigree[,1])==as.character(pedigree[,2]) |
                                 as.character(pedigree[,1])==as.character(pedigree[,3])]
  if(length(which <- which[!is.na(which)])) {
    cat(paste(which,pedigree[,1][which]),sep="\n")
    stop("Individuals appear as their own parent")
  }

  fmode <- FALSE
  col4 <- -1 # used as a flag in sortPed
  col4Name <- NULL
  col4Type <- character(0)
  fg <- as.double(fgen[[2]])

  if(length(fgen[[1]]) > 0) {
    if(length(col4 <- pedigree[,fgen[[1]]]) == 0)
      stop("Argument 'fgen' specifies a non-existent dataframe column.")
    else {
      fmode <- TRUE
      col4Name <- fgen[[1]]
      col4Type <- "fgen"
      w <- is.na(pedigree[,col4Name])
      if(sum(as.numeric(w)) > 0) {
        if(sum(as.numeric(!(is.missing(pedigree[w,2],mv) |
                            is.missing(pedigree[w,3],mv)))) > 0)
          stop("fgen is 'NA' for non-base individuals")
      }
    }
  }
  if(length(gender) > 0) {
    if(fmode) stop("Cannot specify both 'fgen' and 'gender'")
    xlink <- TRUE
    if(length(pedigree[,gender]) == 0)
      stop("Argument 'gender' specifies non-existent dataframe column.")
    else if(!is.factor(pedigree[,gender]))
      pedigree[,gender] <- factor(pedigree[,gender])
    sx <- levels(pedigree[,gender])
    if(length(sx) != 2)
      stop("'gender' must only have two levels.")
    col4 <- as.numeric(pedigree[,gender])
    col4Name <- gender
    col4Type <- "gender"
  }

  lped <- nrow(pedigree)

  ## Get levels ignoring missing value character(s)
  charPed <- c(as.character(pedigree[,1]),
               as.character(pedigree[,2]),
               as.character(pedigree[,3]))
  lvls <- unique(charPed)
  lvls <- lvls[!is.missing(lvls,mv)]
  charPed <- NULL

  nan <- length(lvls)
  ## convert to integer
  memumdad <- names(pedigree)[1:3]
  pedigree <- matrix(c(match(pedigree[,1],lvls,nomatch=0),
                       match(pedigree[,2],lvls,nomatch=0),
                       match(pedigree[,3],lvls,nomatch=0)),nrow=lped,byrow=FALSE)
  xtras <- numeric(0)
  if(lped < nan) {
    xmum <- match(pedigree[,2],pedigree[,1])
    xmum <- unique(pedigree[is.na(xmum),2])
    xmum <- xmum[!is.na(xmum) & xmum > 0]
    xdad <- match(pedigree[,3],pedigree[,1])
    xdad <- unique(pedigree[is.na(xdad),3])
    xdad <- xdad[!is.na(xdad) & xdad > 0]
    xtras <- unique(c(xmum,xdad))
    nxtra <- length(xtras)
    ## check
    if(nxtra != nan-lped)
      stop("programming error in checkPedigree")
    xtra <- matrix(c(xtras,rep(0,2*nxtra)),nrow=nxtra,ncol=3,byrow=FALSE)
    pedigree <- rbind(xtra,pedigree)
    if(fmode) col4 <- c(rep(fg,(nan-lped)),col4)
  }

  ## No need for insertions (now, lped==nan) but "sortPed" puts founders before offspring
  ## stop on error
  params <- list(levels=lvls, fg=fg, col4=col4, verbose=as.numeric(verbose))
  out <- .Call("sortped", params, pedigree, package="pedicure")
  pedigree <- matrix(out[[1]], ncol=3, byrow=FALSE)
  col4 <- out[[2]]
  which <- pedigree > 0
  pedigree[which] <- lvls[pedigree[which]]
  pedigree <- data.frame(pedigree,stringsAsFactors=FALSE)

  if(length(col4Name)) {
    pedigree[[col4Name]] <- col4
    if(fmode) {
      pedigree[is.na(pedigree[,4]),4] <- fg
      x <- pedigree[,4]>1 & (as.character(pedigree[,2]) == "0" | as.character(pedigree[,3]) == "0")
      if(length(x)) pedigree[x,4] <- 1.0 - 0.5^pedigree[x,4]
    }
  }
  names(pedigree) <- c(memumdad,col4Name)
  attr(pedigree,"rowNames") <- pedigree[,1]
  if(length(col4Type))
    attr(pedigree,"col4Type") <- col4Type
  if(length(xtras)) {
    attr(pedigree,"Insertions") <- lvls[xtras]
    if(verbose) {
      print("Insertions")
      print(lvls[xtras])
    } else {
      message(length(xtras), " insertions; see 'Insertions' attribute.\n")
    }
  }
  return(pedigree)
}
#' Calculate an inverse relationship matrix.
#'
#' Generates an inverse relationship matrix in sparse triplet form
#' from a pedigree data frame.
#'
#' Uses the method of \cite{Meuwissen and Luo, 1992} to compute the inverse
#' relationship matrix directly from the pedigree.
#'
#' @param pedigree
#'
#' A data frame where the first three columns correspond to the
#' identifiers for the individual, male parent and female parent,
#' respectively. The row giving the pedigree of an individual must
#' appear before any row where that individual appears as a
#' parent. Founders use 0 (zero) or NA in the parental columns.
#'
#' @param fgen
#'
#' An optional list of length 2 where \code{fgen[[1]]} is a character
#' string naming the column in \code{pedigree} that contains the level
#' of selfing or the level of inbreeding of an individual. In
#' \code{pedigree[,fgen[[1]]]}, 0 indicates a simple cross, 1
#' indicates selfed once, 2 indicates selfed twice, etc. A value
#' between 0 and 1 for a base individual is taken as its inbreeding
#' value. If the pedigree has implicit individuals (they appear as
#' parents but not as individuals), they will be assumed base
#' non-inbred individuals unless their inbreeding level is set with
#' \code{fgen[[2]]}, where \code{0 < fgen[[2]] < 1} is the inbreeding
#' level of such individuals.
#'
#' @param gender
#'
#' An optional character string naming the column of \code{pedigree}
#' that codes for the gender of an individual. \code{pedigree[,
#' gender]} is coerced to a factor and must only have two (arbitrary)
#' levels, the first of which is taken to mean "male". An inverse
#' relationship matrix is formed for the X chromosome as described by
#' \cite{Fernando and Grossman, 1990} for species where the male is XY
#' and the female is XX.
#'
#' @param groups
#'
#' An integer scalar (\emph{g}) indicating genetic groups in the
#' pedigree. The first \emph{g} lines of the pedigree identify the
#' genetic groups (with zero in both the male and female parent
#' columns). All other rows must specify one of the genetic groups as
#' the male or female parent if the actual parent is unknown. The
#' default is \eqn{g=0}.
#'
#' @param groupOffset
#'
#' A numeric scalar \emph{e > 0} added to the diagonal elements of
#' \eqn{A^{-1}} pertaining to \code{groups}, shrinking the group
#' effects by \eqn{e}. When a constant is added, no adjustment of the
#' degrees of freedom is made for genetic groups.  Set to -1 to add no
#' offset but to suppress insertion of constraints where empty groups
#' appear; the empty groups are then not counted in the degress of
#' freedom adjustment. The default is \eqn{e=0}.
#'
#' @param selfing
#'
#' A numeric scalar (\emph{s}) allowing for partial selfing when the
#' third field of \code{pedigree} is unknown. It indicates that
#' progeny from a cross where the male parent is unknown is assumed to
#' be from selfing with probability \emph{s} and from outcrossing with
#' probability \emph{(1-s)}. This is appropriate in some forestry tree
#' breeding studies where seed collected from a tree may have been
#' pollinated by the mother tree or pollinated by some other tree
#' (\cite{Dutkowski and Gilmour, 2001}). Do not use the \code{selfing}
#' argument in conjunction with \code{inBreed} or \code{mgs}.
#'
#' @param inBreed
#'
#' A numeric scalar (default \code{NA}) giving the inbreeding
#' coefficient for \emph{base} individuals. This argument generates
#' the numerator relationship matrix for inbred lines. Each cross is
#' assumed to be selfed several times to stabilize as an inbred line
#' as is usual for cereal crops, for example, before being evaluated
#' or crossed with another line.  Since inbreeding is usually
#' associated with strong selection, it is not obvious that a pedigree
#' assumption of covariance of 0.5 between parent and offspring
#' actually holds. The \code{inBreed} argument cannot be used in
#' conjunction with \code{selfing} or \code{mgs}.
#'
#' @param mgs
#'
#' If \code{TRUE} (default \code{FALSE}), the third identity in the
#' pedigree is the male parent of the female parent (maternal
#' grand-sire) rather than the female parent.
#'
#' @param mv
#'
#' A character vector of missing value indicators; elements of
#' \code{pedigree} that exactly match any of the members of \code{mv}
#' are treated as missing.
#'
#' @param psort
#'
#' If \code{TRUE} (default \code{FALSE}), the pedigree data frame is
#' returned in founder order after any insertions and permutations.
#'
#' @return
#'
#' A three-column matrix with class \code{ginv} holding the lower
#' triangle of the inverse relationship matrix in sparse form. The
#' first 2 columns are the \emph{row} and \emph{column} indices,
#' respectively, and the third column holds the inverse matrix element
#' itself. Sort order is columns within rows, that is, the lower
#' triangle row-wise. This matrix has attributes:
#'
#' \describe{
#'
#' \item{\code{rowNames}}{A character vector of identifiers for the
#' rows of the matrix.}
#'
#' \item{\code{inbreeding}}{A numeric vector containing the
#' inbreeding coefficient for each individual, calculated as
#' \code{diag(A-I)}.}
#'
#' \item{\code{geneticGroups}}{A numeric vector of length 2 containing
#' the \code{groups} and \code{groupOffset} arguments.}
#'
#' \item{\code{logdet}}{The log determinant.}
#'
#' }
#'
#' @examples
#'
#' \dontrun{
#'
#' # Simple pedigree
#'
#' ped <- data.frame(me = c(1,2,3,4,5,6,7,8,9,10),
#'                   dad = c(0,0,0,1,1,2,4,5,7,9),
#'                   mum = c(0,0,0,1,1,2,6,6,8,9))
#' p.ai <- ainv(ped)
#'
#' # Known filial generation
#'
#' pdfg <- data.frame(me = c(1,2,3,4,5,6,7),
#'                    dad = c(0,0,1,1,1,1,1),
#'                    mum = c(0,0,0,2,2,2,2),
#'                    fgen = c(0.8,0.0,2.0,0.0,2.0,3.0,0.0))
#' pdfg.ai <- ainv(pdfg,fgen=list('fgen',0.4))
#' pdfg.mat <- sparse2mat(pdfg.ai)
#' zapsmall(solve(pdfg.mat))
#' zapsmall(cbind(pdfg.a$inbreeding,diag(pdfg.mat)))
#' }
#'
#' @references
#' % bibentry: Dutkowski:2001
#' % bibentry: Fernando:1990
#' % bibentry: Meuwissen:1992
#'
ainv <- function(pedigree, fgen=list(character(0),0.01),
                 gender=character(0), groups=0, groupOffset=0.0,
                 selfing=NA, inBreed=NA, mgs=FALSE,
                 mv=c("NA","0","*"), psort=FALSE)
{
  method <- 0 ## meurwissen
  newped <- checkPedigree(pedigree, fgen, gender, mv)
  if(psort)
    return(newped)
  ##
  nan <- nrow(newped)
    if((!is.na(selfing)) & (!is.na(inBreed)))
    stop("Cannot specify both partial selfing and inbreeding coefficient\n")
  if(is.na(selfing))
    selfing <- 0.0
  ibl <- 1
  if(is.na(inBreed)) {
    inBreed <- 0.0
    ibl <- 0
  }

  fgsx <- rep(0.0,nan)
  fmode <- FALSE
  xlink <- FALSE
  if(ncol(newped) == 4) {
    fgsx <- newped[,4]
    if(attr(newped, "col4Type") == "fgen")
      fmode <- TRUE
    else if(attr(newped, "col4Type") == "gender")
      xlink <- TRUE
  }

  ## Convert to pmat for g5vaig
  lvls <- attr(newped,"rowNames")
  idx <- vector(mode="numeric",length=nan)
  pmat <- matrix(0,nrow=nan,ncol=2)
  vped <- matrix(match(c(newped[,1],newped[,2],newped[,3]),lvls,nomatch=0),
                 ncol=3,byrow=FALSE)

  rowNo <- rep(-1,nan)
  rowNo[vped[,1]] <- seq(1,nan)
  ## Can't happen??
  if(any(rowNo == -1))
    stop("Internal error: pedigree mismatch")
  for(col in 1:2) {
    which <- (vped[,col+1] > 0)
    pmat[which,col] <- rowNo[vped[which,col+1]]
  }
  idx[rowNo] <- seq(1,nan)
  pmat <- t(pmat)

  if(mgs)
    mgs <- 1
  else
    mgs <- 0
  debug <- 0

  params <- list(fmode=as.numeric(fmode), xlink=as.numeric(xlink), fgsx=fgsx,
                 pmat=pmat, nan=nan, groups=groups, groupOffeset=groupOffset,
                 method=method, selfing=selfing, ibl=ibl, mgs=mgs, debug=debug)
  ginv <- .Call("ainverse", params)

  dimnames(ginv)[[2]] <- c("Row","Column","Ainverse")
  inbreeding <- attr(ginv,"inbreeding")
  names(inbreeding) <- lvls
  attr(ginv,"inbreeding") <- inbreeding
  attr(ginv,"rowNames") <- lvls
  attr(ginv,"geneticGroups") <- c(groups, groupOffset)
  class(ginv) <- c(class(ginv),"ginv")

  return(ginv)
}
#' Simulated numerator relationship matrices.
#'
#' Calculates components of the genetic covariance matrix for a given
#' population defined in an ancestral tree.
#'
#' Initially, each founder individual starts with two distinct (and
#' unique) alleles.  For each simulation, the pedigree is traversed
#' and two genes are sampled, one from each parent, and assigned to
#' each individual. If \eqn{f > 0}, the assigned genes are sampled
#' with replacement \eqn{f} times, resulting in homozygosity at the
#' locus if \eqn{f} is sufficiently large. The genes of all
#' individuals are examined pairwise (\eqn{i,j}), and counts of events
#' \eqn{d_{ijs}} contributing to identity state S (see \cite{Lynch and
#' Walsh, 1998}, Chapter 7) are accumulated. If the genes for
#' individual \eqn{i} are identical, then the count of IBD events
#' within the \eqn{i^{th}} individual, \eqn{F_i}, is incremented.
#'
#' The procedure is repeated \code{nsim} (\eqn{N}) times, and the pairwise
#' coefficients \eqn{D_{ijs}} for state S estimated as \eqn{D_{ijs} =
#' d_{ijs}/N}, for \eqn{s=1:9}; the inbreeding coefficients are estimated
#' as \eqn{F_i = F_i/N}. The elements of the various relationship matrices
#' are estimated using the (pairwise) coefficients of the contributing
#' identity states as defined in \cite{Lynch and Walsh (1998)}.
#'
#' @param pedigree
#'
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively. The row
#' giving the pedigree of an individual must appear before any row
#' where that individual appears as a parent. Founders use 0 (zero) or
#' \code{NA} in the parental columns.
#'
#' @param fgen
#'
#' An optional list of length 2 where \code{fgen[[1]]} is a character
#' string naming the column in \code{pedigree} that contains the level
#' of selfing (\eqn{f}) or the level of inbreeding \eqn{b}) of an
#' individual. In \code{pedigree[,fgen[[1]]]}, 0 indicates a simple
#' cross, 1 indicates selfed once, 2 indicates selfed twice, etc. A
#' value between 0 and 1 for a base individual is taken as its
#' inbreeding value. If the pedigree has implicit individuals (they
#' appear as parents but not as individuals), they will be assumed
#' base non-inbred individuals unless their inbreeding level is set
#' with \code{fgen[[2]]} where \code{0<fgen[[2]]<1} is the inbreeding
#' level of such individuals. For use in the simulation, \eqn{f} is
#' determined from \eqn{b} using \eqn{b=1-(0.5)^f} if
#' \code{0<fgen[[2]]<1}.
#'
#' @param nsim
#'
#' The number of traversals of the ancestral tree; default 2000.
#'
#' @param rm
#'
#' A character vector naming the relationship matrices to be returned;
#' the default is to return all simulated matrices. Setting
#' \code{rm=c("A","D")} will return both the additive and dominance
#' relationship matrices, for example; see \cite{Lynch and Walsh
#' (1998)} for details on the components of the genetic covariance
#' matrix.
#'
#' @param form
#'
#' if \code{"full"}, return the complete relationship matrices,
#' otherwise if \code{form="upper"} the upper triangle is returned.
#'
#' @param mv
#'
#' Missing value indicator; elements of \code{pedigree} that exactly
#' match any element of \code{mv} are treated as missing.
#'
#' @references
#'
#' % bibentry: Lynch:1998
#' % bibentry: Verbyla:2011
#'
#'
nrm <- function(pedigree,fgen=list(character(0),0.0),
                nsim=2000, rm=c('A','D','C','Ct','Dh','Di','E'),
                form=c('full','upper'), mv=c("NA","0","*"))
{

  ## which ones to return
  rm.vec <- c('A','D','C','Ct','Dh','Di','E')
  which.rm <- match(rm, rm.vec)
  if(any(is.na(which.rm)))
    stop("'rm' must be one or more of 'A','D','C','Ct','Dh','Di','E'")
  form=match.arg(form)

  newped <- asreml::chkPed(pedigree, fgen, character(0), mv)
  nan <- nrow(newped)
  lvls <- attr(newped,"rowNames")
  ## Call to monte
  for(i in 1:3)
    newped[,i] <- match(newped[,i],lvls,nomatch=0)

  newped <- lapply(newped,function(x) {
  storage.mode(x) <- "double"
  x})
  storage.mode(nsim) <- "integer"

  x <- .Call("monte",newped,nsim)
  ## Return vectors are UPPER triangle row-wise
  out <- vector(mode="list",length=length(rm))
  names(out) <- rm
  if(form=="full") {
    for(i in which.rm) {
      out[[rm.vec[i]]] <- utri2mat(x$RM[,i])
      dimnames(out[[rm.vec[i]]]) <- list(lvls,lvls)
    }
  }
  else {
    for(i in which.rm) {
      out[[rm.vec[i]]] <- x$RM[,i]
      attr(out[[rm.vec[i]]],"rowNames") <- lvls
    }
  }
  out$Fi <- x$Fi
  out
}
#' Remove individuals from an inverse relationship matrix.
#'
#' Deletes individuals from an inverse numerator relationship matrix
#' by sweeping out the respective rows.
#'
#' @param ginv
#'
#' A matrix containing an inverse numerator relationship matrix in
#' three column co-ordinate form with a \code{rowNames} attribute.
#'
#' @param keep
#'
#' A character vector naming the individuals to retain.
#'
#' @return
#'
#' A three column sparse coordinate form matrix in row major order
#' with attribute \code{"INVERSE"} set to \code{TRUE}.
#'
AIsweep <- function(ginv, keep)
{
  ## ginv may be a data frame typically from asreml.Ainverse()
  ##      coerce to a matrix
  ## or a matrix from ainv() or ainverse().

  if(missing(keep))
    stop("Nothing to do ('keep' not set).")
  if(is.null(rn <- attr(ginv,"rowNames")))
    stop("Missing 'rowNames' attribute")
  nan <- length(rn)
  if(is.factor(keep))
    keep <- levels(keep)
  ## coerce to character
  if(!is.character(keep))
    keep <- as.character(keep)
  keep <- match(keep,rn)
  if(any(is.na(keep)))
    stop("There are individuals in 'keep' that are not in 'rowNames'.")

  if(is.data.frame(ginv))
    ginv <- as.matrix(ginv)
  dn <- dimnames(ginv)[[2]]

  ## sort ginv such that those to be omitted are at the end
  what <- c(keep,seq(1,nan)[-keep])
  ginv <- sortGinv(ginv, what, old=TRUE)
  what <- seq(1,length(keep))
  save.rn <- attr(ginv, "rowNames")
  ginv <- cbind(ginv, as.numeric(is.element(ginv[,1],what)))
  attr(ginv, "rowNames") <- save.rn
  storage.mode(ginv) <- "double"

  ## comes back as a matrix
  gi <- .Call("aisweep2", ginv, PACKAGE="pedicure")
  dimnames(gi)[[2]] <- dn
  attr(gi,"rowNames") <- rn[keep]
  attr(gi, "INVERSE") <- TRUE
  gi
}
#'
#' Count descendents.
#'
#' Counts the number of offspring for each individual in a pedigree.
#'
#' @param ped
#'
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively. The row
#' giving the pedigree of an individual must appear before any row
#' where that individual appears as a parent. Founders use 0 (zero) or
#' \code{NA} in the parental columns.
#'
#' @return
#'
#' A numeric vector of length \code{nrow(ped)} containing the number
#' of offspring for each individual.
#'
offspring <- function(ped)
  {
    mmd <- data.frame(id = 1:nrow(ped),
                      dam = match(ped[,2],ped[,1],nomatch = 0),
                      sire = match(ped[,3],ped[,1],nomatch = 0))
    return(pedigree::countOff(mmd))
  }
#' Count generation number.
#'
#' Counts the generation number for each individual in a pedigreee.
#'
#' @param ped
#'
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively. The row
#' giving the pedigree of an individual must appear before any row
#' where that individual appears as a parent. Founders use 0 (zero) or
#' \code{NA} in the parental columns.
#'
#' @return
#'
#' A numeric vector of length \code{nrow(ped)} containing the
#' generation number for each individual.
#'
count_gen <- function(ped)
{
  mmd <- data.frame(id = 1:nrow(ped),
                    dam = match(ped[,2],ped[,1],nomatch = 0),
                    sire = match(ped[,3],ped[,1],nomatch = 0))
  return(pedigree::countGen(mmd))
}
#' Trim a pedigree.
#'
#' Trim a pedigree to a given data frame; non-informative individuals
#' without data are removed.
#'
#' Call \code{chkPed} (library \code{asreml}) to add any missing
#' founders, resolve \code{fgen}, and reorder the pedigree (if
#' necessary) prior to pruning.
#'
#' @param ped
#'
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively. The row
#' giving the pedigree of an individual must appear before any row
#' where that individual appears as a parent. Founders use 0 (zero) or
#' \code{NA} in the parental columns.
#'
#' @param data
#'
#' A data frame with a component named \code{names(ped)[1]} containing
#' the identities of the individuals with data.
#'
#' @param gen
#'
#' The number of generations to keep before those with data. The default
#' is \code{max(count_gen(mmd))}.
#'
#' @return
#'
#' A data frame with the same structure as \code{ped}, with the pruned
#' individuals removed.
#'
trim <- function(ped, data, gen = NULL)
{

  if(is.na(which <- match(names(ped)[1],names(data))))
    stop(paste("Cannot find",names(ped)[1],"in data"))
  if(any(is.na(match(as.character(data[[which]]), as.character(ped[,1])))))
    warning("There are individuals in 'data' that are absent in 'ped'")

  data <- as.numeric(!is.na(match(as.character(ped[,1]),
                                  as.character(data[[which]]))))
  mmd <- data.frame(id = 1:nrow(ped),
                    dam = match(ped[,2],ped[,1],nomatch = 0),
                    sire = match(ped[,3],ped[,1],nomatch = 0))

  if(is.null(gen)) {
    gen <- max(count_gen(mmd))
  }
  what <- pedigree::trimPed(mmd, data, gen)
  return(ped[what,])
}
#' Extract parental records from a (pruned) pedigree.
#'
#' Trim a pedigree to the set of unique identifiers in the parental
#' columns.
#'
#' @param ped
#'
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively; typically
#' the result of a previous call to \code{\link{trim}}.
#'
#' @return
#'
#' A data frame with the same structure as \code{ped}, with the pruned
#' individuals removed.
#'
trim.par <- function(ped)
{
  data <- as.numeric(!is.na(match(as.character(ped[,1]),
                                  unique(c(as.character(ped[,2]),as.character(ped[,3]))))))
  mmd <- data.frame(id = 1:nrow(ped),
                    dam = match(ped[,2],ped[,1],nomatch = 0),
                    sire = match(ped[,3],ped[,1],nomatch = 0))

  gen <- as.integer(max(count_gen(ped)))
  what <- pedigree::trimPed(mmd, data, gen)
  return(ped[what,])
}
#' Dot file reporesenting a pedigree tree.
#'
#' Exports the relationship structure in a relationship matrix as a
#' directed graph in a dot file suitable for plotting.
#'
#' The matrix \code{A} can be decomposed as \code{A=LDU}, where
#' \code{U=t(L)} and \code{inv(L)} represents the directed graph of
#' additive genetic relationships in the pedigree. This matrix can be
#' converted to a graph object in dot file format suitable for
#' plotting by \code{Graphviz}, say.
#'
#' The function calls \code{zapsmall} to round nuisance values to
#' zero; this can be controlled through the \code{"digits"} component
#' of \code{options()}
#'
#' @param A
#'
#' The numerator relationship matrix.
#'
#' @param width
#'
#' Width of the graph node character identifiers. If \code{NA} (the
#' default) the complete identity names are written to the dot file.
#'
#' @param dotfile
#'
#' The output dot file name, default is \code{"dot.dot"}.
#'
#' @return
#'
#' A two column matrix of the graph node relationships.
#'
#' @section Side effects:
#'
#' The \code{dotfile} is written to the working directory.
#'
A.dot <- function(A, width=NA, dotfile="A.dot")
{
  id <- dimnames(A)[[1]]
  if(length(id) == 0)
    id <- as.character(seq(1,nrow(A)))
  if(!is.na(width))
    id <- strtrim(id, width)

  L <- zapsmall(ldu(A)$L)
  tLinv <- t(solve(L))
  which <- (!lower.tri(tLinv,diag=TRUE) & tLinv!=0)
  x <- cbind(row(tLinv)[which],col(tLinv)[which])
  cat("digraph grafname {", "\n", file=dotfile)
  cat(paste(id[x[, 1]], id[x[,2]], sep=" -> "), sep="\n", file=dotfile, append=TRUE)
  cat("}", "\n", file=dotfile, append=TRUE)
  x
}
#' Relationship matrix
#'
#' Calculate the numerator relationship matrix from a pedigree.
#'
#' @param ped
#'
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively. Founders
#' use 0 (zero) or \code{NA} in the parental columns.
#'
#' @param keep
#'
#' A logical vector identifying the rows of \code{ped} to retain.
#'
#' @section Side-effects:
#' A file \code{"A.txt"} is created
#'
#' @return
#'
#' A \code{matrix} object containing the numerator
#' relationship matrix.
#'
amat <- function(ped, keep=rep(TRUE, nrow(ped)))
{
  if(!is.logical(keep))
    stop("keep should be a logical")
  else if(length(keep) != nrow(ped))
    stop("Length of 'keep' should be", nrow(ped))
  #keep <- whi
  (keep)
  mmd <- data.frame(id = 1:nrow(ped),
                    dam = match(ped[,2],ped[,1],nomatch = 0),
                    sire = match(ped[,3],ped[,1],nomatch = 0))
  pedigree::makeA(mmd, keep)
  A <- as.matrix(read.table("A.txt"))
  return(A)
}
#' Dot file representing a pedigree tree.
#'
#' Exports the relationship structure in a pedigree as a directed
#' graph in a dot file suitable for plotting.
#'
#' The resulting dot file can be edited prior to rendering with
#' \emph{Graphviz}, say, or conversion to a graphics file format with
#' the \emph{dot} application (part of the \emph{Graphviz} package). The
#' default draws founder individuals in rectangles, offspring in
#' ellipses, and offspring with a single parent in circles.
#'
#' @param ped
#'
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively. Founders
#' use 0 (zero) or \code{NA} in the parental columns.
#'
#' @param keep
#'
#' A logical vector identifying the rows of \code{ped} to retain. The
#' default is \code{rep(TRUE, nrow(ped))}.
#'
#' @param dotfile
#'
#' The output dot file primary name, default is \code{"ped"}; the
#' suffix \code{".dot"} is appended to the file name. Also used to
#' label the pedigree tree.
#'
#' @param numeric
#'
#' Convert \code{ped} to integer values for graphing. The default is
#' \code{TRUE}.
#'
#' @param url
#'
#' If not \code{""} (the default) then \code{url} is linked to the resulting graph.
#'
#' @param height
#'
#' Node height.
#'
#' @param width
#'
#' Node width
#'
#' @param rotate
#'
#' If \code{rotate=90} landscape mode is selected; the default is
#' \code{0}.
#'
#' @section Side effects:
#'
#' The \code{dotfile} is written to the working directory using
#' \code{sink()}.
#'
ped.dot <- function (ped, keep=rep(TRUE, nrow(ped)), dotfile="ped", numeric=TRUE,
                     url = "", height = 0.5, width = 0.75, rotate = 0)
{
  ped <- ped[keep,]
  ssize <- dim(ped)[1]

  if(numeric) {
    me <- 1:ssize
    mum <- match(ped[,2], ped[,1], nomatch=NA)
    dad <- match(ped[,3], ped[,1], nomatch=NA)
  }
  else {
    me <- ped[, 1]
    mum <- ped[, 2]
    mum[mum == 0] <- NA
    dad <- ped[, 3]
    dad[dad == 0] <- NA
  }

  ## shape: rectangle for founders, circle for one unknown parent, elipse for progeny

  shape <- rep("ellipse", nrow(ped))
  founders <- as.numeric(is.na(mum)) + as.numeric(is.na(dad))
  shape[founders == 2] <- "box"
  shape[founders == 1] <- "circle"

  cat(paste("[", dotfile, "]", sep = ""))
  sink(paste(dotfile, ".dot", sep = ""))
  cat(paste("digraph ped_", dotfile, sep = ""), "{\n")
  cat("fontname=Helvetica;\n")
  cat("label=\"Pedigree:", dotfile, "\" ;\n")
  cat(paste("rotate=", rotate, sep = ""), " ;\n")
  if (url != "")
    cat(paste("URL=\"", url, "\"", sep = ""), " ;\n")

  cat(paste("\"",me,"\" [shape=", shape,
            ",height=", height,
            ",width=", width,
            ",fontname = Helvetica",
            "] ;",sep=""), sep="\n")

  cat(paste("\"",dad,"\"", " -> ", "\"",me,"\"", ";",sep="")[!is.na(dad)], sep="\n")
  cat(paste("\"",mum,"\"", " -> ", "\"",me,"\"", ";",sep="")[!is.na(mum)], sep="\n")
  cat("}\n")
  sink()

  cat("\n")
  invisible()
}
utri2mat <- function(vec, diag=TRUE, cor=TRUE)
{
  ## vec is assumed upper triangle row-wise
  ## diag:
  ##   logigal
  ##     TRUE: us type matrix
  ##     FALSE corg type
  ##   numeric
  ##     diag gives elements of vec that form the diagonal

  makeMat <- function(v,m,diagonal) {
    mat <- matrix(nrow=m,ncol=m)
    mat[lower.tri(mat, diagonal)] <- v
    mat[!lower.tri(mat, TRUE)] <- t(mat)[!lower.tri(mat,TRUE)]
    mat
  }
  if(is.logical(diag)) {
    n <- length(vec)
    dii <- numeric(0)
    if(diag)
      m <- (sqrt(8*n+1)-1)/2
    else {
      m <- (sqrt(8*n+1)+1)/2
      dii <- rep(ifelse(cor,1,0),m)
    }
    mat <- makeMat(vec,m,diag)
  }
  else if(is.numeric(diag)) {
    nd <- length(diag)
    n <- length(vec)-nd
    dii <- {if(nd==1)
              rep(vec[diag],n)  ## corgv type
            else
              vec[diag]}        ## corgh type
    m <- (sqrt(8*n+1)+1)/2
    mat <- makeMat(vec[1:n],m,FALSE)
  }
  if(length(dii))
    diag(mat) <- dii

  return(mat)
}
#' Convert a sparse matrix.
#'
#' Convert a sparse matrix held in three-column coordinate form to a
#' dense matrix.
#'
#' The attribute \code{rowNames} of \code{x} is preserved and returned
#' as the \code{dimnames} attribute of the returned matrix.
#'
#' @param x
#'
#' A matrix (or data frame) holding the lower triangle of a sparse
#' symmetric matrix in coordinate form in row major order; typically
#' the result of a call to \code{ainv}.
#'
#' @return
#'
#' A dense matrix containing \code{x} padded with zeros.
#'
sparse2mat <- function(x)
{
  ## sparse form is 3 col matrix or dataframe
  ## row,col,value
  ## mainly for ginv from Ainverse

  rn <- attr(x,"rowNames")

  nrow <- max(x[,1])
  ncol <- max(x[,2])
  y <- rep(0,nrow*ncol)
  y[(x[,2]-1)*nrow+x[,1]] <- x[,3]
  y[(x[,1]-1)*nrow+x[,2]] <- x[,3]
  a <- matrix(y,nrow=nrow,ncol=ncol,byrow=FALSE)
  dimnames(a) <- list(rn,rn)
  a
}
#' Generalised Cholesky decomposition.
#'
#' Generalised Cholesky decomposition of a numerator relationship
#' matrix.
#'
#' Determines the modified Cholesky decomposition of a
#' positive-definite matrix \eqn{A} given by \eqn{A=LDL^T}, where
#' \eqn{L} is lower triangular with ones on the diagonal, and \eqn{D}
#' is a diagonal matrix. The function uses the \code{MCHOL} function
#' in \cite{McLeod and Holanda Sales, 1983}.
#'
#' @param A
#'
#' The numerator relationship matrix.
#'
#' @return
#'
#' A list with the following components:
#' \describe{
#'  \item{L}{A lower triangular matrix.}
#'  \item{D}{A diagonal matrix.}
#' }
#'
#' @references
#'
#' % bibentry: McLeod:1983
#'
ldu <- function(A)
{
  ans <- as.double(A[!lower.tri(A,diag=FALSE)])
  n <- as.integer(nrow(A))
  det <- as.double(0.0)
  ifault <- as.integer(0)
  out <- .C("mchol",ans=ans,n,det=det,ifault=ifault)
  if(out$ifault == 2)
    warning("A is not positive-definite")
  A[!lower.tri(A,diag=FALSE)] <- out$ans
  A[row(A)>col(A)] <- t(A)[lower.tri(A)]
  A[row(A)<col(A)] <- 0
  D <- diag(A)
  diag(A) <- 1.0
  list(L=A,D=D)
}
#' Convert a dense matrix.
#'
#' Convert a dense matrix to a sparse matrix in three-column
#' coordinate form.
#'
#' The \code{dimnames(X)[[1]]} attribute of \code{x} is preserved and returned
#' as the \code{rowNames} attribute of the returned matrix.
#'
#' @param x
#'
#' A dense matrix padded with zeros where needed.
#'
#' @param rowNames
#'
#' A character vector set as the "rowNames" attribute of the returned
#' sparse matrix. The default is the \code{dimnames} attribute of the
#' leading dimension of \code{x}.
#'
#' @param tol
#' 
#' Matrix elements whose absolute value is less than tol are considered zero.
#' The default is 1e-9.
#' 
#' @return
#'
#' A matrix holding the lower triangle of a sparse symmetric matrix in
#' coordinate form in row major order.
#'
mat2sparse <- function(x, rowNames=dimnames(x)[[1]], tol=1e-9)
{
  which <- (abs(x) > tol & lower.tri(x, diag=TRUE))
  df <- as.matrix(data.frame(row = t(row(x))[t(which)],
                   col = t(col(x))[t(which)],
                   val = t(x)[t(which)]))
  if(is.null(rowNames))
    rowNames <- as.character(1:nrow(x))
  attr(df,"rowNames") <- rowNames
  df
}
sortGinv <- function(ginv, ord, old=FALSE) {
  if(old) {
    tmp <- sparse2mat(ginv)[ord,ord] # convert to dense matrix
    dimnames(tmp) <- list(attr(ginv,"rowNames")[ord],attr(ginv,"rowNames")[ord])
    return(mat2sparse(tmp))
  }
  else { #any better? - saves space
    if(!requireNamespace("Matrix", quietly=TRUE)) {
    stop("Requires package 'Matrix'")
  }
    tmp <- Matrix::sparseMatrix(i=ginv[,1], j=ginv[,2], x=ginv[,3], symmetric=TRUE)[ord,ord]
    tmp <- as(tmp,"dgTMatrix")
    tmp <- cbind(tmp@i+1, tmp@j+1, tmp@x) ## zero based
    tmp <- tmp[(tmp[,1] >= tmp[,2]),] #lower tri
    tmp[base::order(as.numeric(tmp[,1]), as.numeric(tmp[,2])),]
    attr(tmp, "rowNames") <- attr(ginv,"rowNames")[ord]
    return(tmp)
  }
}
#' Realized relationship matrix.
#'
#' Calculates an IBS (possibly) additive relationship matrix.
#'
#' The realized relationship matrix is computed following two
#' pre-processing steps: 1) the matrix \code{x} is scanned for
#' consistency against the arguments \code{min.af} and \code{max.mv},
#' and 2) any missing values are imputed using \code{na.method}.
#' Monomorphic and completely missing markers are removed in step (1).
#' If \code{center=TRUE}, the
#' imputed matrix \eqn{W} is optionally centred, and if
#' \code{scale=TRUE} \eqn{W} is scaled by the reciprocal of the average marker
#' variance (\cite{VanRaden 2008}) before computing
#' \eqn{WW^\prime}. The imputation methods \code{svd} and \code{knn}
#' replace the missing values in \code{x} with their column means
#' before proceeding.
#'
#' The \code{svd} method (\cite{Troyanskaya etal 2001}) uses the
#' first \code{ev} terms of the SVD of \eqn{W} to replace the
#' previously missing values with their regression predictions from
#' the SVD. The procedure iterates until convergence is achieved or
#' \code{maxit} is exceeded. Convergence is declared when \eqn{|RSS_0
#' - RSS_1|/RSS_1 < 0.02} (\cite{Rutkoski 2013}), where \eqn{RSS} is
#' the residual sum of squares between the non-missing values and
#' their predictions from the SVD, and \eqn{RSS_0} and \eqn{RSS_1} are
#' the RSS values from successive iterations.
#'
#'  The \code{knn} procedure (\cite{Troyanskaya etal 2001}) orders the
#' (marker) neighbours of the \eqn{j^{th}} marker according to
#' euclidean distance, replacing the \eqn{i^{th}} missing value with
#' the average of the non-missing values for genotype \eqn{i} from the
#' \eqn{knn} nearest neighbours of marker \eqn{j}.
#'
#' If \code{imputed} is \code{TRUE}, the imputed marker matrix is
#' returned with monomorphic and completely missing markers removed,
#' as well as those screened out by \code{min.af} and \code{max.mv}.
#'
#' If \code{scale} is \code{TRUE}, the reciprocal of the average marker
#' variance is returned in an attribute \code{scale} of the returned
#' matrix.
#'
#' @param x
#'
#' An \ifelse{html}{\out{n x m}}{\eqn{n \times m}} matrix of biallelic marker scores for
#' \eqn{n} individuals and \eqn{m} markers. The scores are assumed
#' coded as \{-1, 0, 1\} with missing values coded as \code{NA}.
#'
#' @param center
#'
#' If \code{TRUE}, the columns of \code{x} (or the imputed matrix
#' \eqn{W}) are centred.
#'
#' @param scale
#'
#' If \code{TRUE} (the default), the realized association matrix
#' \eqn{xx'} (or \eqn{WW'}) is scaled by the reciprocal of
#' the average marker variance (\eqn{v}); that is, \eqn{xx'/v}.
#' If \code{scale} is numeric then the
#' realized association matrix is \eqn{xx' * scale}.
#'
#' @param na.method
#'
#' \describe{
#'
#' \item{mean}{missing values within each marker column are replaced
#'	by the column mean;}
#'
#' \item{knn}{a missing value for a marker is replaced by the
#'	mean of its \code{knn} nearest neighbours with non-missing data
#'	for that individual;}
#'
#' \item{svd}{iteratively impute missing values in a predicted
#'	(current) \eqn{W_1} using the \code{ev} most significant eigenvalues
#'	from a singular value decomposition of the previous \eqn{W_0};}
#'
#' \item{none}{no imputation.}
#' }
#'
#' @param min.af
#'
#' Minimum allele frequency below which markers are removed;
#' monomorphic markers are automatically removed. If not set, the default is
#' \code{1/(2*nrow(x))}.
#'
#' @param max.mv
#'
#' Maximum proportion of missing values allowed for a marker;
#' completely missing markers are removed. If not set, the default is
#' \code{1-1/(2*nrow(x))}.
#'
#' @param ev
#'
#' If \code{na.method = svd}, the number of significant eigenvectors
#' to use when imputing missing values (the default is 10).
#'
#' @param knn
#'
#' If \code{na.method = knn}, the cluster size of similar markers to
#' use when imputing missing values (the default is 10).
#'
#' @param maxit
#'
#' If \code{na.method = svd}, the maximum number of iterations used to
#' impute missing values (the default is 10). The iteration sequence may terminate before
#' \code{maxit} if the convergence criterion is met (see details).
#'
#' @param imputed
#'
#' If \code{TRUE} (the default is \code{FALSE}), a list of length two
#' is returned with components \code{rrm} and \code{inputed}
#' containing the realized relationship and the imputed marker
#' matrices, respectively. If \code{FALSE}, only the relationship matrix
#' is returned.
#'
#' @return
#'
#' A list with components \code{rrm} and \code{imputed}, or just the
#' relationship matrix object depending on the value
#' of the \code{imputed} argument. For \code{method=svd}, the
#' convergence sequence (RSS) is included as attribute \code{rss} of
#' the returned object.
#'
#' @references
#'
#' % bibentry: Rutkoski:2013
#' % bibentry: Troyanskaya:2001
#' % bibentry: VanRaden:2008
#'
#' @examples
#'
#' \dontrun{
#' ## Wheat data from Rutkoski (2013) with 70\% missing
#' data(WW_VersionNA70_rep1)
#' ## the imputed marker matrix
#' W <-  mmat(WW_VersionNA70_rep1,na.method="knn",knn=4, imputed=TRUE)
#' ## the realized relationship matrix
#' WW <-  mmat(WW_VersionNA70_rep1,na.method="knn",knn=4)
#' }
#'
mmat <- function(x, center=FALSE, scale=TRUE, na.method=c("mean","svd","knn","none"),
                 min.af=NULL, max.mv = NULL,
                 ev=10, knn=10, maxit=10, imputed=FALSE) {
  ## x is biallelic marker matrix coded -1,0,1,NA
  ## if(center && scale) form MM' from x using VanRaden, 2008.

  ## mmt:
  ## 1. impute missing values
  ## 2. form MM'

  if(!is.matrix(x))
    stop(substitute(x)," must be a numeric matrix.")
  if(min(x, na.rm=TRUE) < -1 || max(x, na.rm=TRUE) > 1)
    stop("'x' must be coded {-1, 0, 1, NA}")
  if(length(dimnames(x)) == 0) {
    warning(substitute(x)," has no 'dimnames' attribute")
    dimnames(x) <- list(seq(1,nrow(x)),seq(1,ncol(x)))
  }
  if(missing(na.method))
    na.method = "mean"

  if(is.logical(scale)) {
    if(scale)
      scale = 1e37 # flag van raden
    else
      scale = 1.0e0
  }
  n <- nrow(x)
  if (is.null(min.af)) {
    min.af <- 1/(2 * n)
    message("'min.af' set to ", min.af)
  }
  if (is.null(max.mv)) {
    max.mv <- 1 - 1/(2 * n)
    message("'max.mv' set to ",max.mv)
  }
  options <- list(centre=center,
                  scale=scale,
                  na.method=na.method,
                  min.af=min.af,
                  max.mv=max.mv,
                  ev=ev,
                  knn=knn,
                  maxit=maxit,
                  imputed=imputed)
  out <- .Call("mmt", x, options)
  out
}
#' Eigenvalues.
#'
#' Calculate the eigenvalues of a symmetric matrix.
#'
#' @param x
#'
#' A real symmetric matrix.
#'
#' @param k
#'
#' Calculate the k largest eigenvalues and eigenvectors; the default
#' is \code{nrow(x)}.
#'
#' @param tol
#'
#' Error tolerance to which each eigenvalue is required. If
#' \code{tol=-1} (the default), it is set to a machine dependent
#' \emph{safe} minimum.
#'
#' @return
#'
#' \describe{
#'
#' \item{V}{An \eqn{n \times k} numeric matrix of eigenvectors.}
#' \item{lambda}{The vector of eigenvalues of length k.}.
#' }
#'
mmat.ev <- function(x, k = nrow(x), tol=-1.0) {
  ## x is MM'
  ## get eigenvalues
  out <- .Call("mmev", x, k, tol)
  out
}
#' Singular value decomposition.
#'
#' Calculate the singular value decomposition of a rectangular matrix.
#'
#' @param x
#'
#' A real \eqn{m \times n} matrix.
#'
#' @return
#'
#' A list with three components:
#'
#' \describe{
#'
#' \item{u}{An \eqn{m \times m} matrix of left singular vectors.}
#' \item{v}{An \eqn{n \times n} transposed matrix of right singular vectors.}
#' \item{d}{The \code{min(m,n)} vector of singular values.}
#' }
#'
mmat.svd <- function(x) {
  ## x is MM'
  ## get svd

  out <- .Call("mmsvd", x)
  out
}
#' Recode a marker matrix.
#'
#' Recode the matrix of molecular marker scores using a simple
#' mapping.
#'
#' @param x
#'
#' An \eqn{n\times m} matrix of scores for \eqn{m} molecular markers
#' on \eqn{n} individuals; missing values are allowed and retained.
#'
#' @param map
#'
#' A two column matrix where \code{map[,1]} are the existing scores in
#' \code{x}, and \code{map[,2]} are the new scores. The default
#' assumes scores \code{(0,1,2)} are to be recoded as \code{(-1,0,1)}.
#'
#' @param char.convert
#'
#' If \code{TRUE} (the default) and \code{x} (and therefore
#' \code{map}) is of type \code{character}, \code{x} will attempted to
#' be coerced to \code{numeric} after recoding.
#'
#' @param transpose
#'
#' If \code{TRUE} (the default is \code{FALSE}), transpose \code{x} if
#' rows are markers and columns are individuals.
#'
#' @return
#'
#' The recoded (and possibly transposed) matrix with missing values
#' preserved.
#'
recode.mm <- function(x, map=matrix(c(0,1,2,-1,0,1),ncol=2), char.convert=TRUE,
                      transpose=FALSE) {

  ## x is marker matrix
  ## map[,1] is original, map[,2] is recoded

  x <- as.matrix(x)
  if(is.null(dimnames(x)))
    stop(substitute(x)," must be able to be coerced to a matrix with a dimnames attribute")

  if(transpose)
    x <- t(x)

  centre <- function(y, map) {
    ## y is a (marker) column of x
    a <- y
    for(i in 1:nrow(map))
      y[which(a == map[i,1])] <- map[i,2]
    return(y)
  }
  if(is.character(x) && char.convert)
    return(type.convert(apply(x,2,centre,map), as.is=TRUE))
  else
    return(apply(x,2,centre,map))
}
#' Calculate Hamming distances for genetic markers.
#'
#' Calculate the Hamming distance for each marker pair from a matrix
#' of marker scores.
#'
#' For a given pair of markers, the Hamming distance is calculated as
#' the total number of mismatches scaled by n. If \code{na.match =
#' FALSE}, a missing value in either marker will count as a mismatch
#' and contribute to the Hamming distance summation.
#'
#' @param M
#'
#' An \code{n x m} matrix of marker scores for m markers on n individuals. The
#' matrix must be numeric and capable of being coerced to integer
#' values.
#'
#' @param na.match
#'
#' If \code{FALSE} (the default), \code{NA} in either marker counts as
#' a mismatch, that is \code{NAs} match nothing. If \code{TRUE},
#' \code{NA} in either marker counts as a match, that is \code{NAs}
#' match everything.
#'
#' @return
#'
#' A \code{dist} class numeric vector of distances being the strict
#' lower triangle of the distance matrix in column major order.
#'
hdist <- function(M, na.match = FALSE) {
  ## matrix of pairwise hamming distances
  ## M is a numeric matrix that can be coerced to int
  ## NA.match = FALSE NAs match nothing
  ##            TRUE NAs match anything

  if(!is.matrix(M))
    stop("'M' must be a matrix.")
  if(!is.numeric(M))
    stop("'M' must be numeric.")

  sim <- .Call("hamming", M, na.match)

  ## return a symmetric Matrix class object

  class(sim) <- c("dist", class(sim))
  return(sim)
}
#' Equivalence classes for (near) co-incident genetic markers.
#'
#' Identify redundant genetic genetic markers based on the Hamming
#' distance between marker pairs.
#'
#' For a given pair of markers, the Hamming distance is calculated as
#' the total number of mismatches scaled by n. If \code{na.match =
#' FALSE}, a missing value in either marker will count as a mismatch
#' and contribute to the Hamming distance summation.
#'
#' Beginning with column 1 of \code{M}, the \eqn{m(m-1)/2} pairs of
#' markers are scanned sequentially and markers allocated to
#' equivalence classes if the pairwise Hamming distance is less than
#' \code{threshold}. The algorithm uses the method of D. Eardley in
#' Section 8.6 of \cite{Press et al., 2002}. Pairwise Hamming
#' distances are computed on the fly from \code{M} to avoid storing
#' the distance matrix.
#'
#' @param M
#'
#' An \code{n x m} matrix of marker scores for m markers on n individuals. The
#' matrix must be numeric and capable of being coerced to integer
#' values.
#'
#' @param na.match
#'
#' If \code{FALSE} (the default), \code{NA} in either marker counts as
#' a mismatch, that is \code{NAs} match nothing. If \code{TRUE},
#' \code{NA} in either marker counts as a match, that is \code{NAs}
#' match everything.
#'
#' @param threshold
#'
#' Marker pairs whose Hamming distance is less than \code{threshold}
#' (the default is 0.05) are deemed equivalent and assigned to an
#' equivalence class.
#'
#' @return
#'
#' A list of length the number of equivalence classes, where each
#' component is a numeric vector of equivalent marker column
#' numbers. The names of the list are the root marker numbers of the
#' equivalence classes and can be coerced to \code{numeric}.
#'
#' @references
#'
#' % bibentry: Press:2002
#'
equiv.mm <- function(M, na.match = FALSE, threshold = 0.05) {
  ## matrix of pairwise hamming distances
  ## M is a numeric matrix that can be coerced to int
  ## NA.match = FALSE NAs match nothing
  ##            TRUE NAs match anything
  ## threshold: keep if hamming > threshold

  if(!is.matrix(M))
    stop("'M' must be a matrix.")
  if(!is.numeric(M))
    stop("'M' must be numeric.")

  eqv <- .Call("equivmkr", M, na.match, threshold)

  ## return a list of equivalence classes.

  eqcl <- lapply(unique(eqv), function(x, u) {which(u == x)}, eqv)
  names(eqcl) <- unlist(lapply(eqcl, function(x)x[1]))

  return(eqcl)
}
#' Linkage disequilibrium
#'
#' Calculate the LD coefficient \eqn{r} from genotype marker scores.
#'
#' This function accessess a modification to the EM method of
#' \cite{Excoffier and Slatkin, 1995}, and the approximation of
#' \cite{Rogers and Huff, 2009}. The approximation works with pairs of
#' biallelic loci, and the EM method has been modified to also deal
#' only with pairs of loci.
#'
#' The advantages of the approximate method are its ease of coding and
#' speed improvement over the EM optimisation method, see \cite{Rogers
#' and Huff, 2009} for details. This implementation uses a simple high
#' level multithreading strategy for both methods, so timing is
#' proportional to the number of cores. 
#'
#' The underlying compiled code is from 
#' \url{https://github.com/alanrogers/covld} with minor modifications
#' to capture printed results in \pkg{R}.
#'
#' @param x
#'
#' An \eqn{n \times m} matrix of biallelic marker scores for \eqn{n}
#' individuals and \eqn{m} markers. The scores are assumed coded as
#' \{-1, 0, 1\} representing genotypes \code{aa}, \code{aA} and
#' \code{AA}, respectively, with missing values coded as \code{NA}.
#'
#' @param method
#'
#' If method is \code{"rh"} (the default), \eqn{r} is calculated for
#' each pair of markers in \code{x} using the approximate method of
#' \cite{Rogers and Huff}, otherwise a modified EM algorithm is used.
#'
#' @param verbose
#'
#' If \code{TRUE} and \code{method="em"}, any irregularities in
#' calculating the haplotype frequencies are reported. Warning: this
#' may generate considerable console output and would be most useful
#' in debugging a limited set of markers. The default is \code{FALSE}.
#'
#' @return
#'
#' The \eqn{m(m+1)/2} LD \eqn{r_{ij}} values stored as a dense
#' symmetric \code{"dspMatrix"} class Matrix; the diagonal elements
#' are set to 1.
#'
#' @references
#'
#' % bibentry: Rogers:2009
#' % bibentry: Excoffier:1995
#'
ld <- function(x, method=c("rh", "em"), verbose = FALSE)
{
  if(!is.matrix(x))
    stop(substitute(x)," must be a numeric matrix.")
  if(min(x, na.rm=TRUE) < -1 || max(x, na.rm=TRUE) > 1)
    stop("'x' must be coded {-1, 0, 1, NA}")
  if(length(dimnames(x)) == 0) {
    warning(substitute(x)," has no 'dimnames' attribute")
    dimnames(x) <- list(seq(1,nrow(x)),seq(1,ncol(x)))
  }
  ## Check for imputed or other values
  uu <- unique(x)
  uu <- uu[!is.na(uu)]
  if(any(!is.element(uu, c(-1,0,1))))
    stop("Elements of 'x' not in the set (-1,0,1)")

  ## recode x as 0,1,2 since the underlying C uses these values as subscripts
  x <- recode.mm(x, map=matrix(c(-1,0,1,0,1,2),ncol=2))
  ## set missing to -1
  x[is.na(x)] <- -1
  n <- dim(x)[2]
  storage.mode(x) <- "integer"
  method <- match.arg(method)
  ld_stats <- .Call("ld_driver", x, method, verbose)
  ld_stats[is.nan(ld_stats)] <- NA
  return(
    new("dspMatrix",
        uplo="U",
        Dim=as.integer(c(n,n)),
        Dimnames=list(dimnames(x)[[2]], dimnames(x)[[2]]),
        x=ld_stats)
  )
}
wordsize <- function(n,k,with.replacement=FALSE) {
  ## k = nCr or nPr
  ## find r
  if(with.replacement) {
    prod <- 1
    r <- 0
    while(prod < k) {
      r <- r+1
      prod <- prod * (n-r+1)
    }

  }
  else {
    r <- 1
    while(choose(n,r) < k) {
      r <- r+1
    }
  }
  return(r)
}
aka <- function(k) {
  ## use letters
  ## no replacement (combn)
  r <- wordsize(length(letters), k)
  words <- apply(combn(sample(letters), r),2,function(x)paste(x,collapse=""))[1:k]
  return(words)
}
is.missing <- function(pcol,mv) {
  isna <- is.na(pcol) | is.element(as.character(pcol),mv)
  return(isna)
}
#' Mask individual identifiers
#'
#' Replace the identifiers of individuals in the pedigree with
#' auto generated names using alphabetic combinations.
#'
#' @param pedigree
#' A data frame with (at least) three columns that correspond to the
#' individual, male parent and female parent, respectively.
#'
#' @param data
#' The names column in data order (with dups); the default is the
#' first column of the pedigree.
#'
#' @param mv
#' missing value indicator.
#'
#' @return
#' A list with components \code{pedigree} and \code{data} containing
#' the pedigree with renamed individuals and their aliases, respectively.
#'
mask <- function(pedigree, data=pedigree[,1],
                 mv=c("0","*"," ")) {
  ## pedigree is a me,mum,dad,... data frame
  ## data is the names column in data order (with dups)

  if(!is.data.frame(pedigree))
    stop("pedigree must be a data frame")
  lped <- nrow(pedigree)

  ## Get levels ignoring missing value character(s)
  lvls <- unique(c(as.character(pedigree[,1]),
                   as.character(pedigree[,2]),
                   as.character(pedigree[,3])))
  lvls <- lvls[!is.missing(lvls,mv)]
  data <- as.character(data)
  if(any(is.na(ln <- match(data, lvls))))
    stop("Names in data but not in pedigree")

  ## convert to integer
  p123 <- matrix(c(match(pedigree[,1],lvls,nomatch=0),
                   match(pedigree[,2],lvls,nomatch=0),
                   match(pedigree[,3],lvls,nomatch=0)),nrow=lped,byrow=FALSE)
  alias <- aka(length(lvls))
  w <- p123 > 0
  p123[w] <- alias[p123[w]]
  pedigree[,1:3] <- data.frame(p123, stringsAsFactors=FALSE)
  return(list(pedigree=pedigree, data=alias[ln]))
}
