\name{findFluxGeneExpr}
\alias{findFluxGeneExpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to find minimum set of required genes to get a given flux distribution
}
\description{
     given a flux distribution, use gpr to formulate a MILP to find 
     minimal set of required genes such that any GPR rule of a reaction 
     carrying a nonzero flux must be evaluated to TRUE. 
}
\usage{
findFluxGeneExpr(model, fluxes, threshold = 1e-06, 
      lpdir = SYBIL_SETTINGS("OPT_DIRECTION"), 
      solver = SYBIL_SETTINGS("SOLVER"), method = SYBIL_SETTINGS("METHOD"), 
      solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM"), verboseMode = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An object of class \code{\link{modelorg}}.
}
  \item{fluxes}{
  The flux distribution that we want to find the minimal set of genes able to produce it.
}
  \item{threshold}{
   Threshold in flux default value 1e-6. Fluxes with magnitude less than 
    threshold are considered zeros.
}
  \item{lpdir}{
  Character value, direction of optimisation.  Can be set to \code{"min"} or
    \code{"max"}.\cr
    Default: \code{SYBIL_SETTINGS("OPT_DIRECTION")}.
}
  \item{solver}{
Single character value.  The solver to use.  See
    \code{\link{SYBIL_SETTINGS}} for possible values.\cr
    Default: \code{SYBIL_SETTINGS("SOLVER")}.
}
  \item{method}{
 Single character value.  The optimization algorithm to use.  Possible values
    depend on the setting in \code{solver}.  See \code{\link{SYBIL_SETTINGS}} for
    possible values.\cr
    Default: \code{LP_METHOD(SYBIL_SETTINGS)}.
}
  \item{solverParm}{
    A data frame containing parameters for the specified solver.\cr
    Default: SOLVER_CTRL_PARM(SYBIL_SETTINGS).
}
\item{verboseMode}{
    An integer value indicating the amount of output to stdout:
    0: nothing, 1: status messages, 2: like 1 plus with more details,
    3: generates files of the LP problem.\cr
    Default: \code{2}.
}
}
\value{
return list of genes with State: "ON" if it is required, "OFF": when it is not required.
}
\author{Abdelmoneim Amer Desouki}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{modelorg}},  
\code{\link{optimizeProb}}, 
\code{\link{eFBA_gene}}
}

\examples{
\dontrun{

 data(iAF1260)
 model=iAF1260
 allgenes=allGenes(model)
 exoff=cbind(LOCUS=allgenes,State=rep(0,length(allgenes)));

 testgpr=(gpr(model)!="")
 table(testgpr)
 

 slvr="glpkAPI"
 seroff = eFBA_gene(model, exoff,solver=slvr,verbose=3,
          testgpr=testgpr,Tf=0.0001)
 mnflx=optimizeProb(model, algorithm = "mtf",solver=slvr);
 sfba=optimizeProb(model,solver=slvr)
 
 
 eFBAg_flx=seroff$rxn[, "newFlux"]
 geFBA=findFluxGeneExpr(model,fluxes(sfba)[fldind(sfba)],
         threshold=0.0001,verboseMode=4)
 gemnFBA=findFluxGeneExpr(model,fluxes(mnflx)[fldind(mnflx)],
         threshold=0.0001,verboseMode=4)
 
 geEFBAg=findFluxGeneExpr(model,as.numeric(seroff$rxn[, "newFlux"]),
        threshold=0.0001,verboseMode=4)
 
 # The minimum number of genes required to get each of the three fluxes
 print(cbind(EFBA_reqGenes=sum(geEFBAg$State=="ON"),
         minTotFlx_reqGene=sum(gemnFBA$State=="ON"),
       FBA_reqGenes=sum(geFBA$State=="ON")) )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ flux distribution }
\keyword{ gene expression }% __ONLY ONE__ keyword per line
