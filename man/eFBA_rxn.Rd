\name{eFBA_rxn}

\alias{eFBA_rxn}

\encoding{utf8}

\title{ Function: eFBA}
\description{
This function performs an expression based flux balance analysis.
It takes into account the expression status of genes. First step calculate FB as the max biomass using standard FBA
use FB=cTx as a new constraint in a MILP with a new ojective function
     Mininmize: Sum(expr(gpr(r)) != flux(r)  
  where expr(r)=1 if gpr rule is ON (use boolean logic to evaluate rules depending on given gene state)
         flux(r)=0 if reaction r has no flux (Threshold Tf), 0 otherwise.

}

\usage{
eFBA_rxn(model, expressionData,
                  Tf =0.01,  
            pct_objective=100,
  lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
  solver = SYBIL_SETTINGS("SOLVER"),
      method = SYBIL_SETTINGS("METHOD"),
  solverParm=data.frame(CPX_PARAM_EPRHS=1e-6),testgpr=NULL,
              verboseMode = 2)
}

\arguments{
  \item{model}{   An object of class \code{\link{modelorg}}.}
  \item{expressionData}{
 a dataframe (gene expression data): gene ID, and expression status (0 :OFF or 1: ON)
}
  \item{pct_objective}{
  Biomass will be garnteed to be at least this value multiplied by max biomass calculated using standard FBA.
  Values are 0 to 100.
}
  \item{Tf}{
   Threshold in flux default value 1e-5. Fluxes with magnitude less than Tf are considered zeros.
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

\item{testgpr}{Boolean flag containing flags of whether a reaction is considered or not} 
  \item{verboseMode}{
    An integer value indicating the amount of output to stdout:
    0: nothing, 1: status messages, 2: like 1 plus with more details,
    3: generates files of the LP problem.\cr
    Default: \code{2}.
}

}
\details{
formulates the following LP:
    the problem: minimize:
         
                   |      |     
                S  |  0   |  = b
                   |      |     
              ------------------
                   |      |
               c^T |  0   | = FB
                   |      |
              ------------------
                   |      |
                -1 | -ub  | <= 0    ] 
                   |      |           > -ub*yi <= vi <= ub * yi
                +1 | -ub  | <= 0    ]
                   |      |
              ------------------
         lb   wt_lb|  0   |
         ub   wt_ub|  1   |
                   |      |
         obj    0  |  2*ei-1    where ei =0 gene i NOT expressed/ 1 otherwise


}
\value{
returns a list containing slots:
	origfluxes: original fluxes calculted by standard FBA.
	fluxes: new fluxes that better matches gene expression. 
	gene_stat: state of genes (0 or 1) according to calculated fluxes
	stat: solution status returned from solver
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
	## The examples here require the package glpkAPI to be
	## installed. If that package is not available, you have to set
	## the argument 'solver' (the default is: solver = "glpkAPI").

	## load the example data set
	data(iAF1260)
        model=iAF1260
	allgenes=allGenes(model)
	##set all genes to OFF
	ex=cbind(LOCUS=allgenes,State=rep(0,length(allgenes)));
	##set state of genes as required
	
	## Optionally optimize only for a set of reactions
	testgpr=rep(FALSE,react_num(model))
	testgpr[react_id(model) \%in\% c('R_PYK','R_FLDR','R_FEROpp')]=TRUE
	 
	 # read or set the expression status of genes : 0 not expressed, 
	 #     1 Expressed, NA: unknown or ignored
	 ex_on=ex
	 ex_on[ex_on[,1] \%in\% c('b0123','b1854','b0684','b3924'),2]=1
	 
	 ex_off=ex# all genes are not expressed
	
	 slvr="glpkAPI"
	# run two times, one with genes ON, the other with genes OFF
	 ser_on = eFBA_rxn(model, ex_on,solver=slvr,verbose=3,testgpr=testgpr,Tf=0.001);
	 ser_off = eFBA_rxn(model, ex_off,solver=slvr,verbose=3,testgpr=testgpr,Tf=0.001)
	
	 #print results
	 cbind(rxn=c('R_PYK','R_FLDR','R_FEROpp'),geneONfluxes=ser_on$rxn[testgpr, "newFlux"],
	            geneOFFfluxes=ser_off$rxn[testgpr, "newFlux"],
	            ruleStateON=ser_on$rxn[testgpr,"expr"],ruleStateOFF=ser_off$rxn[testgpr,"expr"],
	            FBA_flux=ser_off$rxn[testgpr,"origFlux"])
	
	#Count difference between these two solutions
	sum(abs( as.numeric(ser_on$rxn[, "newFlux"])-as.numeric(ser_off$rxn[, "newFlux"]))>0.001)


}% end dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FBA }
\keyword{ EFBA reaction }% __ONLY ONE__ keyword per line
