\name{FECorr}

\alias{FECorr}

\encoding{utf8}

\title{ Function: FECorr: Flux Expression Correlation}
\description{
 This function uses FVA under different conditions to find fluxes that linearly correlates to corresponding gene expression. 
}

\usage{
FECorr (model, nCond, initCond, geneExpressionData=NULL,RuleExpressionData=NULL,
				  pct_objective=100,
				  selected_rxns=NULL,
				  only_identified_rules=FALSE,
				  minExprFoldChange=0,
				 lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
				  solver = SYBIL_SETTINGS("SOLVER"),
			      method = SYBIL_SETTINGS("METHOD"),
				  solverParm=data.frame(CPX_PARAM_EPRHS=1e-6),
               verboseMode = 2)
}

\arguments{
  \item{model}{   An object of class \code{\link{modelorg}}.}
 \item{nCond}{	Number of conditions (FBA problems to be solved)}
 \item{geneExpressionData}{ a data frame: geneID,Cond_id, ExpressionVal
 		column rows are genes and column j+1 is representing gene expression under condition j
 		}
 \item{RuleExpressionData}{rxn_id,cond_id, ExpressionVal}		
 \item{initCond}{rxn_id,cond_id,lb,ub,objcoef : lower and upper bounds of rxns under different conditions that represents 
        the available nutrients under these conditions
}
 \item{pct_objective}{
  Biomass will be garnteed to be at least this value multiplied by max biomass calculated using standard FBA.
  Values are 0 to 100.
}
 \item{selected_rxns}{optional parameter used to select a set of reactions not all, Boolean with the same length react_id(model)}
 \item{only_identified_rules}{ ignore rxns containing genes with unidentified expression}
 \item{minExprFoldChange}{can be used to consider only genes with a significant change in expression level
     (i.e min(expression(gene))*minExprFoldChange*2 must be less than or equal to max(expression(gene)))}

 \item{lpdir}{  Character value, direction of optimisation.  Can be set to \code{"min"} or
    \code{"max"}.\cr
    Default: \code{SYBIL_SETTINGS("OPT_DIRECTION")}.
}
 \item{solver}{
       Single character string giving the solver package to use.  See
       \code{\link{SYBIL_SETTINGS}} for possible values.\cr
       Default: \code{SYBIL_SETTINGS("SOLVER")}.
     }
     \item{method}{
       Single character string giving the method the desired solver has to use.
       \code{\link{SYBIL_SETTINGS}} for possible values.\cr
       Default: \code{SYBIL_SETTINGS("METHOD")}.
     }
     \item{solverParm}{
       A named data frame or list containing parameters for the specified
       solver.  Parameters can be set as data frame or list:
       \code{solverParm = list(parm1 = val1, parm2 = val2)} with \code{parm1} and
       \code{parm2} being the names of two different parameters and \code{val1}
       and \code{val2} the corresponding values. For possible parameters and
       values see the documentation of the used solver package
       (e.g. \pkg{glpkAPI}).\cr
       Default: \code{SYBIL_SETTINGS("SOLVER_CTRL_PARM")}.
    }
 
  \item{verboseMode}{
    An integer value indicating the amount of output to stdout:
    0: nothing, 1: status messages, 2: like 1 plus with more details,
    3: generates files of the LP problem.\cr
    Default: \code{2}.
}
}
\details{
Main steps
    1- Run FVA for all conditions, exclude rxns fixed in all conditions
    2- Identify ruleExpression for set of rxns remaining from 1,
    3- Fit Expr to FVA range.
    4- Run findMDCFlux to find closest genome-scale flux
    5- Recalculate correlation: Posterior, iFlux, ruleExpr

}
\value{
returns a list containing slots:
	 geneID,slope,intercept,base_level: OGOR: one gene one rxn
	 iFlux: new fluxes calculated at all the given conditions which considers gene expression data to get a linear fit between 
	 gene expression and fluxes. It is a data frame containing the following columns 
	 	rxn_id: reactionId in model,cond_id:,lb,ub,objCoef,xpc_flux,fva_min,fva_max,RuleExprVal,iflux

}

\author{Abdelmoneim Amer Desouki}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{modelorg}},  
\code{\link{optimizeProb}}, 
\code{\link{gene2Rule}}
}
\examples{
\dontrun{
	library(sybil)
	data(iAF1260)
	model= iAF1260
	# trivial test 1, linear levels for an unbounded rxn
	ncnd=3
	

	rxn="R_ADK1"
	nc=length(react_id(model))
	initCond=cbind(rxn_id=react_id(model),cond_id=rep(1,nc),
	         lb=lowbnd(model),ub=uppbnd(model),obj=obj_coef(model))
	initCond=rbind(initCond,cbind(rxn_id=react_id(model),
	        cond_id=rep(2,nc),lb=lowbnd(model),ub=uppbnd(model),
	        obj=obj_coef(model)))
	initCond=rbind(initCond,cbind(rxn_id=react_id(model),
	     cond_id=rep(3,nc),lb=lowbnd(model),ub=uppbnd(model),
		 obj=obj_coef(model)))
	
	cnds=(1:3)
	gprExp=cbind(rxn_id=rxn,cond_id=1,expr_val=2)
	gprExp=rbind(gprExp,cbind(rxn_id=rxn,cond_id=2,expr_val=4))
	gprExp=rbind(gprExp,cbind(rxn_id=rxn,cond_id=3,expr_val=6))
	
	fcflx=FECorr(model,nCond=ncnd,initCond=initCond,
	       RuleExpressionData=gprExp,selected_rxns=(react_id(model)==rxn),
		   verboseMode=4);
	fcflx[[2]][fcflx[[2]][,2]==rxn,]
	
	
	cor(as.numeric(fcflx[[2]][fcflx[[2]][,2]==rxn,"expr_val"]),
	               as.numeric(fcflx[[2]][fcflx[[2]][,2]==rxn,"iflx"]))


}% end dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FBA }
\keyword{ gene expression }
\keyword{ EFBA reaction }% __ONLY ONE__ keyword per line
