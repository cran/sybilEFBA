\name{gene2Rule}
\alias{gene2Rule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get reaction expression state from gene expression
}
\description{
 function to get state of expression of GPR rules from gene expression data

 # rules in brief:
 #1-Complexes: average, 2-isoenzymes: sum
 #3-multifunctioning: divide by count
 

 N.B.: GPR rules should be in form Sum-of-products (AND to OR) 

}
\usage{
gene2Rule(model, geneExpr, selected_rxns = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ An object of class \code{\link{modelorg}}.
}
  \item{geneExpr}{
  a data frame: geneID, expr_val for each gene.
}
  \item{selected_rxns}{
  optional parameter to select only a set of reactions.
}
}

\value{
 return list with main slot:

  ruleExpr: rxn_id,expr_val 
}

\author{Abdelmoneim Amer Desouki}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{modelorg}},  
\code{\link{FECorr}}
}
\examples{
\dontrun{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, geneExpr, selected_rxns = NULL) 
{
    if (length(selected_rxns) == 0) {
        ugpr = as.vector(unique(gpr(model)))
    }
    else {
        ugpr = as.vector(unique(gpr(model)[selected_rxns]))
    }
    ugpr = ugpr[ugpr != ""]
    gprExpr = NULL
    for (v_rule in ugpr) {
        rl = gsub("\\)", " ) ", v_rule)
        rl = gsub("\\(", " ( ", rl)
        pr = lapply(strsplit(unlist(strsplit(rl, " or ")), " and "), 
            function(x) gsub("[() ]", "", x))
        expr_val = 0
        for (p in 1:length(pr)) {
            gene_ind = match(pr[[p]], geneExpr$geneID)
            if (length(gene_ind) < length(pr[[p]])) {
                warning(sprintf("Rule \%s containing gene names not in geneID list, 
				term no: \%d term: \%s ", 
                  v_rule, p, pr[[p]][1]))
            }
            else {
                expr_val = expr_val + mean(geneExpr[gene_ind, 
                  "expr_val"])
            }
        }
        cnt = sum(gpr(model) == v_rule)
        gprExpr = rbind(gprExpr, cbind(rxn_id = react_id(model)[gpr(model) == 
            v_rule], expr_val = expr_val/cnt, gpr = v_rule, cnt = rep(cnt, 
            cnt)))
    }
    return(gprExpr)
  }
}%dont run
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
