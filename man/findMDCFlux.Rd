\name{findMDCFlux}
\alias{findMDCFlux}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
find Manhaten-Distance-Closest Flux
}
\description{
 Given a wildtype flux (some fluxes can be absent i.e NA), find FBA solution of network given by model,
 such that Sum(|v_i-v_wt|) is minimal (like lMOMA but with option to exclude some reactions)
}
\usage{
findMDCFlux(model, wtflux, objVal = NA, pct_objective=100, 
        lpdir = SYBIL_SETTINGS("OPT_DIRECTION"), 
       solver = SYBIL_SETTINGS("SOLVER"), method = SYBIL_SETTINGS("METHOD"), 
       solverParm = data.frame(CPX_PARAM_EPRHS = 1e-06), verboseMode = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{model}{   An object of class \code{\link{modelorg}}.}
 
  \item{wtflux}{desired flux distribution, when wtflux is NA its constraint won't be included in LP. 
}
\item{pct_objective}{
  Biomass will be garnteed to be at least this value multiplied by max biomass calculated using standard FBA.
  Values are 0 to 100.
}
  \item{objVal}{
lower bound of objective function, if not set it will be set to maximum biomass.
}
 
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

\value{
return list with slot mdcflx containing the new calculated fluxes and status returned from solver.
}

\author{Abdelmoneim Amer Desouki}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{modelorg}},  
 \code{\link{lmoma}} 
}
\examples{
\dontrun{
## The function is currently defined as
function (model, wtflux, objVal = NA, lpdir = SYBIL_SETTINGS("OPT_DIRECTION"), 
    solver = SYBIL_SETTINGS("SOLVER"), method = SYBIL_SETTINGS("METHOD"), 
    solverParm = data.frame(CPX_PARAM_EPRHS = 1e-06), verboseMode = 2) 
{
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (is.na(objVal)) {
        sol = optimizeProb(model, solver = solver, method = method, 
            solverParm = solverParm)
        objVal = lp_obj(sol)
    }
    out <- FALSE
    nc <- react_num(model)
    nr <- met_num(model)
    nd <- sum(!is.na(wtflux))
    absMAX <- SYBIL_SETTINGS("MAXIMUM")
    nRows = nr + 2 * nd + 1
    nCols = nc + 2 * nd
    LHS <- Matrix::Matrix(0, nrow = nRows, ncol = nCols, sparse = TRUE)
    LHS[1:nr, 1:(nc)] <- S(model)
    ii = matrix(c((nr + 1):(nr + nd), which(!is.na(wtflux))), 
        ncol = 2)
    LHS[ii] <- 1
    if (nd == 1) {
        LHS[(nr + 1):(nr + nd), (nc + 1):(nc + nd)] <- 1
    }
    else {
        diag(LHS[(nr + 1):(nr + nd), (nc + 1):(nc + nd)]) <- 1
    }
    ii = matrix(c((nr + nd + 1):(nr + 2 * nd), which(!is.na(wtflux))), 
        ncol = 2)
    LHS[ii] <- -1
    if (nd == 1) {
        LHS[(nr + nd + 1):(nr + 2 * nd), (nc + nd + 1):(nc + 
            2 * nd)] <- 1
    }
    else {
        diag(LHS[(nr + nd + 1):(nr + 2 * nd), (nc + nd + 1):(nc + 
            2 * nd)]) <- 1
    }
    LHS[(nr + 2 * nd + 1), 1:nc] <- obj_coef(model)
    if (verboseMode > 2) 
        print(sprintf("nrows:\%d, ncols=\%d, nd=\%d,nc=\%d,nr=\%d", 
            nRows, nCols, nd, nc, nr))
    lower <- c(lowbnd(model), rep(0, 2 * nd))
    upper <- c(uppbnd(model), rep(absMAX, 2 * nd))
    rlower <- c(rep(0, nr), wtflux[!is.na(wtflux)], -wtflux[!is.na(wtflux)], 
        objVal)
    rupper <- c(rep(0, nr), rep(absMAX, 2 * nd + 1))
    cobj <- c(rep(0, nc), rep(1, 2 * nd))
    cNames = paste(c(rep("x", nc), rep("dp", nd), rep("dn", nd)), 
        c(1:nc, which(!is.na(wtflux)), which(!is.na(wtflux))), 
        sep = "_")
    switch(solver, glpkAPI = {
        out <- vector(mode = "list", length = 5)
        prob <- glpkAPI::initProbGLPK()
        rtype <- c(rep(glpkAPI::GLP_FX, nr), rep(glpkAPI::GLP_LO, 
            2 * nd))
        if (lpdir == "max") {
            rtype <- c(rtype, glpkAPI::GLP_LO)
        } else {
            rtype <- c(rtype, glpkAPI::GLP_UP)
        }
        TMPmat <- as(LHS, "TsparseMatrix")
        out[[1]] <- glpkAPI::addRowsGLPK(prob, nrows = nRows)
        outj <- glpkAPI::addColsGLPK(prob, ncols = nCols)
        mapply(setColNameGLPK, j = c(1:nCols), cname = cNames, 
            MoreArgs = list(lp = prob))
        glpkAPI::setObjDirGLPK(prob, glpkAPI::GLP_MIN)
        out[[2]] <- glpkAPI::loadMatrixGLPK(prob, length(TMPmat@x), 
            TMPmat@i + 1, TMPmat@j + 1, TMPmat@x)
        out[[3]] <- glpkAPI::setColsBndsObjCoefsGLPK(prob, c(1:nCols), 
            lower, upper, cobj)
        out[[4]] <- glpkAPI::setRowsBndsGLPK(prob, c(1:nRows), 
            rlower, rupper, rtype)
        parm <- sapply(dimnames(solverParm)[[2]], function(x) eval(parse(text = x)))
        val <- solverParm[1, ]
        if (method == "interior") {
            glpkAPI::setInteriorParmGLPK(parm, val)
            out[[5]] <- TRUE
        } else {
            glpkAPI::setSimplexParmGLPK(parm, val)
            out[[5]] <- TRUE
        }
        if (verboseMode > 2) {
            fname = format(Sys.time(), "glpk_ManhatDist_\%Y\%m\%d_\%H\%M.lp")
            print(sprintf("write problem: \%s/\%s", getwd(), fname))
            writeLPGLPK(prob, fname)
            print("Solving...")
        }
        lp_ok <- glpkAPI::solveSimplexGLPK(prob)
        lp_obj <- glpkAPI::getObjValGLPK(prob)
        lp_stat <- glpkAPI::getSolStatGLPK(prob)
        if (is.na(lp_stat)) {
            lp_stat <- lp_ok
        }
        lp_fluxes <- glpkAPI::getColsPrimGLPK(prob)
    }, cplexAPI = {
        out <- vector(mode = "list", length = 4)
        prob <- openProbCPLEX()
        out <- setIntParmCPLEX(prob$env, CPX_PARAM_SCRIND, CPX_OFF)
        chgProbNameCPLEX(prob$env, prob$lp, "ManhatenDist cplex")
        setObjDirCPLEX(prob$env, prob$lp, CPX_MIN)
        rtype <- c(rep("E", nr), rep("G", 2 * nd), "G")
        TMPmat <- as(LHS, "TsparseMatrix")
        out[[1]] <- newRowsCPLEX(prob$env, prob$lp, nRows, rlower, 
            rtype)
        out[[2]] <- newColsCPLEX(prob$env, prob$lp, nCols, cobj, 
            lower, upper)
        out[[3]] <- chgCoefListCPLEX(prob$env, prob$lp, length(TMPmat@x), 
            TMPmat@i, TMPmat@j, TMPmat@x)
        parm <- sapply(dimnames(solverParm)[[2]], function(x) eval(parse(text = x)))
        out[[4]] <- setDblParmCPLEX(prob$env, parm, solverParm)
        if (verboseMode > 2) {
            fname = format(Sys.time(), "Cplex_ManhatDist_\%Y\%m\%d_\%H\%M.lp")
            print(sprintf("write problem: \%s/\%s", getwd(), fname))
            writeProbCPLEX(prob$env, prob$lp, fname)
            print("Solving...")
        }
        lp_ok <- lpoptCPLEX(prob$env, prob$lp)
        lp_obj <- getObjValCPLEX(prob$env, prob$lp)
        lp_stat <- getStatCPLEX(prob$env, prob$lp)
        if (is.na(lp_stat)) {
            lp_stat <- lp_ok
        }
        lp_fluxes <- getProbVarCPLEX(prob$env, prob$lp, 0, nCols - 
            1)
    }, {
        wrong_type_msg(solver)
    })
    optsol <- list(ok = lp_ok, obj = lp_obj, stat = lp_stat, 
        fluxes = lp_fluxes, mdcflx = lp_fluxes[1:nc], wtflx = wtflux)
    return(optsol)
  }
}
}%dontrun
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manhaten distance closest flux }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
