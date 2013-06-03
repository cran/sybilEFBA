################################################
# Function: eFBA_rxn
#
# Performs a expression based flux balance analysis
# 
# Take into account the expression status of genes.
#  first step indicate FB as the max biomass using simpleFBA
# add FB*pct_objective>=Transpose(c) x as a new constraint in a MILP with a new ojective function
#      Mininmize: Sum(expr(g)<>flux(g)  // log2 expr level can be used as scaling
#   where expr(g)=1 if g expressed value>T1, else 0
#          flux(g)=0 if (all reactions catalyzed by g) have no flux(Threshold Tf and sum of all), 0 otherwise.
############
# irrev reactions: TWO identifying constraint
#		Eps*f <= x <= ub*f  : f binary flag =0 when x=0, 1 o.w.
#
# rev FOUR identifying constraints 2 flags(3 states of each flux and corresponding gene state) (0,d):0, (1,0):+ve, (1,1):-ve
#		-ub*f <= x <= ub*f
#		-My+Eps*f -  x  <= 0  ; y=1 when x<0 ,y=0 when x>0
# 4-  -M(1-y)+Eps*f + x <= 0

# TODO: add options to turn soft constraints to be obligatory
#		add parameter pct of objective to get as much considence as possible
 
eFBA_rxn <- function (model, expressionData,
                  Tf =1e-5,#SYBIL_SETTINGS("TOLERANCE"),  # threshold on flux 
				  pct_objective=100,
				   lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
				   solver = SYBIL_SETTINGS("SOLVER"),
				    method = SYBIL_SETTINGS("METHOD"),
				                                 #solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM")
				  solverParm=data.frame(CPX_PARAM_EPRHS=1e-6),
                  testgpr=NULL,
                  verboseMode = 2){
#PARAMETERS:
#===========
# model                 Sybil model structure (class modelorg)
# expressionData        class gxd (gene expression data): gene ID, and expression state
#						0: gene is off, 1 gene is expressed , NA: ignored/unknown			
# Tf                    Threshold value for flux(g) if(flux>Tf) it will be considered ON else it will be considered OFF
# testgpr				Boolean vector: containing flags of 

# RETURN
# optsol class
# Second obj function optimal value
##--------------------------------------------------------------------------##
## getRxnEqn: returns the equation for a given rxn, used for output
getRxnEqn=function(rxn=1){
	m=S(model)[,rxn]
	# input ==> output
	mcf=ifelse(abs(m)>1,paste("(",abs(m),")",sep=""),"") #!! <>1
	eqn=paste(gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m<0)],sep=""))),
		ifelse(react_rev(model)[rxn],"<==>","-->"),
	    gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m>0)],sep=""))) ,sep=" ")
	
return(eqn)
}
##--------------------------------------------------------------------------##
adjustTol <- function(val, tol = 0.00001) { #SYBIL_SETTINGS("TOLERANCE")

    if (!is(val, "numeric")) {
        stop( c("Argument val has to be numeric!") )
    }

    floorVal <- floor(val/tol)*tol

    return(floorVal)
}
##--------------------------------------------------------------------------##
 # check prerequisites 
    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }

##-------------------------------------------------------------------------------##

	# # Run FBA
		orig_sol = sybil::optimizeProb(model,solver=solver);
		FB =  lp_obj(orig_sol)*pct_objective/100# $obj;  ##objvalue sol.f

		origFlux=fluxes(orig_sol)[fldind(orig_sol)];#orig_sol$flux

	    if ( length(checkSolStat(lp_stat(orig_sol),solver))!=0 ){## checkSolStat
			print('No feasible solution - for FBA problem\n');
			break;
	    }
	  if (verboseMode > 3) {
	  	print(orig_sol);
	  	}
		
##-------------------------------------------------------------------------------##
# get flux status Flux(g)
# 
#	exprStatus=ifelse(expressionData$level<Te,0,1);
#   genflux as integer variables yi


#formulate new problem eFBA problem
    #1- add new n variables integer vars ,
    
    #2- add new constraint for FB
    #3- set new obj function if expr(gi)=0 then add yi to obj else add -yi. 
    #4- Add new n constraints(identifying yi's):  -Bi*yi <= vi <= Bi*yi

        
        
        #  the problem: minimize:
        #  
        #            |      |     
        #         S  |  0   |  = b
        #            |      |     
        #       ------------------
        #            |      |
        #        c^T |  0   | = FB
        #            |      |
        #       ------------------
        #            |      |
        #         -1 | -ub  | <= 0    ] 
        #            |      |           } -ub*yi <= vi <= ub * yi
        #         +1 | -ub  | <= 0    ]
        #            |      |
        #       ------------------
        #  lb   wt_lb|  0   |
        #  ub   wt_ub|  1   |
        #            |      |
        #  obj    0  |  2*ei-1    where ei =0 gene i NOT expressed/ 1 otherwise

        nc     <- react_num(model)
        nr     <- met_num(model)
        if(length( testgpr)==0){
                  gpr_rxn_ind=(gpr(model)!="" ) # only rxns having GPR rule
         }else { gpr_rxn_ind=testgpr;}
		#
			print(sprintf("%s: Calculate GPR state ....",format(Sys.time(), "%d-%m-%Y %X")))
            x <- logical(length(allGenes(model)));
		   	gi=expressionData[,1] %in% allGenes(model);##expressionData$Locus
		   	gj=match(expressionData[,1],allGenes(model) )  ;
		   	x <- rep(NA,length(x));  # missing genes will be ignored and not included in objective function or as 0 state.
			x[gj]=ifelse(expressionData[,2][gi]==0,FALSE,ifelse(expressionData[,2][gi]==1,TRUE,NA)); #NA undefined
			rxnStatus =rep(0,nc);# 0:will not be included in objective, -1: ON(Max), 1 OFF (Min)
		   	for (i in 1:nc){
		    		# test if gene is important or not in identifying rule state
		       		if(gprRules(model)[i]!="") {
						tmp=eval(parse(text = gprRules(model)[i]));
						rxnStatus[i]=ifelse(is.na(tmp),0,ifelse(tmp,-1,1));#coefficient of rxn in objective
						#=0 when no rule exists 
					}
			}
			print("state after GPR state evaluation:");
			print(table(rxnStatus))
		rxnStatus[!gpr_rxn_ind]=0;# exclude rxns not in testgpr
		print("state after testgpr:");
			print(table(rxnStatus))
		# test flux variabilty of rxns before adding constraints
		# if a reaction fixed or having a range different from gene state it will be ignored
			print(sprintf("%s: Calculate FVA ....",format(Sys.time(), "%d-%m-%Y %X")))
            rxnfva=(rxnStatus!=0)
			fv <- fluxVar(model,react=which(rxnfva),solver=solver)
			fv_min=lp_obj(fv)[1:sum(rxnfva)];fv_max=lp_obj(fv)[(sum(rxnfva)+1):(2*sum(rxnfva))]
		
		#print(length(fv
		rxnStatus[react_pos(react(fv))][abs(fv_min-fv_max)<Tf]=0;# exclude rxns with fixed values
		#state fixed: exclude rxns that must be ON (state can't be changed
		rxnStatus[react_pos(react(fv))][(fv_min>Tf) | (fv_max < -Tf)]=0;# cutoff
		print(sprintf("Number of rxns in FVA: %d, Number of rxns found fixed: %d,
		\n rxns that are always ON: %d ",sum(rxnfva),sum(abs(fv_min-fv_max)<Tf),sum(((fv_min>Tf) | (fv_max < -Tf)) & abs(fv_min-fv_max)>=Tf) ))
    	print("state after FVA:");
		print(table(rxnStatus))
		gpr_rxn_ind[rxnStatus==0]=F;
		rxnstname=ifelse(rxnStatus==0,"No rule/Unk",ifelse(rxnStatus==-1,"ON","OFF"));
		print(table(rxnstname))
		
		nGpr=sum(gpr_rxn_ind);
        is_irrev=(gpr_rxn_ind & lowbnd(model)>=0);
		is_irrev_on=(is_irrev & rxnStatus==-1);
		is_irrev_off=(is_irrev & rxnStatus==1);
        nIrrev=sum(is_irrev);
		nIrrev_on=sum(is_irrev_on);
		nIrrev_off=sum(is_irrev_off);
		
		is_rev=(gpr_rxn_ind & lowbnd(model)<0);
        nRev=sum(is_rev);
        is_rev_on=(is_rev & rxnStatus==-1);
		is_rev_off=(is_rev & rxnStatus==1);
        nRev_on=sum(is_rev_on);
		nRev_off=sum(is_rev_off);
		
		print(table(rxnStatus))
		print("gpr_rxn_ind final:")
		print(table(gpr_rxn_ind))
		# N.B: Eps minimal flux should be less than all fluxes and more than tolerance in constraints !!!
		print(sprintf("default tolerance:%f",Tf));
        Eps=Tf/3;#min(abs(origFlux[origFlux!=0]))*1e-1; # minimal flux, when it is equal to Tf, it was ignored
		# if(solverParm>Eps)solverParm[1,1]=Eps*1e-2; # important when usin CPLEX solver
		#if(solver="cplexAPI"){			Eps=1e-15;			#Tf=1		}
        INF=max(uppbnd(model));
        # ---------------------------------------------
        # constraint matrix
        # ---------------------------------------------

       # the initial matrix dimensions
       # LHS <- as.matrix.csr(0, nrow = nr+2*nIrrev+4*nRev+1, ncol = nc+nIrrev+2*nRev   )
       LHS <- Matrix::Matrix(0, nrow = nr+nIrrev+2*nRev+1, ncol = nc+nIrrev+nRev_off+2*nRev_on   )

       # rows for the initial S
       LHS[1:nr,1:nc] <- S(model)

       # fix the value of the objective function
       LHS[(nr+1),1:nc] <- obj_coef(model)

       # rows for the identifying constraint of flux variables
       #  Eps*f<=x<=ub*f
       # irrev( Eps*f<= x <=ub*f ) : 1-  Eps*f -  x  <= 0 : if state is ON
	#if(nIrrev>0){  # update 22/4/2013: use one variable and one constraint in case of irreversible rxn
		 if(nIrrev_on>0){
				 ii=matrix(c((nr+2)   :(nr+nIrrev_on+1)  ,which(is_irrev_on ) ),ncol=2)
			 
				 LHS[ii ] <- -1;
				 if(nIrrev_on==1){
					LHS[((nr+2)   :(nr+nIrrev_on+1)) , (nc+1) ]= Eps
				 }else{ 
					diag(LHS[((nr+2)   :(nr+nIrrev_on+1)) , ((nc+1):(nc+nIrrev_on)) ] ) = Eps;
				 }
		}
		   #           2-  -ub*f + x <= Eps : if state is OFF
		 if(nIrrev_off>0){
			 ii=matrix(c((nr+nIrrev_on+2)   :(nr+nIrrev+1)  , which(is_irrev_off) ),ncol=2)
				 LHS[ii ] <- 1;
			 if(nIrrev_off==1){
				LHS[((nr+nIrrev_on+2)   :(nr+nIrrev+1))  ,    ((nc+nIrrev))]   = -INF;
			 }else{
					diag ( LHS[((nr+nIrrev_on+2)   :(nr+nIrrev+1))  ,    ((nc+nIrrev_on+1):(nc+nIrrev))]  ) = -INF;
				}
		 }
    #}
  #  -ub*f<=x<=ub*f
       # rev : 1-  -ub*f -  x  <= Eps
	   print(nRev);
  #if(nRev>0){
	if(nRev_off>0){# Rev/OFF state
		#- c1offp: + x_1 - 1000 f_1 <= Eps
		ii=matrix(c((nr+nIrrev+2)   :(nr+nIrrev+nRev_off+1)  ,which(is_rev_off) ),ncol=2)
		LHS[ii ] <- -1;
		diag(LHS[((nr+nIrrev+2)   :(nr+nIrrev+nRev_off+1)) , ((nc+nIrrev+1):(nc+nIrrev+nRev_off)) ] ) = -INF;
       # c1offn: - x_1 - 1000 f_1 <= 0.000002
		print("  #            2-  c1offn: - x_1 - 1000 f_1 <= Eps ")
		ii=matrix(c((nr+nIrrev+nRev_off+2)   :(nr+nIrrev+2*nRev_off+1)  , which(is_rev_off) ),ncol=2)
		 LHS[ii ] <- -1;
		diag ( LHS[((nr+nIrrev+nRev_off+2)   :(nr+nIrrev+2*nRev_off+1))  ,    ((nc+nIrrev+1):(nc+nIrrev+nRev_off))]  ) = -INF;
	}
    
	print(c("ON:",nRev_on));
	if(nRev_on>0){		# state ON
		#print("# rev : 3-  -My+Eps*f -  x  <= 0  ; y=1 when x<0 ,y=0 when x>0");
		#Rev: ON bj coef -1 (Max), y mandatory to disable one constraint according to sign(x)
		print("#c1onp:  Rev/ON 1-:  - x_1 + Eps f_1 - 1000 y_1 <= 0");

		ii=matrix(c((nr+nIrrev+2*nRev_off+2)   :(nr+nIrrev+2*nRev_off+nRev_on+1)  ,which(is_rev_on) ),ncol=2)
		LHS[ii ] <- -1;
		diag(LHS[((nr+nIrrev+2*nRev_off+2)   :(nr+nIrrev+2*nRev_off+nRev_on+1)) , ((nc+nIrrev+nRev_off+1):(nc+nIrrev+nRev)) ] ) = Eps;
		diag(LHS[((nr+nIrrev+2*nRev_off+2)   :(nr+nIrrev+2*nRev_off+nRev_on+1)) , ((nc+nIrrev+nRev+1):(nc+nIrrev+nRev+nRev_on)) ] ) = -INF;
       
		print("#c1onn:    x_1 + 0.000002 f_1 + 1000 y_1 <= 1000");
		#print("#             4-  -M(1-y)+Eps*f + x <= 0 ")
		ii=matrix(c((nr+nIrrev+2*nRev_off+nRev_on+2)   :(nr+nIrrev+2*nRev+1)  ,which(is_rev_on) ),ncol=2)
		LHS[ii ] <- 1;  # x
		diag(LHS[((nr+nIrrev+2*nRev_off+nRev_on+2)   :(nr+nIrrev+2*nRev+1)) , ((nc+nIrrev+nRev_off+1):(nc+nIrrev+nRev)) ] ) = Eps; #f
		diag(LHS[((nr+nIrrev+2*nRev_off+nRev_on+2)   :(nr+nIrrev+2*nRev+1)) , ((nc+nIrrev+nRev+1):(nc+nIrrev+nRev+nRev_on)) ] ) = INF;#y
    }
	   # ---------------------------------------------
       # lower and upper bounds
       # ---------------------------------------------
	
	lower  <- c(lowbnd(model), rep(0, nIrrev+nRev+nRev_on ))
	upper  <- c(uppbnd(model), rep(1, nIrrev+nRev+nRev_on ))
	#RHS 1:nr+1 : as FBA model , FB,
	rlower <- c(rep(0,nr), FB,rep(-INF,nIrrev),rep(-INF,nRev),rep(-2*INF,nRev_on))	# will not be used in identifying constraints
	rupper <- c(rep(0,nr), FB,rep(0 ,nIrrev),rep(Eps ,2*nRev_off),rep(0 ,nRev_on),rep(INF ,nRev_on))
#rupper <- c(rep(0,nr), FB,rep(0 ,nIrrev),rep(Eps ,nIrrev),rep(Eps ,2*nRev),rep(0 ,nRev),rep(INF ,nRev))


       # ---------------------------------------------
       # objective function
       # ---------------------------------------------
       # Notes: 1-use reaction status instead of gene status
       # Minimize(ifelse(expr(rxn(i),-yi,yi)     math. |x-y|=x-y  when x>=y and y-x when y>=x 
       # then the difference will be Sum(rxn(i))+obj 
       
	    cobj <- c(rep(0, nc), rxnStatus[gpr_rxn_ind],rep(0, nRev_on))

        # ---------------------------------------------
        cNames=paste(c(rep("x",nc),rep("Irv",nIrrev),rep("Rev",nRev),rep("y",nRev_on)),
		 		 		               c(1:nc,which(is_irrev),which(is_rev),which(is_rev_on)),sep="_" ) ;
        # ---------------------------------------------
        # build problem object
        # ---------------------------------------------
        # can it be done independent from solver!!
		
        switch(solver,
            # ----------------------- #
            "glpkAPI" = {
				print(sprintf("%s : creating MILP using GLPK solver....",format(Sys.time(), "%d-%m-%Y %X")))
                out <- vector(mode = "list", length = 4)
                prob <- initProbGLPK();# new problem
               
                out[[1]] <-  addRowsGLPK(prob, nrows=nr+nIrrev+2*nRev +1)
				outj <- addColsGLPK(prob, ncols=nc+nIrrev+nRev+nRev_on )
				#setColNameGLPK(prob,c(1:(nc+nIrrev+2*nRev )),cNames );
				mapply(setColNameGLPK, j = c(1:(nc+nIrrev+nRev+nRev_on )), cname = cNames, MoreArgs = list(lp = prob));
				setObjDirGLPK(prob, GLP_MIN);
                ## note: when FX or LO value taken from lb and when UP take from UP
               rtype <- c(rep(GLP_FX, nr),GLP_LO, rep(GLP_UP, nIrrev+2*nRev  ))# objective >=
	       
	        # set the right hand side Sv = b

				out[[4]] <- setRowsBndsGLPK(prob, c(1:(nr+nIrrev+2*nRev +1)), lb=rlower, ub=rupper,type=rtype )
        
	        # add upper and lower bounds: ai <= vi <= bi
                cc <- setColsBndsObjCoefsGLPK(prob, c(1:(nc+nIrrev+nRev+nRev_on )),   lower,   upper,  cobj   )
                if (verboseMode > 2) { print(cc);	}		
				#nzLHS=nzijr(LHS);
        	#print(nzLHS);
				TMPmat <- as(LHS, "TsparseMatrix")
                #cc <- loadMatrixGLPK(prob,   nzLHS$ne,   nzLHS$ia,  nzLHS$ja,  nzLHS$ar     )
				cc <- loadMatrixGLPK(prob,length(TMPmat@x),ia  = TMPmat@i + 1,ja  = TMPmat@j + 1,ra  = TMPmat@x)
                if (verboseMode > 3) {  print(c("Load matrix GLPK return code:",cc));}# technical msgs
                ctype <- c(rep(GLP_CV, nc), rep(GLP_BV,nIrrev+nRev+nRev_on ));
				setColsKindGLPK(prob,c(1:(nc+nIrrev+nRev+nRev_on )),ctype);

               if (verboseMode > 2) {                      
					fname=format(Sys.time(), "glpk_eFBA_%Y%m%d_%H%M.lp");
					print(sprintf("Writing problem to file: %s/%s ...",getwd(),fname));
                	writeLPGLPK(prob,fname);
                	print(format(Sys.time(), "Testing time : %Y%m%d %X Solving..."));
                }
				
				## Solve
            setMIPParmGLPK(PRESOLVE,GLP_ON);
			lp_ok=solveMIPGLPK(prob);
			return_codeGLPK(lp_ok);
			lp_stat=mipStatusGLPK(prob);
			status_codeGLPK(lp_stat);
			lp_obj=mipObjValGLPK(prob);
			colst=mipColsValGLPK(prob);
			newFlux=colst
	       ## --- 
	        #newFlux=adjustTol(newFlux,tol=Tf);
	        newFlux=newFlux[1:nc];
	        newStat=ifelse(abs(newFlux)>Tf,1,0);
            },
            # -------------------------------------------------------------------------------------------- #
           "cplexAPI" = {
                print(sprintf("%s : creating MILP using CPLEX solver....",format(Sys.time(), "%d-%m-%Y %X")))
				out <- vector(mode = "list", length = 3)
				prob <- openProbCPLEX()
				out <- setIntParmCPLEX(prob$env, CPX_PARAM_SCRIND, CPX_OFF)
		                
                chgProbNameCPLEX(prob$env, prob$lp, "eFBA rxn cplex");
                #E,L,"G"
                rtype <- c(rep("E",nr),"G", rep("L", nIrrev+2*nRev  ))
                setObjDirCPLEX(prob$env, prob$lp, CPX_MIN);

                out[[1]] <- newRowsCPLEX(prob$env, prob$lp,
                                         nrows=nr+nIrrev+2*nRev +1, rhs=rupper, sense=rtype)

                out[[2]] <- newColsCPLEX(prob$env, prob$lp,
                                         nc+nIrrev+nRev+nRev_on , obj=cobj, lb=lower, ub=upper,cnames=cNames)
				
				print(sprintf("%s : step 2: nzijr....",format(Sys.time(), "%d-%m-%Y %X")))
                        # constraint matrix
				TMPmat <- as(LHS, "TsparseMatrix")
				out[[3]] <- chgCoefListCPLEX(prob$env, prob$lp,#lp@oobj@env, lp@oobj@lp,
                                   nnz = length(TMPmat@x),
                                   ia  = TMPmat@i,
                                   ja  = TMPmat@j,
                                   ra  = TMPmat@x);
							
#                nzLHS=nzijr();#LHS    
#print(sprintf("%s : step 3: chngCoef....",format(Sys.time(), "%d-%m-%Y %X")))				

				ctype<-c(rep('C', nc), rep('B',nIrrev+nRev+nRev_on ));
				status = copyColTypeCPLEX (prob$env, prob$lp, ctype);
				check <- setObjDirCPLEX(prob$env, prob$lp, CPX_MIN);
		
                if (verboseMode > 2) {                      
                      fname=format(Sys.time(), "Cplex_eFBA_%Y%m%d_%H%M.lp");
					  print(sprintf("Writing problem to file: %s/%s  ...",getwd(),fname));
					  writeProbCPLEX(prob$env, prob$lp, fname);
		       }
			   
			   #set precision 
                parm <- sapply(dimnames(solverParm)[[2]],
                                     function(x) eval(parse(text = x)))
			    if(solverParm>Eps)   
					solverParm[1,1]=Tf/10;
				print(c("solverparam",solverParm));
			   out[[4]]  <- cplexAPI::setDblParmCPLEX(prob$env, parm, solverParm);
				print(unlist(out))
			print(sprintf("%s : Solving MILP using CPLEX solver....",format(Sys.time(), "%d-%m-%Y %X")))
#		       print("Solving...");		 
			   lp_ok     <- mipoptCPLEX(prob$env, prob$lp);
	           print(lp_ok);
				sol=solutionCPLEX(prob$env, prob$lp);
	       
				if (verboseMode > 3) {print(sol);}	
               
               lp_obj=sol$objval;
               lp_stat   <- getStatCPLEX(prob$env, prob$lp)
              
                colst=sol$x;
                
               newFlux=sol$x;#adjustTol(sol$x,tol=Tf);
               newFlux=newFlux[1:nc];#May cause problems if Auxiliary variable were at the begining
               newStat=ifelse(abs(newFlux)>Tf,1,0);
               # when using binary variables as indicators: .floorValues(.ceilValues(newFlux[(nc+1):(2*nc)],tol=Tf/10),tol=Tf);
           },
           # ----------------------- #
           {
               wrong_type_msg(solver)
            }
        )
			print(sprintf("%s : Preparing return values ....",format(Sys.time(), "%d-%m-%Y %X")))
        
 	       print(sprintf("Number of Rxns with gene ON=%d ",length(rxnStatus[rxnStatus==-1]) ));

	       print(sprintf("Rxns with gene OFF=%d",length(rxnStatus[rxnStatus==1]) ));  
	       print(sprintf("Rxns with no rules=%d ",length(rxnStatus[rxnStatus==0]) ));
	       print(sprintf("Total number of rxns: %d",length(rxnStatus) ));
	       print(sprintf("The differnece is: %.0f ", (lp_obj+length(rxnStatus[rxnStatus==-1])) ));

  		#excReact = findExchReact(model)[1];# 1 is position
		#excReactPos=react_pos(excReact$exchange);
		excReact = findExchReact(model);
        excReactPos=(react_id(model) %in% react_id(excReact));#excReact$exchange

                is_excReact=((c(1:length(react_id(model)))) %in% excReactPos);
                
                origStat=ifelse(abs(origFlux)>Tf,1,0);
               
               rxn_st<-cbind(gpr=gpr(model),react=react_id(model),reactName=react_name(model),is_excReact=is_excReact,
                                 expr=rxnstname,lb=lowbnd(model),
                              ##           eqns=sapply(c(1:length(react_id(model))),function(x)getRxnEqn(x)) ,
                                          origFlux,origStat, newFlux,newStat,difference=abs(origStat-newStat));
   flxst=rep(NA,nc);
   flxst[which(is_irrev)]=colst[(nc+1):(nc+nIrrev)];
   flxst[which(is_rev)]=colst[(nc+nIrrev+1):(nc+nIrrev+nRev)];
   
  remove(prob);
# ----------------------------------------  ***  -------------------------------- #
	optsol <- list(        ok = lp_ok,
			       obj = FB,
			       stat = lp_stat,
			       fluxes = newFlux,
			       origfluxes=origFlux,
			       rxn=rxn_st,
			       diff_obj=lp_obj,
			       colst=colst,
			       flxst=flxst,
			       cNames=cNames
			  )
	return(optsol);
}
#update one constraint
# if(nIrrev>0){
		 # ii=matrix(c((nr+2)   :(nr+nIrrev+1)  ,which(is_irrev) ),ncol=2)
		 # LHS[ii ] <- -1;
		 # if(nIrrev==1){
			# LHS[((nr+2)   :(nr+nIrrev+1)) , ((nc+1):(nc+nIrrev)) ]
		 # }else{ 
			# diag(LHS[((nr+2)   :(nr+nIrrev+1)) , ((nc+1):(nc+nIrrev)) ] ) = Eps;
		 # }
		  ##           2-  -ub*f + x <= Eps : if state is OFF
		 # ii=matrix(c((nr+nIrrev+2)   :(nr+2*nIrrev+1)  , which(is_irrev) ),ncol=2)
			 # LHS[ii ] <- 1;
		 # diag ( LHS[((nr+nIrrev+2)   :(nr+2*nIrrev+1))  ,    ((nc+1):(nc+nIrrev))]  ) = -INF;
    # }
