################################################
# Function: eFBA
#
# Performs an expression based flux balance analysis
# 
# Take into account the expression status of genes
#  first step ind FB as the max biomass using simpleFBA
# use FB=cTx as a new constraint in a MILP with a new ojective function
#      Mininmize: Sum(expr(g)<>flux(g)  // log2 expr level can be used as scaling
#   where expr(g)=1 if g expressed value>T1, else 0
#          flux(g)=0 if (all reactions catalyzed by g) have no flux(Threshold T2 and sum of all), 0 else.
############

eFBA_gene <- function (model, expressionData,
                  Tf =1e-5,#SYBIL_SETTINGS("TOLERANCE"),  # threshold on flux
				  pct_objective=100,
 		           lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
			   solver = SYBIL_SETTINGS("SOLVER"),
			    method = SYBIL_SETTINGS("METHOD"),
					 #solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM")
			  solverParm=data.frame(CPX_PARAM_EPRHS=1e-6),
			  testgpr=NULL,
                  verboseMode = 2, ...){
#PARAMETERS:
#===========
# model                 Sybil model structure (class modelorg)
# expressionData        a dataframe (gene expression data): gene ID, and expression level
# Tf                    Threshold value for flux(g)

# RETURN
# optsol class
# Second obj function optimal value

## getRxnEqn: returns the equation for a given rxn, used for output
getRxnEqn=function(rxn=1){
	m=S(model)[,rxn]
	# input ==> output
	mcf=ifelse(abs(m)>1,paste("(",abs(m),")",sep=""),"")
	eqn=paste(gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m<0)],sep=""))),
		ifelse(react_rev(model)[rxn],"<==>","-->"),
	    gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m>0)],sep=""))) ,sep=" ")
	
return(eqn)
}
##--------------------------------------------------------------------------##
 # check prerequisites 
    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }
    
#require(sybilCPLEX)
##-------------------------------------Prepare Problem object --------------------##
# get OptObj instance
#lp <- sybil::sysBiolAlg(model, algorithm = "fba", ...)

			 
##-------------------------------------------------------------------------------##


	# # Run FBA
          orig_sol = sybil::optimizeProb(model,solver=solver);# solver and method already in lp
		FB =  lp_obj(orig_sol)*pct_objective/100# $obj;  ##objvalue sol.f
	    if ( length(checkSolStat(lp_stat(orig_sol),solver))!=0 ){## checkSolStat
			print('No feasible solution - for initial FBA problem\n');
			break;
	    }
	  if (verboseMode > 3) {
	  	print(orig_sol);
	  	}
	   #if (length(max_growth_rate)>0 ) return (orig_sol);
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

        nc     <- react_num(model)
        nr     <- met_num(model)
        ng    <- length(allGenes(model))
        geneNames=allGenes(model)
        if(length( testgpr)==0){
                  gpr_rxn_ind=(gpr(model)!="" ) # only rxns having GPR rule
         }else { gpr_rxn_ind=testgpr;}
		
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
		
	    gpr_rxn=sum(gpr_rxn_ind);
		rules=gpr(model)[gpr_rxn_ind]
         
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


        # ---------------------------------------------
        # constraint matrix
        # ---------------------------------------------

       # the initial matrix dimensions  variables:   nc:rates, nc:flux(r)(0,1) , ng
       #sum((gpr(model)!="") & !react_rev(model)) + sum((gpr(model)!="") & react_rev(model)) * 2;
      
       #LHS <- as.matrix.csr(0, nrow = nr+gpr_rxn*3+1, ncol = (nc+gpr_rxn+ng))
       LHS <- Matrix::Matrix(0, nrow = nr+gpr_rxn*3+1, ncol = (nc+gpr_rxn+ng), sparse = TRUE)

       # rows for the initial S
       LHS[1:nr,1:nc] <- S(model)

       # fix the value of the objective function
       LHS[(nr+1),1:nc] <- obj_coef(model)

       # rows for the identifying constraint of flux variables
       #diag is wrong : it cant be seq, only the first 810 will be used regardless of gpr existence
	if(gpr_rxn>0){

       ii=matrix(c((nr+2)   :(nr+gpr_rxn+1)  ,((1:nc)[gpr_rxn_ind] )),ncol=2)
       LHS[ii ]<-1
		
		if(gpr_rxn>1){
			diag(LHS[(nr+2)   :(nr+gpr_rxn+1)  ,((1:nc)[gpr_rxn_ind] )   ]) <- -1
			diag(LHS[(nr+gpr_rxn+2)   :(nr+2*gpr_rxn+1)  ,((1:nc)[gpr_rxn_ind] )    ]) <- 1

			diag(LHS[(nr+2)   :(nr+gpr_rxn+1)  ,(nc+1)  :(nc+1+gpr_rxn)]) <- -uppbnd(model)[gpr_rxn_ind]
			diag(LHS[(nr+gpr_rxn+2)   :(nr+2*gpr_rxn+1)  ,(nc+1)  :(nc+1+gpr_rxn) ]) <- -uppbnd(model)[gpr_rxn_ind]
		}else{# diag function fails when it is one row
				LHS[(nr+2)   :(nr+gpr_rxn+1)  ,((1:nc)[gpr_rxn_ind] ) ] = -1
				LHS[(nr+gpr_rxn+2)   :(nr+2*gpr_rxn+1)  ,((1:nc)[gpr_rxn_ind] )    ] <- 1
                LHS[(nr+2)   :(nr+gpr_rxn+1)  ,(nc+1)  :(nc+1+gpr_rxn)] <- -uppbnd(model)[gpr_rxn_ind]
			    LHS[(nr+gpr_rxn+2)   :(nr+2*gpr_rxn+1)  ,(nc+1)  :(nc+1+gpr_rxn) ] <- -uppbnd(model)[gpr_rxn_ind]
		}
	}
       # ---------------------------------------------
            # lower and upper bounds for COLUMNS(variables) and ROWS(constraints)
       # ---------------------------------------------
     
     	lower  <- c(lowbnd(model), rep(0, gpr_rxn), rep(0, ng))
     	upper  <- c(uppbnd(model), rep(1, gpr_rxn), rep(1, ng))
     	#RHS 1:nr+1 : as FBA model , FB,      ?? why -2*ub!!
		#rhs(model) : removed from model, 30/3/2013
     	#rlower <- c(rhs(model), FB, -2*uppbnd(model)[gpr_rxn_ind],-2*uppbnd(model)[gpr_rxn_ind],rep(0, gpr_rxn))  # ,rep(0,ngpr)
    	#rupper <- c(rhs(model), FB, rep(0, 3*gpr_rxn))  #,rep(0,ngpr)
        rlower <- c(rep(0,nr), FB, -2*uppbnd(model)[gpr_rxn_ind],-2*uppbnd(model)[gpr_rxn_ind],rep(0, gpr_rxn))  # ,rep(0,ngpr)
    	rupper <- c(rep(0,nr), FB, rep(0, 3*gpr_rxn))  #,rep(0,ngpr)
        
		gprtype=rep("E",gpr_rxn);
     # constraints of gpr: NOT is not considered
	message("start gpr....");
	 lastRow=dim(LHS)[1];

	 row_i=nr+2*gpr_rxn+1;

	for(i in 1:gpr_rxn){  # i : is the rxn number
		 rl=rules[i];
# 				# search for brackets : add aux variable for each bracket  in the same way as above
				# Consider only SUM of PRODUCT (AND to [OR]), without NOT, one level
		row_i=row_i+1;
		print(c(i,rl))# may make a problem if gene name contains 'or' like YOR123
		#replace )or with ) or & or( with or (
		rl=gsub("\\)"," ) ",rl)# 
		rl=gsub("\\("," ( ",rl)# 

		pr=lapply(strsplit(unlist(strsplit(rl," or "))," and "),function(x) gsub("[() ]","",x))
		if( length(pr)==1) {# no OR (only one term) 
			LHS[row_i,nc+gpr_rxn+which(geneNames %in% pr[[1]])]=1;
			LHS[row_i,nc+i]=-length(pr[[1]]);  # this is for the rxn
			#rlower[row_i]=-0.1; 
			rupper[row_i]=length(pr[[1]] )-1;#+0.1
			if( length(pr[[1]]) >1)  gprtype[i]="R"; # Range
		 } else {# first the sum row
		     LHS[row_i,nc+i]=length(pr);  # this is for the rxn
		     #rlower[row_i]=-0.1; # put to zero 
		     rupper[row_i]=length(pr)-1;	
		     gprtype[i]="R";
		     for( p in 1:length(pr)){
                     	 if( length(pr[[p]])==1 ){ 
			              LHS[row_i,nc+gpr_rxn+which( geneNames %in% pr[[p]] ) ]=-1;
			             }else{       # add auxiliary variable : add row for it then column update new row and SUM row
							  crow=Matrix::Matrix(0, nrow = 1, ncol =dim(LHS)[2])
							  LHS=rBind(LHS,crow)
							  rlower=c(rlower,0); 				  rupper=c(rupper,0);
							  lastRow=lastRow+1;
							  ccol=Matrix::Matrix(0, nrow = lastRow, ncol =1)
							  LHS=cBind(LHS,ccol) # new column
							  lower=c(lower,0);	  upper=c(upper,1);
							  LHS[lastRow,dim(LHS)[2]]=-length(pr[[p]]);  
							  LHS[lastRow,nc+gpr_rxn+which(geneNames %in% pr[[p]])]=1;
							  rupper[lastRow]=length(pr[[p]])-1;
							  LHS[row_i,dim(LHS)[2]]=-1;  # add the Aux variable to Sum row
							}
			}#for p
		}#if len
       }# for i

       message("gpr ... ... OK")
       # ---------------------------------------------
       # objective function
       # ---------------------------------------------
       # Notes: 1-use gene status
       # Min(ifelse(expr(rxn(i),-yi,yi)     math. |x-y|=x-y  when x>=y and y-x when y>=x 
       # then the difference will be Sum(rxn(i))+obj 
       #
       
        
    num_constr=dim(LHS)[1];
	num_var=dim(LHS)[2];
	aux_cnt=num_constr-(nr+3*gpr_rxn+1);
	#print(paste("no cons:",num_constr,"no var:",num_var));
    cobj <- c(rep(0, nc+gpr_rxn),ifelse(is.na(x),0,ifelse(x,-1,1)) ,rep(0, num_var-(nc+gpr_rxn+ng))) # NA->0 

        # ---------------------------------------------
        # build problem object
        # ---------------------------------------------
	switch(solver,
            # ----------------------- #
            "glpkAPI" = {
                out <- vector(mode = "list", length = 4)
                prob <- initProbGLPK();# new problem
               
                out[[1]] <-  addRowsGLPK(prob, nrows=num_constr)
		outj <- addColsGLPK(prob, ncols=num_var)
		setColNameGLPK(prob,c(1:(num_var)),paste(c(rep("x",nc),rep("rxn",gpr_rxn),rep("g",ng),rep("Aux",(num_var-(nc+gpr_rxn+ng)))),
		               c(1:nc,which(gpr_rxn_ind),1:ng,1:(num_var-(nc+gpr_rxn+ng))),sep="" ) );
		setObjDirGLPK(prob, GLP_MIN);
                ## note: when FX or LO value taken from lb and when UP take from UP and when R
               rtype <- c(rep(GLP_FX, nr+1), rep(GLP_UP, 2*gpr_rxn), ifelse(gprtype=="E",GLP_FX,GLP_DB),rep(GLP_FX,aux_cnt))
	       
	        # set the right hand side Sv = b

	       out[[4]] <- setRowsBndsGLPK(prob, c(1:num_constr), lb=rlower, ub=rupper,type=rtype )
        
	        # add upper and lower bounds: ai <= vi <= bi
                cc <- setColsBndsObjCoefsGLPK(prob, c(1:num_var),   lower,   upper,  cobj   )
                #print(cc);			
        	#nzLHS=nzijr(LHS);
            #cc <- loadMatrixGLPK(prob,   nzLHS$ne,   nzLHS$ia,  nzLHS$ja,  nzLHS$ar     )
			TMPmat <- as(LHS, "TsparseMatrix")
            cc <- loadMatrixGLPK(prob,length(TMPmat@x),ia  = TMPmat@i + 1,ja  = TMPmat@j + 1,ra  = TMPmat@x)
                
            ctype <- c(rep(GLP_CV, nc), rep(GLP_BV,num_var-nc)); # number of integer variable
		    setColsKindGLPK(prob,c(1:(num_var)),ctype);

                if (verboseMode > 2) {                      
				        fname=format(Sys.time(), "glpk_eFBA_%Y%m%d_%H%M.lp");
					print(sprintf("writing problem to: %s/%s...",getwd(),fname));
		                	writeLPGLPK(prob,fname);
		                	print("Solving...");
                }
                
                ## Solve
                setMIPParmGLPK(PRESOLVE,GLP_ON);
		lp_ok=solveMIPGLPK(prob);
		 if (verboseMode > 2) {
		 	print(return_codeGLPK(lp_ok));
		 	}
		lp_stat=mipStatusGLPK(prob);
		 if (verboseMode > 2) {
		 	print(status_codeGLPK(lp_stat));
		 	}
		lp_obj=mipObjValGLPK(prob);
		colst=mipColsValGLPK(prob);
		newFlux=colst
	       ## --- 
	        newFlux=sybil:::.floorValues(newFlux,tol=Tf);
	        sol_geneStat=ifelse(newFlux[(nc+gpr_rxn+1) : (nc+gpr_rxn+ng)]==1,"ON","OFF");
	        newFlux=newFlux[1:nc];
	        newStat=ifelse(abs(newFlux)>Tf,1,0);

            },
            # ----------------------- ----------------------- ------------------------- ----------------------- ----------------------#
           "cplexAPI" = {
                 out <- vector(mode = "list", length = 3)
				 prob <- openProbCPLEX()
				 out <- setIntParmCPLEX(prob$env, CPX_PARAM_SCRIND, CPX_OFF)
		                
                chgProbNameCPLEX(prob$env, prob$lp, "eFBA gene cplex");
                # when R: rngval+rhs   to rhs   (rngval<0)
                rtype <- c(rep("E",nr+1),  rep("L", 2*gpr_rxn),gprtype,rep("R",aux_cnt))
                prob$lp<- initProbCPLEX(prob$env)
  		chgProbNameCPLEX(prob$env, prob$lp, "eFBA");
                 out[[1]] <- newRowsCPLEX(prob$env, prob$lp,
                                         nrows=num_constr, rhs=rupper, sense=rtype,rngval=rlower-rupper)

                 cNames=paste(c(rep("x",nc),rep("rxn",gpr_rxn),rep("g",ng),rep("Aux",(num_var-(nc+gpr_rxn+ng)))),
		 		               c(1:nc,which(gpr_rxn_ind),1:ng,1:(num_var-(nc+gpr_rxn+ng))),sep="" ) ;
                 out[[2]] <- newColsCPLEX(prob$env, prob$lp,
                                         num_var, obj=cobj, lb=lower, ub=upper,cnames=cNames)
          
				TMPmat <- as(LHS, "TsparseMatrix")
				out[[3]] <- chgCoefListCPLEX(prob$env, prob$lp,#lp@oobj@env, lp@oobj@lp,
                                   nnz = length(TMPmat@x),
                                   ia  = TMPmat@i,
                                   ja  = TMPmat@j,
                                   ra  = TMPmat@x);				
		  #    nzLHS=nzijr(LHS);        
                # out[[3]] <- chgCoefListCPLEX(prob$env, prob$lp,
                                             # nzLHS$ne,
                                             # nzLHS$ia  - 1,
                                            # nzLHS$ja - 1,
                                             # nzLHS$ar)
          
                if(num_var>nc){ ctype<-c(rep('C', nc), rep('B',num_var - nc));# integer variables
				}else{ctype<-rep('C', nc);}
				status = copyColTypeCPLEX (prob$env, prob$lp, ctype);
				check <- setObjDirCPLEX(prob$env, prob$lp, CPX_MIN);
		
	         if (verboseMode > 2) {                      
                   fname=format(Sys.time(), "Cplex_eFBA_gene_%Y%m%d_%H%M.lp");
	 			    print(sprintf("writing problem to: %s/%s ...",getwd(),fname));
		 			writeProbCPLEX(prob$env, prob$lp,fname);
		       }
              #   --------------------------------------------------------------
              print("Solving...");
	       lp_ok     <- mipoptCPLEX(prob$env, prob$lp);
	       print(lp_ok);
	       sol=solutionCPLEX(prob$env, prob$lp);
                if (verboseMode > 2) {
                	print(sol)
                	}
              if(is(sol)[1]=="cpxerr"){
					print(sol);
					stop("Execution Terminated! error in MILP")
              }
              lp_obj=sol$objval;
               lp_stat   <- getStatCPLEX(prob$env, prob$lp)
              
               colst=sol$x;
                 
               #have flux and gene ON /flux and gene OFF?
               newFlux=sybil:::.floorValues(sol$x,tol=Tf);
               newFlux=newFlux[1:nc];
               sol_geneStat=ifelse(sol$x[(nc+gpr_rxn+1):(nc+gpr_rxn+ng)]==1,"ON","OFF");
               newStat=ifelse(abs(newFlux)>Tf,1,0);
               # when using binary variables as indicators: .floorValues(.ceilValues(newFlux[(nc+1):(2*nc)],tol=Tf/10),tol=Tf);
           },
           # ----------------------- #
           {
               wrong_type_msg(solver)
            }
        )
        
 	       print(sprintf("Number of Rxns with gene ON=%d ",length(rxnStatus[rxnStatus==-1]) ));

	       print(sprintf("Rxns with gene OFF=%d",length(rxnStatus[rxnStatus==1]) ));  
	       print(sprintf("Total Nr of selected rules=%d ",length(rxnStatus[rxnStatus==0]) ));
	       print(sprintf("Total number of rxns: %d",length(rxnStatus) ));
	       print(sprintf("The differnece is: %.0f ", (lp_obj+length(rxnStatus[rxnStatus==-1])) ));

 		origFlux=sybil:::.floorValues(fluxes(orig_sol)[fldind(orig_sol)],tol=Tf);
 		
 		#excReact = findExchReact(model)[1];# 1 is position
		#excReactPos=react_pos(excReact$exchange);
        excReact = findExchReact(model);
        excReactPos=(react_id(model) %in% react_id(excReact));#excReact$exchange

        is_excReact=((c(1:length(react_id(model)))) %in% excReactPos);
               
                origStat=ifelse(abs(origFlux)>Tf,1,0);
               rxnstname=ifelse(rxnStatus==0,"No rule/Unk",ifelse(rxnStatus==-1,"ON","OFF"));
               rxn_st <- cbind(gpr=gpr(model),react=react_id(model),reactName=react_name(model),is_excReact=is_excReact,
                                 expr=rxnstname,lb=lowbnd(model),
                                       #  eqns=sapply(c(1:length(react_id(model))),function(x)getRxnEqn(x)) ,
                                          origFlux,origStat, newFlux,newStat,difference=abs(origStat-newStat));
 
               gene_st <- cbind(geneLocus=allGenes(model),inputStat=x,SolStat=sol_geneStat)
#               write.csv(rxn_st,"rst.csv");
#              write.csv(gene_st,"gene.csv");
               
#               output_rep=table(list(expr=rxnstname,orig=origStat,newsol=newStat));
#               write.csv(output_rep,"eFBA_report.csv");

 
  remove(prob);
#----------------------------------------  ***  --------------------------------#
	optsol <- list(        ok = lp_ok,
			       obj = FB,
			       stat = lp_stat,
			       fluxes = newFlux,
			       origfluxes=origFlux,
			       rxn=rxn_st,
			       new_objFn=lp_obj,
			       gene_stat=sol_geneStat,
			       colst=colst
			  )
	return(optsol);
}



#if(solver=="cplex"){
#   writeProbCPLEX(prob$env, prob$lp,"e:\\sybil\\Ec_eFBA_f1.lp");
#   }
#   else{
#   writeLPGLPK(prob,"E:\\Sybil\\Ec_FBAGLPK.lp");
#   }
#----------------------------------------  ***  --------------------------------#
#check <- loadProblemDataEFBA(lp, model, FB, expressionData=expressionData, nCols=react_num(model)*2, 
#                                                 nRows= met_num(model)+1+react_num(model)*2);
# Test gpr
# rl="(yor123 and abcd)or yfl012 or (abc and (xyz or wxc))"
# rl=gsub("\\)"," ) ",rl)# 
# rl=gsub("\\("," ( ",rl)# 
# rl
# pr=lapply(strsplit(unlist(strsplit(rl," or "))," and "),function(x) gsub("[() ]","",x))
# pr
