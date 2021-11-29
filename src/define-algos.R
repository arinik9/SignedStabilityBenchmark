#################################################################
# Correlation Clustering (CC) problem
#################################################################


# ===============================================================
# Exact Approach: ExCC
# ===============================================================

NB.THREAD = 5

COR.CLU.ExCC <- "ExCC"
COR.CLU.ExCC.ENUM.ALL <- "ExCC-all" # TODO
ExCC.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.ExCC)
ExCC.JAR.PATH = paste(ExCC.LIB.FOLDER,"ExCC.jar",sep="/") # gaia cluster - CERI
#CPLEX.BIN.PATH = "/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/"
CPLEX.BIN.PATH = "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/"
ExCC.MAX.TIME.LIMIT = 3600

ExCCAll.MAX.NB.SOLS = 1000
ExCCAll.MAX.TIME.LIMIT = 3600*12 # 12 hours


# ===============================================================
# Exact Approach: EnumCC
# ===============================================================

ENUMCC = "EnumCC"
ENUMCC.LIB.FOLDER = file.path(LIB.FOLDER,ENUMCC)
ENUMCC.JAR.PATH = paste(ENUMCC.LIB.FOLDER,paste0(ENUMCC,".jar"),sep="/") # gaia cluster - CERI
ENUMCC.JAR.PATH = paste(ENUMCC.LIB.FOLDER,paste0("EnumCC.jar"),sep="/") # gaia cluster - CERI
#ENUM.POPULATE.CC.JAR.PATH = paste(ENUM.POPULATE.CC.LIB.FOLDER,paste0("MyPopulateCC_CplexFast_NEW.jar"),sep="/") # gaia cluster - CERI

ENUMCC.MAX.NB.SOLS = 50000 #10000 # we know that our method is not efficient for very large number of solutions
ENUMCC.MAX.TIME.LIMIT = 3600*12 # 12 hours


# ===============================================================
# Exact Approach: RNSCC (which is part of EnumCC)
# ===============================================================

CoNSCC = "CoNS"
RNSCC = "RNSCC"
RNSCC.JAR.PATH = paste(ENUMCC.LIB.FOLDER,paste0(RNSCC,".jar"),sep="/")



# ===============================================================
# Heuristic Approach: ILS and GRASP
# ===============================================================
CODE.COR.CLU.ILS = "ILS"
COR.CLU.ILS <- "ILS"
COR.CLU.ILS.CC <- paste(COR.CLU.ILS,"CC",sep="-")
CODE.COR.CLU.GRASP = "GRA"
COR.CLU.GRASP <- "GRASP"
COR.CLU.GRASP.CC <- paste(COR.CLU.GRASP,"CC",sep="-")

CODE.COR.CLU.VNS = "VNS" # Levorato et al
COR.CLU.VNS = "VNS"
COR.CLU.VNS.CC = paste(COR.CLU.VNS,"CC",sep="-")
COR.CLU.VNS.RCC = paste(COR.CLU.VNS,"RCC",sep="-")

CODE.COR.CLU.VOTE.BOEM = "V-B"
COR.CLU.VOTE.BOEM = "VOTE-BOEM"
COR.CLU.VNS.VOTE.BOEM.CC = paste(COR.CLU.VOTE.BOEM,"CC",sep="-")
COR.CLU.VNS.VOTE.BOEM.RCC = paste(COR.CLU.VOTE.BOEM,"RCC",sep="-")

ILS.GRASP.LIB.FOLDER = file.path(LIB.FOLDER,paste(COR.CLU.ILS,COR.CLU.GRASP,sep="-"))




# ===============================================================
# Heuristic Approach: NIFTY
# ===============================================================
CODE.COR.CLU.NIFTY = "NIFTY"
COR.CLU.NIFTY = "NIFTY"


# ===============================================================
# Heuristic MLMSB
# ===============================================================
CODE.COR.CLU.MLMSB = "MLMSB"
COR.CLU.MLMSB = "MLMSB"


# ===============================================================
# Heuristic Simulated Annealing
# ===============================================================
CODE.COR.CLU.SA.CC = "SA-CC"
COR.CLU.SA = "SA"
COR.CLU.SA.CC = "SA-CC"
SA.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.SA)
SA.CC.JAR.PATH = paste(SA.LIB.FOLDER,"SACC.jar",sep="/")

# ===============================================================
# Heuristic Tabu Search (Brusco & Doreian)
# ===============================================================
CODE.COR.CLU.TS.CC = "TS-CC"
COR.CLU.TS = "TS"
COR.CLU.TS.CC = "TS-CC"
TS.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.TS)
TS.CC.JAR.PATH =  paste(TS.LIB.FOLDER,"TSCC.jar",sep="/")

# ===============================================================
# Heuristic VNS (Brusco & Doreian)
# ===============================================================
CODE.COR.CLU.VNS.CC = "Brusco-VNS-CC"
COR.CLU.Brusco.VNS = "Brusco-VNS"
COR.CLU.Brusco.VNS.CC = "Brusco-VNS-CC"
Brusco.VNS.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.Brusco.VNS)
Brusco.VNS.CC.JAR.PATH =  paste(Brusco.VNS.LIB.FOLDER,"Brusco-VNSCC.jar",sep="/")


# ===============================================================
# LPMP: Message Passing, GAEC+KLj and ICP (Lange et al)
# ===============================================================

CODE.COR.CLU.ICP.GAEC.KLj.CC = "ICP-GAEC-KLj-CC"
ICP.GAEC.KLj.LIB.FOLDER = file.path(LIB.FOLDER,"LPMP")
ICP.GAEC.KLj.EXECUTABLE.PATH = paste(ICP.GAEC.KLj.LIB.FOLDER,"multicut_cycle_packing",sep="/")

CODE.COR.CLU.GAEC.KLj.CC = "GAEC-KLj-CC"
GAEC.KLj.EXECUTABLE.PATH = paste(ICP.GAEC.KLj.LIB.FOLDER,"gaec_kernighan_lin",sep="/")

CODE.COR.CLU.GAEC.KLj.CC = "GAEC-KLj-CC"
GAEC.KLj.EXECUTABLE.PATH = paste(ICP.GAEC.KLj.LIB.FOLDER,"gaec_kernighan_lin",sep="/")

CODE.COR.CLU.MP.GAEC.KLj.CC = "MP-GAEC-KLj-CC"
MP.GAEC.KLj.EXECUTABLE.PATH = paste(ICP.GAEC.KLj.LIB.FOLDER,"multicut_message_passing_odd_wheel",sep="/")




#############################################################################################
#
#############################################################################################
get.k.from.algo.name = function(token, part.folder){
	# example1: kFrom=(ALGO.NAME):k+1
	# example2: k=?
	k = NA
	
	tmp2 <- strsplit(x=token, split="=", fixed=TRUE)[[1]] # c("kFrom=ILS-RCC:k+1", "..")
	if(tmp2[1] == "k"){ # the 'k' is already provided by user
		k = tmp2[2]
	} else if(tmp2[1] == "kFrom"){ # the 'k' is not provided by user, but we get this info through the result of an algo
		tmp3 <- strsplit(x=tmp2[2], split=":", fixed=TRUE)[[1]] # c("ILS-RCC","k+1")
		list.chars = strsplit(tmp3[1],'')[[1]]
		# remove the paranthesis located in the begining and the end
		k.from.algo.name <- paste(list.chars[2:(length(list.chars)-1)], collapse="") # except the 1st and the last character which are paranthesis
		
		partition.file = file.path(part.folder,paste0(k.from.algo.name,"-membership.txt"))
		
		
		if(!file.exists(partition.file)){
			tlog("........Partition file ",partition.file," not found")
			return(-1)
		}
		else{
			partition <- as.matrix(read.table(partition.file))
			k = length(unique(partition))
		}
		
		if(grepl('\\+', tmp3[2])){ # if it contains the '+' sign ==> 2 possiblities: "k" or "k+<NUMBER>" like "k+2"
			tmp4 <- strsplit(x=tmp3[2], split="+", fixed=TRUE)[[1]][2]
			k.incr.nbr = as.integer(tmp4)
			k = k + k.incr.nbr
		}
	}
	
	return(k)
}


#############################################################################################
# This method is especially designed when the "k" parameter is needed for execution of the algo.
# But, it handles with the other case where there is not any "k" value needed
# There are 2 input options when user needs to provide "k":
# 1): kFrom=(ALGO.NAME):k+1			====> take the k value of an algo
# 2) k=<NUMBER>                     ====> take the provided k value
#############################################################################################
break.down.algo.code = function(algo.name){
	# break down the specified code (short name)
	# ================================================================================================================
	tmp=NA
	if(grepl('\\(', algo.name) && grepl('\\)', algo.name)){ # if "(" and ")" paranthesis symbols are used
		list.chars = strsplit(algo.name, "")[[1]] # get a list of characters of the string 'alog.name'
		start.paranthesis = which(list.chars=="(")
		end.paranthesis = which(list.chars==")")
		# substiture the content of the paranthesis part with the word 'REPLACE'.
		# Thus, it will be easy to split the 'algo.name' by the symbol "_"
		algo.name2 = paste(c(list.chars[1:(start.paranthesis-1)], "REPLACE", list.chars[(end.paranthesis+1):length(list.chars)]), collapse="")
		tmp <- strsplit(x=algo.name2, split="_", fixed=TRUE)[[1]]
		
		# find the index whose item contains "REPLACE"
		indx = which(grepl("REPLACE", tmp))
		tmp[indx] = gsub("REPLACE", paste(list.chars[start.paranthesis:end.paranthesis], collapse=""), tmp[indx]) 
		
	} else {
		tmp <- strsplit(x=algo.name, split="_", fixed=TRUE)[[1]]
	}
	
	# At the end, the result of the example "ILS-RCC_kFrom=(ILS-CC_l1_a1_g0_p3_t3_i10):k+1_l1_a1_g0_p3_t3_i10"
	#	[1] "ILS-RCC"                              
	#	[2] "kFrom=(ILS-CC_l1_a1_g0_p3_t3_i10):k+1"
	#	[3] "l1"                                   
	#	[4] "a1"                                   
	#	[5] "g0"                                   
	#	[6] "p3"                                   
	#	[7] "t3"                                   
	#	[8] "i10"
	# ================================================================================================================
	
	return(tmp) # list of parameters
}









############################################################################
# It reads a .G graph file, and returns the contents as a data frame object.
#
# network.path: the file path which stores the .G graph file
#
############################################################################
read.graph.ils.file.as.df = function(network.path){
	# skip the first line bc it does not contain graph info
	df = read.table(
			file=network.path, 
			header=FALSE, 
			sep="\t", 
			skip=1, 
			check.names=FALSE
	)
	# df$V1: vertex1
	# df$V2: vertex2
	# df$V3: weight
	return(df)
}



############################################################################
#  It reads a .G graph file, and returns the contents as a igraph graph object.
#  To handle isolated nodes, first we had to find the max vertex id.
#  Then, indicate explicitely vertices ids in graph.data.frame()
#
# network.path: the file path which stores the .G graph file
#
############################################################################
read.graph.ils = function(network.path){
	df = read.graph.ils.file.as.df(network.path)
	
	edg.list = df[,c(1, 2)]
	max.v.id = max(unique(c(edg.list[,1], edg.list[,2])))
	
	g <- graph.data.frame(edg.list, vertices=seq(0,max.v.id), directed=FALSE)
	cat("max id: ",max.v.id, "\n")
	E(g)$weight = df[, 3]
	# V(g)$id = seq(0,max.v.id)
	# V(g)$id = seq(1,max.v.id+1)
	
	return(g)
}


############################################################################
# It writes the graph object into a file
#
# graph: the igraph graph object
# file path: the file path which will store the graph content in the format of .G graph file
############################################################################
write.graph.ils = function(graph, file.path){
# export using a format compatible with pILS
	t <-get.edgelist(graph=graph)
	t =  matrix(as.integer(t), nrow(t), ncol(t))
	if(t[1,1] == 1)
		t <- t - 1	# start numbering nodes at zero
	
	t <- cbind(t,E(graph)$weight)		# add weights as the third column
	
	write.table(data.frame(vcount(graph),ecount(graph)), file=file.path, append=FALSE, sep="\t", row.names=FALSE, col.names=FALSE) # write header
	write.table(t, file=file.path, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE) # write proper graph
}









#############################################################################################
# 
#############################################################################################
get.ExCC.code <- function(enum.all)
{
	result <- COR.CLU.ExCC
	if(enum.all)
		result <- COR.CLU.ExCC.ENUM.ALL
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.ExCC.code <- function(enum.all, formulation.type="", triangle.constr.reduced.ILP.model=FALSE)
{
	result <- COR.CLU.ExCC
	if(enum.all)
		result <- COR.CLU.ExCC.ENUM.ALL
	
	if(formulation.type != ""){
		result <- paste0(result,"-",formulation.type)
		
		if(formulation.type == "vertex" && triangle.constr.reduced.ILP.model == TRUE)
			result <- paste0(result,"_","reduced")
		else if(formulation.type == "vertex" && triangle.constr.reduced.ILP.model == FALSE)
			result <- paste0(result,"_","full")
		
	}
	return(result)
}





#############################################################################################
# 
#############################################################################################
get.ExCC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	is.cp = "true" # in any case, use cutting plane approach
	is.enum.all = "false"
	# tilim = 3600 # 1 hour
	formulation.type = ""
	ILP.form = "full"
	triangleIneqReducedForm = "false"
	lazyCB = "false"
	
	if(length(strsplit(x=algo.name, split="-", fixed=TRUE)[[1]])>1){
		base.algo.name <- strsplit(x=algo.name, split="-", fixed=TRUE)[[1]][1]
		str.params <- gsub(paste0(base.algo.name,"-"),"",algo.name)
		
		if(base.algo.name == COR.CLU.ExCC.ENUM.ALL){
			is.enum.all = "true"
		}
		
		if(length(strsplit(x=str.params, split="_", fixed=TRUE)[[1]])>1){ # Fv: ILP based on vertex-pair
			formulation.type <- strsplit(x=str.params, split="_", fixed=TRUE)[[1]][1]
			ILP.form = strsplit(x=str.params, split="_", fixed=TRUE)[[1]][2]
			
			if(ILP.form == "full")
				triangleIneqReducedForm = "false"
			else if(ILP.form == "reduced")
				triangleIneqReducedForm = "true"
			
		} else { # Fe: ILP based on edges
			formulation.type <- str.params
			if(formulation.type == "edge"){
				lazyCB = "true"
			}
		}
		
	}
	else {
		if(algo.name == COR.CLU.ExCC.ENUM.ALL){
			is.enum.all = "true"
		}
	}
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	input.file.for.g = file.path(input.folder, graph.name)
	g = read.graph.ils(input.file.for.g)
	
	
	isReducedForm = FALSE
	if(triangleIneqReducedForm == "true")
		isReducedForm = TRUE
	ExCC.folder = get.ExCC.code(FALSE,formulation.type,isReducedForm)
	
	initSolutionFilePath = file.path(out.folder,"..","..",ExCC.folder,SIGNED.UNWEIGHTED.FILE,"membership0.txt")
	if(!file.exists(initSolutionFilePath))
		initSolutionFilePath="''"
	
	cmd = "NONE"
	if(is.enum.all == "false"){
		# An example:
		# java -Djava.library.path=/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/
		# -DinFile="in/""$name" -DoutDir="out/""$modifiedName" -DenumAll=false -Dcp=true -DMaxTimeForRelaxationImprovement=20
		# -DuserCutInBB=false -DinitSolutionFilePath="$initSolutionFilePath" -DLPFilePath="$LPFilePath"
		# -DonlyFractionalSolution=false -DfractionalSolutionGapPropValue=-1.0 -DnbThread=2 -Dverbose=true -Dtilim=200 -jar exe/ExCC.jar
		
		# -------------------------------------------------------------------------
		# TODO handle this in a better way, for instance, use small values for sparse networks
		# TODO write a function called 'estimate.max.time.for.relaxation.improvement(..)'
		maxTimeForRelaxationImprovement = "600"
		#if(vcount(g)>39 && vcount(g)<=50)
		#	maxTimeForRelaxationImprovement = "4500" # 1h15m
		#else if(vcount(g)>50)
		#	maxTimeForRelaxationImprovement = "10000" # 166 mins 
		# -------------------------------------------------------------------------
		
		cmd = 
			paste(
				"java",		
				paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
				paste0("-DformulationType='",formulation.type,"'"),
				paste0("-DinFile=", input.file),
				paste0("-DoutDir=", out.folder),
				paste0("-Dcp=",is.cp),
				paste0("-DenumAll=",is.enum.all),
				paste0("-Dtilim=",ExCCAll.MAX.TIME.LIMIT),
				paste0("-DtilimForEnumAll=",-1),
				paste0("-DsolLim=",1),
				paste0("-DMaxTimeForRelaxationImprovement=",maxTimeForRelaxationImprovement),
				paste0("-DlazyCB=",lazyCB),
				"-DuserCutCB=false",
				paste0("-DinitMembershipFilePath=",initSolutionFilePath),
				"-Dverbose=true",
				paste0("-DnbThread=", NB.THREAD),
				"-DLPFilePath=''",
				"-DonlyFractionalSolution=false",
				"-DfractionalSolutionGapPropValue=0.01",
				paste0("-DtriangleIneqReducedForm=",triangleIneqReducedForm),
				"-jar",
				ExCC.JAR.PATH,
				sep=" "
			)
		
	} else { # if(is.enum.all == "true")
		# An example:
		# java -DinFile="in/""$name" -DoutDir="out/""$modifiedName" -DenumAll=true
		# -Dcp=false -DinitSolutionFilePath="$initSolutionFilePath" -DLPFilePath="$LPFilePath"
		# -DnbThread=2 -Dverbose=true -Dtilim=-1 -DtilimForEnumAll=60 -DsolLim=100 -jar exe/ExCC.jar
	
		
		LP.filepath = file.path(out.folder,"..","..",ExCC.folder,SIGNED.UNWEIGHTED.FILE,paste0("strengthedModel_",formulation.type,".lp"))
		print(LP.filepath)
		
		if(file.exists(LP.filepath)) {
			cmd = 
				paste(
					"java",		
					paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
					paste0("-DformulationType='",formulation.type,"'"),
					paste0("-DinFile=", input.file),
					paste0("-DoutDir=", out.folder),
					paste0("-Dcp=","false"),
					paste0("-DenumAll=","true"),
					paste0("-Dtilim=",-1),
					paste0("-DtilimForEnumAll=",ExCCAll.MAX.TIME.LIMIT), # We take this time limit account when OneTreeCC found a first solution
					paste0("-DsolLim=",ExCCAll.MAX.NB.SOLS),
					paste0("-DMaxTimeForRelaxationImprovement=","-1"), # no specific time limit for this phase, use the default one
					paste0("-DlazyCB=",lazyCB),
					"-DuserCutCB=false",
					paste0("-DinitMembershipFilePath=",initSolutionFilePath),
					"-Dverbose=true",
					paste0("-DnbThread=", NB.THREAD),
					paste0("-DLPFilePath='",LP.filepath,"'"),
					"-DfractionalSolutionGapPropValue=0.01",
					"-DonlyFractionalSolution=false",
					paste0("-DtriangleIneqReducedForm=",triangleIneqReducedForm),
					"-jar",
					ExCC.JAR.PATH,
					sep=" "
				)
		}
	}
	
	print(cmd)
	return(cmd)
}




#############################################################################################
# 
#############################################################################################
get.EnumCC.code <- function(maxNbEdit, formulation.type="", triangle.constr.reduced.ILP.model=FALSE)
{
	result <- paste0(ENUMCC,"-maxNbEdit",maxNbEdit)
	
	if(formulation.type != ""){
		result <- paste0(result,"-",formulation.type)
		
		if(formulation.type == "vertex" && triangle.constr.reduced.ILP.model == TRUE)
			result <- paste0(result,"_","reduced")
		else if(formulation.type == "vertex" && triangle.constr.reduced.ILP.model == FALSE)
			result <- paste0(result,"_","full")
	}
	
	return(result)
}



#############################################################################################
# 
#############################################################################################
get.EnumCC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	print(algo.name)
	base.algo.name <- strsplit(x=algo.name, split="-", fixed=TRUE)[[1]][1]
	remaining.algo.name <- gsub(paste0(base.algo.name,"-"),"",algo.name)
	#print(remaining.algo.name)
	
	base.algo.name <- strsplit(x=remaining.algo.name, split="-", fixed=TRUE)[[1]][1]
	maxNbEdit = as.integer(gsub("maxNbEdit","",base.algo.name))
	
	str.params <- strsplit(x=remaining.algo.name, split="-", fixed=TRUE)[[1]][2]
	#print(str.params)
	
	
	formulation.type = ""
	ILP.form = "full"
	triangleIneqReducedForm = "false"
	lazyCB = "false"
		
	if(length(strsplit(x=str.params, split="_", fixed=TRUE)[[1]])>1){ # Fv: ILP based on vertex-pair
		formulation.type <- strsplit(x=str.params, split="_", fixed=TRUE)[[1]][1]
		ILP.form = strsplit(x=str.params, split="_", fixed=TRUE)[[1]][2]
		
		if(ILP.form == "full")
			triangleIneqReducedForm = "false"
		else if(ILP.form == "reduced")
			triangleIneqReducedForm = "true"
		
	} else { # Fe: ILP based on edges
		formulation.type <- str.params
		if(formulation.type == "edge"){
			lazyCB = "true"
		}
	}
		
	
	isReducedForm = FALSE
	if(triangleIneqReducedForm == "true")
		isReducedForm = TRUE
	ExCC.folder = get.ExCC.code(FALSE,formulation.type,isReducedForm)
	
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	
	cmd = "NONE"
	LP.filepath = file.path(out.folder,"..","..",ExCC.folder,SIGNED.UNWEIGHTED.FILE,paste0("strengthedModel_",formulation.type,".lp"))
	print(LP.filepath)
	if(file.exists(LP.filepath)) {
		cmd = 
				paste(
						"java",		
						paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
						paste0("-DinFile=", input.file),
						paste0("-DoutDir=", out.folder),
						paste0("-DLPFilePath=", LP.filepath),
						paste0("-DinitMembershipFilePath=", file.path(out.folder,"..","..",ExCC.folder,SIGNED.UNWEIGHTED.FILE,"membership0.txt")),
						paste0("-DnbThread=",NB.THREAD),
						paste0("-DmaxNbEdit=",maxNbEdit),
						paste0("-DsolLim=",ENUMCC.MAX.NB.SOLS),
						paste0("-Dtilim=",ENUMCC.MAX.TIME.LIMIT),
						paste0("-DJAR_filepath_RNSCC=",paste(ENUMCC.LIB.FOLDER,paste0("RNSCC.jar"),sep="/")),
						paste0("-DtriangleIneqReducedForm=",triangleIneqReducedForm),
						paste0("-DlazyCB=",lazyCB),
						paste0("-DuserCutCB=","false"),
						paste0("-DformulationType=",formulation.type),
						"-jar",
						ENUMCC.JAR.PATH,
						sep=" "
				)
	}
	
	print(cmd)
	return(cmd)
}


#############################################################################################
# 
#############################################################################################
get.CoNSCC.code <- function(maxNbEdit, pruning.without.MVMO)
{
	result <- paste0(CoNSCC,"-maxNbEdit",maxNbEdit)
	if(pruning.without.MVMO)
		result <- paste0(result,"-pruningWithoutMVMO")
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.RNSCC.code <- function(maxNbEdit, pruning.without.MVMO)
{
	result <- paste0(RNSCC,"-maxNbEdit",maxNbEdit)
	if(pruning.without.MVMO)
		result <- paste0(result,"-pruningWithoutMVMO")
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.RNSCC.command <- function(algo.name, sol.lim, input.folder, out.folder, graph.name)
{
	print(algo.name)
	base.algo.name <- strsplit(x=algo.name, split="-", fixed=TRUE)[[1]][1]
	params.str <- gsub(paste0(base.algo.name,"-"),"",algo.name)
	print(params.str)
	
	pruningWithoutMVMO = "false"
	maxNbEdit = 1
	
	params.str.list <- unlist(strsplit(x=params.str, split="-", fixed=TRUE))
	if(length(params.str.list) == 1)
		maxNbEdit = as.integer(gsub("maxNbEdit","",params.str))
	else #  length(params.str.list) == 2
	{
		maxNbEdit = as.integer(gsub("maxNbEdit","",params.str.list[1]))
		pruningWithoutMVMO = "true"
	}
	
	graph.name = paste0("signed-unweighted", ".G")
	network.path = file.path(input.folder,graph.name)
	
	cmd = 
		paste(
			"java",
			paste0("-DinitMembershipFilePath=", file.path(out.folder,paste0(MBRSHP.FILE.PREFIX,"0.txt"))),
			paste0("-DallPreviousResultsFilePath=", file.path(out.folder,"allResults.txt")), # it should be an empty file at startup
			paste0("-DinputFilePath=", network.path),
			paste0("-DoutDir=", out.folder),
			paste0("-DmaxNbEdit=", maxNbEdit),
			"-Dtilim=-1",
			paste0("-DsolLim=",sol.lim),
			"-DnbThread=1",
			paste0("-DisBruteForce=",pruningWithoutMVMO),
			"-DisIncrementalEditBFS=true",
			"-jar",
			RNSCC.JAR.PATH,
			sep=" "
		)
	
	print(cmd)
	return(cmd)
}






#############################################################################################
# Returns the code (short name) for the Iterated Local Search (ILS) partioning method. See 
# the algorithm documentation for more details.
#
# rcc: whether to solve the correlation clustering (FALSE) or relaxed CC problem (TRUE).
# l: neighborhood size (during the local search).
# k: number of clusters (max number for RCC)
# k.from: the algo name from which we use the 'k' value. for RCC problem
# rel.k.val: the relative 'k' value. It is used with "k.from" parameter. Ex: "k", "k+2", etc.
# alpha: randomness factor.
# gain: 0 for min imbalance, 
#       1 for max modularity gain function, 
#       2 for max negative modularity gain function, 
#       3 for max positive-negative modularity gain function, 
#       4 for max pos-neg mod gain function II, 
#       5 for max pos-neg mod gain function III
# perturbation: maximal level of perturbation.
# vns: enbale VNS metaheuristic
#
# returns: the short name corresponding to the ILS method with the specified parameters.
#############################################################################################
get.ils.code <- function(l, alpha, rcc, gain, perturbation, time.limit=3600, iter.nbr=10,
		vns=FALSE, k=NA, k.from=NA, rel.k.val=NA)
{	result <- NA
	
	if(rcc){
		# example1: ILS-RCC_kFrom=ILS-RCC:k+1_....
		# example2: ILS-RCC_k=2_....
		
		if(vns){
			result <- COR.CLU.VNS.RCC
		} else {
			
			result <- COR.CLU.ILS.RCC
			if(is.na(k))
				result <- paste0(result,"_kFrom=(",k.from,"):",rel.k.val)
			else
				result <- paste0(result,"_k=",k)
		}
	}
	else{
		if(vns){
			result <- COR.CLU.VNS.CC
		} else {
			result <- COR.CLU.ILS.CC
		}
		
	}
	
	
	result <- paste0(result,"_l",l)
	result <- paste0(result,"_a",alpha)
	result <- paste0(result,"_g",gain)
	result <- paste0(result,"_p",perturbation)
	result <- paste0(result,"_t",time.limit)
	result <- paste0(result,"_i",iter.nbr)
	
	return(result)
}






#############################################################################################
# Returns the inline command for the Iterated Local Search (ILS) partioning method. See 
# the algorithm documentation for more details.
#
# algo.name: short code associated to the algorithm.
# input.folder: relative path to the folder containing targeted graph file.
# out.folder: relative path to the folder in which the produced files will be placed.
# k: nb cluster to be detected for RCC (when RCC is enbaled)
#
# returns: the command allowing to invoke the program externally.
#############################################################################################
get.ils.command <- function(algo.name, input.folder, out.folder, graph.name)
{	
	# break down the specified code (short name)
	tmp = break.down.algo.code(algo.name)
	
	
	params <- c()
	k=NA
	next.param.indx=NA
	# =======================================================================================
	# rcc flag
	# example1: ILS-RCC_kFrom=(ILS-CC_l1_a1_g0_p3_t3_i10):k+1_l1_a1_g0_p3_t3_i10
	# example2: ILS-RCC_k=2_l1_a1_g0_p3_t3_i10
	algo.name <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][1]
	rcc.flag <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][2]
	if(rcc.flag=="RCC"){
		params["rcc"] <- 1
		k = get.k.from.algo.name(tmp[2], out.folder)
		next.param.indx = 3
	}
	else{
		params["rcc"] <- 0
		next.param.indx = 2
	}
	
	# VNS meta heuristic
	if(algo.name == COR.CLU.VNS)
		params["vns"]=1
	else if(algo.name == COR.CLU.ILS)
		params["vns"]=0
	
	# =======================================================================================
	
	
	
	for(s in tmp[next.param.indx:length(tmp)])
	{	params <- c(params,substr(s,2,nchar(s)))
		names(params)[length(params)] <- substr(s,1,1)
	}
	
	
	# init
	input.file <- file.path(input.folder, graph.name)
	output.folder <- file.path(out.folder)
	command.folder <- file.path(ILS.GRASP.LIB.FOLDER)
	result <- file.path(command.folder, "graspcc")
	#result <- paste0("mpirun -n 1 ",result)
	
	
	
	# build command
	result <- paste0(result, " --vns ",params["vns"])
	result <- paste0(result, " --rcc ",params["rcc"])
	if(rcc.flag=="RCC")
		result <- paste0(result, " --k ",k)
	result <- paste0(result, " --alpha ",params["a"])
	result <- paste0(result, " --iterations ",params["i"])
	result <- paste0(result, " --neighborhood_size ",params["l"])
	result <- paste0(result, " --time-limit ",params["t"])
	result <- paste0(result, " --input-file \"",input.file,"\"")
	result <- paste0(result, " --output-folder \"",output.folder,"\"")
	result <- paste0(result, " --gain-function-type ",params["g"])
	result <- paste0(result, " --strategy ","ILS") 
	result <- paste0(result, " --perturbationLevelMax ",params["p"])
	
	return(result)
}


#############################################################################################
# Returns the code (short name) for the Grasp partioning method. See the algorithm documentation
# for more details.
#
# rcc: whether to solve the correlation clustering (FALSE) or relaxed CC problem (TRUE).
# l: neighborhood size (during the local search).
# k: number of clusters (max number for RCC)
# alpha: randomness factor.
# gain: 0 for min imbalance, 
#       1 for max modularity gain function, 
#       2 for max negative modularity gain function, 
#       3 for max positive-negative modularity gain function, 
#       4 for max pos-neg mod gain function II, 
#       5 for max pos-neg mod gain function III
# perturbation: maximal level of perturbation.
#
# returns: the short name corresponding to the Grasp method with the specified parameters.
#############################################################################################
get.grasp.code <- function(rcc, l, k, alpha, gain, time.limit, iter.nbr, oneOptNeig)
{	
	# TODO k can be specified for GRASP as well?
	
	result <- COR.CLU.GRASP
	if(!oneOptNeig && iter.nbr == -1)
		result <- COR.CLU.VOTE.BOEM
	
	if(rcc)
		result <- paste0(result,"-RCC")
	else
		result <- paste0(result,"-CC")
	
	
	result <- paste0(result,"_l",l)
	result <- paste0(result,"_a",alpha)
	result <- paste0(result,"_g",gain)
	result <- paste0(result,"_t",time.limit)
	result <- paste0(result,"_i",iter.nbr)
	result <- paste0(result,"_n", oneOptNeig)
	
	return(result)
}


#############################################################################################
# Returns the inline command for the Greedy Randomized Adaptive Search Procedure (Grasp) partioning 
# method. See the algorithm documentation for more details.
#
# algo.name: short code associated to the algorithm.
# input.folder: relative path to the folder containing targeted graph file.
# out.folder: relative path to the folder in which the produced files will be placed.
# time.limit: maximum duration of the processing.
# iter.nbr: maximum number of iterations of the processing.
#
# returns: the command allowing to invoke the program externally.
#############################################################################################
get.grasp.command <- function(algo.name, input.folder, out.folder, graph.name)
{	
	# break down the specified code (short name)
	tmp = break.down.algo.code(algo.name)
	
	params <- c()
	k=NA
	next.param.indx=NA
	# =======================================================================================
	# rcc flag
	# example1: GRASP-RCC_kFrom=(ILS-CC_l1_a1_g0_p3_t3_i10):k+1_l1_a1_g0_p3_t3_i10
	# example2: GRASP-RCC_k=2_l1_a1_g0_p3_t3_i10
	algo.name <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][1]
	rcc.flag <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][2]
	if(rcc.flag=="RCC"){
		params["rcc"] <- 1
		k = get.k.from.algo.name(tmp[2], out.folder)
		next.param.indx = 3
	}
	else{
		params["rcc"] <- 0
		next.param.indx = 2
	}
	# =======================================================================================
	
	for(s in tmp[next.param.indx:length(tmp)])
	{	params <- c(params,substr(s,2,nchar(s)))
		names(params)[length(params)] <- substr(s,1,1)
	}
	
	
	# init
	input.file <- file.path(input.folder, graph.name)
	output.folder <- file.path(out.folder)
	command.folder <- file.path(ILS.GRASP.LIB.FOLDER)
	result <- file.path(command.folder, "graspcc")
	# result <- paste0("mpirun -n 1 ",result)
	
	
	# build command
	result <- paste0(result, " --alpha ",params["a"])
	result <- paste0(result, " --iterations ",params["i"])
	result <- paste0(result, " --neighborhood_size ",params["l"])
	result <- paste0(result, " --rcc ",params["rcc"])
	if(rcc.flag=="RCC")
		result <- paste0(result, " --k ",k)
	
	result <- paste0(result, " --time-limit ",params["t"])
	result <- paste0(result, " --input-file \"",input.file,"\"")
	result <- paste0(result, " --output-folder \"",output.folder,"\"")
	result <- paste0(result, " --gain-function-type ",params["g"])
	result <- paste0(result, " --strategy ","GRASP")
	result <- paste0(result, " --firstImprovementOnOneNeig ",params["n"])
	
	return(result)
}



#############################################################################################
# Currently, there is not any available implementation for RCC problem
#
# linkage.criteria: criteria from which we decide if we merge 2 clusters or not during the method
# add.cannot.link.constraints: it is not related to CC. Basically, if an negative large weighted edge exists, 
#                               the nodes sharing this edge can not be in the same cluster
#############################################################################################
get.NIFTY.code <- function(method)
{
	result <- COR.CLU.NIFTY
	result <- paste0(result,"-CC")
	result <- paste0(result,"_",method)
	
	return(result)
}

#############################################################################################
# 
#############################################################################################
get.NIFTY.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	graph.name = paste0(unlist(strsplit(graph.name,"\\."))[1],".graphml")
	method = unlist(strsplit(algo.name,"_"))[2]
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	
	# /opt/anaconda3/envs/GASP/bin/python run_NIFTY.py 
	#   "/home/nejat/eclipse/workspace-neon/CoH-SosoCC/in/random-networks/n=20_l0=4_dens=0.1250/propMispl=0.1000/propNeg=0.2000/network=1"
	#   "signed-unweighted.graphml" "out" --method "greedy-additive + cgc-qpbo"
	
	cmd = 
			paste(
					"/opt/anaconda3/envs/GASP/bin/python",	# do not care about the name of the anaconda env (it took time to make somethng work in their env..)
					paste("lib/NIFTY/run_NIFTY.py"),
					paste0("'",input.folder,"'"),
					paste0("'",graph.name,"'"),
					paste0("'",out.folder,"'"),
					paste0("--method ","'",method,"'"),
					sep=" "
			)
	
	print(cmd)
	
	return(cmd)
}



#############################################################################################
# Currently, there is not any available implementation for RCC problem
#
# linkage.criteria: criteria from which we decide if we merge 2 clusters or not during the method
# add.cannot.link.constraints: it is not related to CC. Basically, if an negative large weighted edge exists, 
#                               the nodes sharing this edge can not be in the same cluster
#############################################################################################
get.MLMSB.code <- function(trial, time.limit)
{
	# trial: independent runs
	result <- COR.CLU.MLMSB
	result <- paste0(result,"-CC")
	result <- paste0(result,"_r",trial)
	result <- paste0(result,"_t",time.limit)
	
	return(result)
}




#############################################################################################
# Currently, there is not any available implementation for RCC problem
#
#############################################################################################
get.SA.CC.code <- function(time.limit)
{
	result <- COR.CLU.SA.CC
	#result <- paste0(result,"-CC")
	result <- paste0(result,"_t",time.limit)
	
	return(result)
}


#############################################################################################
# Currently, there is not any available implementation for RCC problem
#
#############################################################################################
get.TS.CC.code <- function(time.limit)
{
	result <- COR.CLU.TS.CC
	#result <- paste0(result,"-CC")
	result <- paste0(result,"_t",time.limit)
	
	return(result)
}


#############################################################################################
# Currently, there is not any available implementation for RCC problem
#
#############################################################################################
get.Brusco.VNS.CC.code <- function(time.limit)
{
	result <- COR.CLU.Brusco.VNS.CC
	#result <- paste0(result,"-CC")
	result <- paste0(result,"_t",time.limit)
	
	return(result)
}



#############################################################################################
# Currently, there is not any available implementation for RCC problem
#
#############################################################################################
get.ICP.GAEC.KLj.CC.code <- function(time.limit)
{
	result <- CODE.COR.CLU.ICP.GAEC.KLj.CC
	#result <- paste0(result,"-CC")
	result <- paste0(result,"_t",time.limit)
	
	return(result)
}



#############################################################################################
# Currently, there is not any available implementation for RCC problem
#
#############################################################################################
get.GAEC.KLj.CC.code <- function(time.limit)
{
	result <- CODE.COR.CLU.GAEC.KLj.CC
	#result <- paste0(result,"-CC")
	result <- paste0(result,"_t",time.limit)
	
	return(result)
}



#############################################################################################
# Currently, there is not any available implementation for RCC problem
#
#############################################################################################
get.MP.GAEC.KLj.CC.code <- function(time.limit)
{
	result <- CODE.COR.CLU.MP.GAEC.KLj.CC
	#result <- paste0(result,"-CC")
	result <- paste0(result,"_t",time.limit)
	
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.MLMSB.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	l = unlist(strsplit(algo.name,"_"))
	tmp = l[2]
	trial = as.numeric(gsub("r","",tmp))
	tmp = l[3]
	tilim = as.numeric(gsub("t","",tmp))
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	# matlab  -nosplash -nodesktop -r "Trial=10;inputFilePath='n=24_l0=4_dens=0.1250/propMispl=0.2000/propNeg=0.3000/network=1/signed-unweighted.G';
	#                                   outputDir='.'; Memetic; exit;"
	cmd = 
			paste(
					"/media/nejat/8AD01AB7D01AAA09/MATLAB/bin/matlab -nosplash -nodesktop -r",	# do not care about the name of the anaconda env (it took time to make somethng work in their env..)
					paste0('"Trial=',trial,';'),
					paste0("inputFilePath=",input.file,";"),
					paste0("outputDir='",out.folder,"';"),
					paste0("addpath('",LIB.FOLDER,"');"),
					paste0("addpath('",file.path(LIB.FOLDER,"MLMSB"),"');"),
					paste0("TILIM=",tilim,";"),
					paste0('Memetic; exit;"'),
					sep=" "
			)
	
	print(cmd)
	
	return(cmd)
}


#############################################################################################
# 
#############################################################################################
get.SA.CC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	# -DinputFilePath="in/EGFR_symm.G" -DoutDir="out/EGFR_symm-SA" -jar lib/SACC.jar"
	cmd = 
			paste(
					"java",		
					paste0("-DinputFilePath=", input.file),
					paste0("-DoutDir=", out.folder),
					"-jar",
					SA.CC.JAR.PATH,
					sep=" "
			)
	
	print(cmd)
	
	return(cmd)
}



#############################################################################################
# 
#############################################################################################
get.TS.CC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	tmp = unlist(strsplit(algo.name,"_"))[2]
	tilim = as.numeric(gsub("t","",tmp))
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	# -DinputFilePath="in/EGFR_symm.G" -DoutDir="out/EGFR_symm-TS" -Dtilim=100 -jar lib/TSCC.jar"
	cmd = 
			paste(
					"java -Xmx1024m",		
					paste0("-DinputFilePath=", input.file),
					paste0("-DoutDir=", out.folder),
					paste0("-Dtilim=",tilim),
					"-jar",
					TS.CC.JAR.PATH,
					sep=" "
			)
	
	print(cmd)
	
	return(cmd)
}


#############################################################################################
# 
#############################################################################################
get.Brusco.VNS.CC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	tmp = unlist(strsplit(algo.name,"_"))[2]
	tilim = as.numeric(gsub("t","",tmp))
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	# -DinputFilePath="in/EGFR_symm.G" -DoutDir="out/EGFR_symm-VNS" -Dtilim=100 -jar lib/Brusco-VNSCC.jar"
	cmd = 
			paste(
					"java -Xmx1024m",		
					paste0("-DinputFilePath=", input.file),
					paste0("-DoutDir=", out.folder),
					paste0("-Dtilim=",tilim),
					"-jar",
					Brusco.VNS.CC.JAR.PATH,
					sep=" "
			)
	
	print(cmd)
	
	return(cmd)
}




#############################################################################################
# 
#############################################################################################
get.GAEC.KLj.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	tmp = unlist(strsplit(algo.name,"_"))[2]
	tilim = as.numeric(gsub("t","",tmp))
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	temp.file = paste("'", input.folder, "/", "temp.txt", "'", sep="")
	
	
	# ./gaec_kernighan_lin n40.txt res1.txt
	cmd = 
			paste(
					GAEC.KLj.EXECUTABLE.PATH,		
					temp.file,
					file.path(out.folder, paste0(ALGO.RESULT.FILE.PREFIX,0,".txt")),
					tilim,
					sep=" "
			)
	
	
	
	print(cmd)
	
	return(cmd)
}



#############################################################################################
# 
#############################################################################################
get.ICP.GAEC.KLj.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	tmp = unlist(strsplit(algo.name,"_"))[2]
	tilim = as.numeric(gsub("t","",tmp))
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	# ./test_multicut_cycle_packing n40.txt res1.txt 0.5
	cmd = 
			paste(
					ICP.GAEC.KLj.EXECUTABLE.PATH,		
					input.file,
					file.path(out.folder, paste0(ALGO.RESULT.FILE.PREFIX,0,".txt")),
					"0.5",
					tilim,
					sep=" "
			)
	
	print(cmd)
	
	return(cmd)
}


#############################################################################################
# 
#############################################################################################
get.MP.GAEC.KLj.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	tmp = unlist(strsplit(algo.name,"_"))[2]
	tilim = as.numeric(gsub("t","",tmp))
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	#  ./multicut_message_passing_text_input_parallel net1.txt res2.txt 1
	cmd = 
			paste(
					MP.GAEC.KLj.EXECUTABLE.PATH,
					paste0("-i", " ", input.file),
					paste0("-o", " ", file.path(out.folder, paste0(ALGO.RESULT.FILE.PREFIX,0,".txt"))),
					"--multicutRounding gaec_kernighan_lin",
					paste0("--timeout", " ", tilim),
					sep=" "
			)
	
	print(cmd)
	
	return(cmd)
}




#############################################################################################
#
#############################################################################################
prepare.algo.output.filename = function(part.folder, algo.name, g.name){
    
    if(algo.name == COR.CLU.ExCC)
    {
#        ExCC.output.file <- file.path(part.folder, "ExCC-result.txt")
#        id=0
#        file.rename(from=ExCC.output.file, to=file.path(part.folder, paste0(ALGO.RESULT.FILE.PREFIX,id,".txt")))
    }
	else if(startsWith(algo.name,ENUMCC)){
		# do nothing
	}
	else if(startsWith(algo.name,RNSCC)){
		# do nothing
	}
	else if(startsWith(algo.name,COR.CLU.ILS) || startsWith(algo.name,COR.CLU.GRASP))
	{
		# identify the latest folder
		temp.folder <- file.path(part.folder, g.name)
		details <- file.info(list.dirs(temp.folder))
		if(nrow(details)==0)
			stop("load.mestrado.partition: No subfolder was found in folder \"", temp.folder ,"\", cannot load the partition.")
		details <- details[with(details, order(as.POSIXct(mtime))), ]
		temp.folders <- rownames(details)
		folder.name <- file.path(temp.folders[length(temp.folders)])
		
		# set up the file name
		tmp <- strsplit(x=algo.name, split="_", fixed=TRUE)[[1]]
		rcc.flag <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][2]
		f.name=NA
		if(rcc.flag=="RCC")
			f.name <- "rcc-result.txt"
		else
			f.name <- "cc-result.txt"
		file.name <- file.path(folder.name, f.name)
		
		# retain only the partiton info file by changing its name
		id=0
		file.copy(from=file.name, to=file.path(part.folder, paste0(ALGO.RESULT.FILE.PREFIX,id,".txt")))
		
		# possibly remove the original algorithm files
		unlink(x=temp.folder, recursive=TRUE)
	}
	else if(startsWith(algo.name,COR.CLU.NIFTY)){
		output.file <- file.path(part.folder, "membership.txt")
		id=0
		file.rename(from=output.file, to=file.path(part.folder, paste0(ALGO.RESULT.FILE.PREFIX,id,".txt")))
	}
	else if(startsWith(algo.name,COR.CLU.MLMSB)){
		output.file <- file.path(part.folder, "membership.txt")
		id=0
		file.rename(from=output.file, to=file.path(part.folder, paste0(ALGO.RESULT.FILE.PREFIX,id,".txt")))
	}
	else if(startsWith(algo.name,COR.CLU.SA.CC)){
		# do nothing
	}
	else if(startsWith(algo.name,COR.CLU.TS.CC)){
		# do nothing
	}
	else if(startsWith(algo.name,COR.CLU.Brusco.VNS.CC)){
		# do nothing
	}
	else if(startsWith(algo.name,CODE.COR.CLU.GAEC.KLj.CC)){
		# do nothing
	}
	else if(startsWith(algo.name,CODE.COR.CLU.ICP.GAEC.KLj.CC)){
		# do nothing
	}
	else if(startsWith(algo.name,CODE.COR.CLU.MP.GAEC.KLj.CC)){
		# do nothing
	}
    else {
        # TODO
        print("!!!!!!!!!! in TODO")
        # do nothing
    }

    
}




#############################################################################################
# Returns the full name based on the normalized (short) name. Note that for parameterized 
# algorithms, this will just return a clean version of the short name, since it contains 
# the parameter values.
#
# algo.names: short names of the considered algorithms.
#
# returns: the corresponding full names, to be used in plots for instance.
#############################################################################################
get.algo.names <- function(algo.names)
{	result <- c()
	
	for(algo.name in algo.names)
	{
		# parameters
		result <- c(result, gsub(pattern="_", replacement=" ", x=algo.name, fixed=TRUE))
	}
	
	return(result)
}



#############################################################################################
# Returns the inline command for the specified algorithm. The "..." parameters are fetched
# to the algorithm-specific function.
#
# algo.name: short code associated to the algorithm, including its parameter values.
#
# returns: the command allowing to invoke the program externally.
#############################################################################################
get.algo.commands <- function(algo.names, ...)
{	result <- c()

    # substring(x, 1, nchar(prefix)) == prefix
    
    for(algo.name in algo.names)
    {	
        if(startsWith(algo.name,COR.CLU.ExCC))
            result <- c(result, get.ExCC.command(algo.name, ...))
		else if(startsWith(algo.name,ENUMCC))
			result <- c(result, get.EnumCC.command(algo.name, ...))
		else if(startsWith(algo.name,RNSCC))
			result <- c(result, get.RNSCC.command(algo.name, sol.lim=-1, ...))
		else if(startsWith(algo.name,CoNSCC))
			result <- c(result, get.RNSCC.command(algo.name, sol.lim=2, ...))
		if(startsWith(algo.name,COR.CLU.ILS))
			result <- c(result, get.ils.command(algo.name, ...))
		else if(startsWith(algo.name,COR.CLU.GRASP))
			result <- c(result, get.grasp.command(algo.name, ...))
		else if(startsWith(algo.name,COR.CLU.NIFTY))
			result <- c(result, get.NIFTY.command(algo.name, ...))
		else if(startsWith(algo.name,COR.CLU.MLMSB))
			result <- c(result, get.MLMSB.command(algo.name, ...))
		else if(startsWith(algo.name,COR.CLU.SA.CC))
			result <- c(result, get.SA.CC.command(algo.name, ...))
		else if(startsWith(algo.name,COR.CLU.TS.CC))
			result <- c(result, get.TS.CC.command(algo.name, ...))
		else if(startsWith(algo.name,COR.CLU.Brusco.VNS.CC))
			result <- c(result, get.Brusco.VNS.CC.command(algo.name, ...))
		else if(startsWith(algo.name,CODE.COR.CLU.GAEC.KLj.CC))
			result <- c(result, get.GAEC.KLj.command(algo.name, ...))
		else if(startsWith(algo.name,CODE.COR.CLU.ICP.GAEC.KLj.CC))
			result <- c(result, get.ICP.GAEC.KLj.command(algo.name, ...))
		else if(startsWith(algo.name,CODE.COR.CLU.MP.GAEC.KLj.CC))
			result <- c(result, get.MP.GAEC.KLj.command(algo.name, ...))
    }
    
    return(result)
}
