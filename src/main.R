
# =============================================================================
# VARIABLES
#   ==> REmark 1: Do not forget to set CPLEX.BIN.PATH correclty in src/define-algos.R
#   ==> Remark 2: PROP.NEGS should be set to 'NA' when DENSITY=1.0
#   ==> Remark 3: CORE.PART.THRESHOLD should be set to 1. Because, when it is less than 1,
#                  there might be multiple way of building core part.
#   ==> Remark 4: It is the responsability of the user who will ensure if the RAM requirement
#                  of his/her system is ok for large graphs, because Cplex may require 
#                  lots of RAM for graphs whose size is larger than 28.
# =============================================================================


## libraries for parallel processing
#library(foreach)
#library(doParallel)

source("src/define-imports.R")


#######################################################################
# STEP 0: Set parameter values
#######################################################################=

BENCHMARK.GRAPH.SIZES = c(20,30) #c(30,40,50)
BENCHMARK.DENSITY = 0.25
BENCHMARK.L0.VALUES = c(3,4)
BENCHMARK.PROP.NEGS = c(0.3,0.5) #c(0.3,0.5,0.7)
BENCHMARK.FORCE = FALSE
BENCHMARK.PLOT.FORMAT = JUST.PLOT


plot.format <- c( # ==========> it is not taken into account everywhere !! TODO
		#PLOT.AS.PDF
		#PLOT.AS.JPEG
		#PLOT.AS.PNG
		JUST.PLOT
)

FORCE = FALSE

keep.algo.log.files = TRUE




#######################################################################
# STEP 1: Generate perfectly balanced signed networks and place them into 'in' folder
#######################################################################

generate.all.perfectly.balanced.networks(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, BENCHMARK.PROP.NEGS)




#######################################################################
# STEP 2: Perturb in many ways the perfectly balanced signed networks
#			and obtain new signed networks
#######################################################################

create.all.extreme.imbalanced.networks(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, BENCHMARK.PROP.NEGS, BENCHMARK.FORCE, BENCHMARK.PLOT.FORMAT)
perform.all.benchmark(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, BENCHMARK.FORCE)
collect.all.benchmark.results(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, BENCHMARK.FORCE)

plot.benchmark.results(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, 3, BENCHMARK.FORCE)
plot.benchmark.results(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, 4, BENCHMARK.FORCE)



############################################################
# Heuristics
#############################################################

#UNUSED.COR.CLU.HEURISTIC.ALGOS = c(
#    get.ils.code(l=1, alpha=0.4, gain=0, perturbation=3, time.limit=3600, iter.nbr=10, rcc=FALSE, vns=TRUE), # VNS metaheursitic
##    get.grasp.code(rcc=FALSE, l=1, k=NA, alpha=0.8, gain=0, time.limit=3600, iter.nbr=-1, oneOptNeig=0), # VOTE-BOEM
##    ###get.NIFTY.code(method="greedy-additive"),
##    ###get.NIFTY.code(method="cgc-qpbo"),
##    ###get.NIFTY.code(method="greedy-additive + cgc-qpbo"),
#    #######get.ZONOCC.code(rank=3)
##    ######## get.SPONGE.sym.CC.code()
#)		
COR.CLU.HEURISTIC.ALGOS = c(
		get.ils.code(l=1, alpha=0.4, gain=0, perturbation=3, time.limit=30, iter.nbr=10, rcc=FALSE),
		get.grasp.code(rcc=FALSE, l=1, k=NA, alpha=0.8, gain=0, time.limit=30, iter.nbr=400, oneOptNeig=1),
		get.MLMSB.code(trial=1,time.limit=30),
		get.SA.CC.code(time.limit=30),
		get.TS.CC.code(time.limit=30),
		get.Brusco.VNS.CC.code(time.limit=30),
		get.ICP.GAEC.KLj.CC.code(time.limit=30), # error with large graphs
		get.GAEC.KLj.CC.code(time.limit=30),
		get.MP.GAEC.KLj.CC.code(time.limit=30)
)
HEURISTIC.REPETITIONS = seq(1, 3, by=1)


partition.networks.with.heuristics(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, 
		COR.CLU.HEURISTIC.ALGOS, HEURISTIC.REPETITIONS, BENCHMARK.FORCE)

collect.heuristic.statistics(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, 
		COR.CLU.HEURISTIC.ALGOS, BENCHMARK.FORCE)

	
	