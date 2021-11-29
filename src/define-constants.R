

#############################################################################################
# Folders
#############################################################################################
# main folder
MAIN.FOLDER <- "."
#MAIN.FOLDER <- "/home/vlabatut/eclipse/remworkspaces/Networks/NetVotes"
# external libraries folder
LIB.FOLDER <- file.path(MAIN.FOLDER,"lib")
# general input folder
IN.FOLDER <- file.path(MAIN.FOLDER,"in")
#IN.FOLDER <- file.path(MAIN.FOLDER,"in-k=5-n20_50")
#IN.FOLDER <- file.path(MAIN.FOLDER,"in-k=2-n20_50")
	# folder containing signed random networks
    NETWORKS.FOLDER <- file.path(IN.FOLDER,"random-networks")
 	BENCHMARK.NETWORKS.FOLDER <- file.path(IN.FOLDER,"random-networks")
	REAL.NETWORKS.FOLDER <- file.path(IN.FOLDER,"real-networks")
# general ouput folder
OUT.FOLDER <- file.path(MAIN.FOLDER,"out")   # ====> TODO "out-random-networks", "out-real-networks"

	# folder containing the partitions corresponding to the document-wise networks
#	PARTITIONS.FOLDER <- file.path(OUT.FOLDER,"partitions")
	PARTITIONS.FOLDER <- file.path(OUT.FOLDER,"partitions")
    BENCHMARK.ANALYSIS.FOLDER = file.path(OUT.FOLDER,"benchmark-analysis")
    BENCHMARK.PARTITIONS.FOLDER = file.path(BENCHMARK.ANALYSIS.FOLDER,"partitions")
	OUTPUT.CSV.FOLDER = file.path(BENCHMARK.ANALYSIS.FOLDER,"csv")

GRAPH.FILENAME = "signed-unweighted"
SIGNED.UNWEIGHTED.FILE = "signed-unweighted" # undirected as well
SIGNED.WEIGHTED.FILE = "signed-weighted" # undirected as well

PLOT.AS.PDF = "PDF"
PLOT.AS.JPEG = "JPEG"
PLOT.AS.PNG = "PNG"
JUST.PLOT = "NA"


ALGO.RESULT.FILE.PREFIX = "sol"
MBRSHP.FILE.PREFIX = "membership"
# script filename
RECORD.MEM.INFO.SCRIPTNAME = "record-mem.sh"

USED.MEM.VAL.FILENAME = "used-ram-memory.txt" #warning: be consistent with 'record-mem.sh' for the filename

EXEC.TIME.FILENAME = "exec-time.txt"


#############################################################################################
# CSV Column names
#############################################################################################

# colnames used for each network
EXEC.TIME.COL.NAME = "exec time (s)"
NB.SOL.COL.NAME = "nb solution"
NB.OPT.SOL.COL.NAME = "nb optimal solution"
NB.CLU.COL.NAME = "nb cluster"
CLUSTER.COL.NAME.PREFIX = "clu"
IMB.COUNT.COL.NAME = "imbalance count"
IMB.PERC.COL.NAME = "imbalance percentage"
IMB.OPTIMALITY.PERC.COL.NAME = "imbalance optimality (%)"
MEAN.COL.NAME = "mean"
SD.COL.NAME = "sd"
ID.COL.NAME = "id"
FREQ.COL.NAME = "frequency"
ASSOCIATION.COL.NAME = "association"


#############################################################################################
# Some parameters used for evaluating the partition results
#############################################################################################
EVAL.EXEC.TIME.FILENAME = "exec-time"


# post processing of heuristic solutions
HEURISTIC.SOLS.SUMMARY.FILENAME = "heuristic-solutions-summary"
EVAL.CLOSENESS.ASS.INFO.FILE.PREFIX = "closeness-association-info"


