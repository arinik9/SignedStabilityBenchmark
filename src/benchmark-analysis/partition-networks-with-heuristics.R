# TODO: Add comment
# 
# Author: arinik9
###############################################################################


##################################################################
#
#
##################################################################
partition.networks.with.heuristics = function(graph.sizes, d, l0.values, heuristic.algos, 
		heuristic.repetitions, force)
{
	prop.mispl.int = 0
	prop.mispl.ext = 0
	
	tlog("starts fast jump model stats")
	for(n in graph.sizes){
		for(l0 in l0.values){
			tlog(8, "performing benchmark => n: ", n, ", l0: ",l0)
			
			in.folder = get.benchmark.input.network.folder.path(n, l0, d, NA, NA, NA)
			dirs = list.dirs(path = in.folder, full.names = FALSE, recursive = FALSE)
			print(dirs)
			
			for(dir in dirs){
				# dir = 'propMisplInt=0.0000_propMisplExt=0.0000'
				parts = unlist(strsplit(dir, "_"))
				prop.mispl.int = as.numeric(unlist(strsplit(parts[1],"="))[2])
				prop.mispl.ext = as.numeric(unlist(strsplit(parts[2],"="))[2])
				in.folder = get.benchmark.input.network.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, NA)
				dirs2 = list.dirs(path = in.folder, full.names = FALSE, recursive = FALSE)
				
				for(dir2 in dirs2){
					prop.neg = as.numeric(unlist(strsplit(dir2,"="))[2])
					#print(force)
					
					graph.desc.name = "signed-unweighted" # SIGNED.UNWEIGHTED.FILE
					graph.name = paste0(graph.desc.name, ".G")
					in.folder = get.benchmark.input.network.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg)
					print(in.folder)
					
					for(algo.name in heuristic.algos){
					
						for(rep.no in heuristic.repetitions){
							part.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, rep.no=rep.no)
							if(!dir.exists(part.folder))
								dir.create(part.folder, showWarnings = FALSE, recursive = TRUE)
							
							graph.name = paste0(graph.desc.name, ".G")
							temp.file = file.path(in.folder, "temp.txt")
							
							# ------------------------------------
							if(startsWith(algo.name,CODE.COR.CLU.GAEC.KLj.CC) || 
									startsWith(algo.name,CODE.COR.CLU.ICP.GAEC.KLj.CC) ||
									startsWith(algo.name,CODE.COR.CLU.MP.GAEC.KLj.CC) ){
								input.file = file.path(in.folder, graph.name)
								
								con <- file(input.file, "r")
								lines <- readLines(con)
								close(con)
								lines[1] = "MULTICUT"
								for(i in 2:length(lines))
									lines[i] = gsub("\t"," ",lines[i])
								write.table(x=paste(lines, collapse="\n"), file=temp.file, row.names=F, col.names=F, quote=FALSE)
								
								graph.name = "temp.txt"
							}
							# ------------------------------------
							
							#apply.partitioning.algorithm(part.folder, in.folder, algo.name, graph.name, plot.formats)
							process.partitioning.algorithm(part.folder, in.folder, algo.name, graph.name, 
									FALSE, NA, NA, force)
			
							# ------------------------------------
							if(startsWith(algo.name,CODE.COR.CLU.GAEC.KLj.CC) || 
									startsWith(algo.name,CODE.COR.CLU.ICP.GAEC.KLj.CC))
							{
								if(file.exists(temp.file))
									unlink(temp.file)
							}
							# ------------------------------------

						}
						
						part.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, rep.no=NA)
						process.repetitions(part.folder, heuristic.repetitions, force)
					}
				}
			}
		}
	}
}