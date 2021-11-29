# TODO: Add comment
# 
# Author: arinik9
###############################################################################




##################################################################
#
#
##################################################################
collect.heuristic.statistics = function(graph.sizes, d, l0.values, heuristic.algos, force)
{
	imb.results = c()
	row.names.desc = c()
	
	prop.mispl.int = 0
	prop.mispl.ext = 0
	
	tlog("starts collecting")
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
					network.path = file.path(in.folder, graph.name)
					g = read.graph.ils(network.path)
						
						
					e.part.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext,
							prop.neg, COR.CLU.ExCC, rep.no=NA)
					membrshp.filepath = file.path(e.part.folder,"membership0.txt")
					membrshp = read.table(membrshp.filepath)$V1
					opt.imb = compute.imbalance.from.membership(g, membrshp, output.type = "count")
						
					
					curr.imb.results = c()
					for(h.algo.name in heuristic.algos){
						
						h.part.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, h.algo.name, rep.no=NA)
						membrshps = load.membership.files(h.part.folder)
						best.imb = Inf
						for(membrshp in membrshps){
							curr.imb = compute.imbalance.from.membership(g, membrshp, output.type = "count")
							if(curr.imb < best.imb)
								best.imb = curr.imb
						}
						curr.imb.results = c(curr.imb.results, best.imb)
					}
					
					imb.results = rbind(imb.results, c(opt.imb,curr.imb.results))
					row.names.desc = c(row.names.desc, in.folder)
				}
			}
		}
	}
	
	colnames(imb.results) = c(COR.CLU.ExCC,heuristic.algos)
	rownames(imb.results) = row.names.desc
	
	if(length(imb.results)>0){
		if(!dir.exists(OUTPUT.CSV.FOLDER))
			dir.create(OUTPUT.CSV.FOLDER, recursive=FALSE, showWarnings=FALSE)
		csv.path = file.path(OUTPUT.CSV.FOLDER,paste0("heuristic-imb-results_d=",d,".csv"))
		write.csv(file=csv.path, x=imb.results)
	}
	
}