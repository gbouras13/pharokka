#!/usr/bin/env python3
import sys
from modules import input_commands




if __name__ == "__main__":
    print("test")
    args = input_commands.get_input()
    tmp = args.infile
    print(tmp)
    sys.exit("phrokka has finished")  # pragma: no cover




# #--------------------------------------------------------------------------------------------------#
# #                               FILE INPUT                                                         #
# #--------------------------------------------------------------------------------------------------#

# my_contigs = file_handling.read_fasta(args.infile);
# if not my_contigs:
# 	sys.stdout.write("Error: no sequences found in infile\n")
# 	sys.exit()

# #--------------------------------------------------------------------------------------------------#
# #                               MAIN ROUTINE                                                       #
# #--------------------------------------------------------------------------------------------------#
# for id, seq in my_contigs.items():

# 	#-------------------------------Find the ORFs----------------------------------------------#
# 	my_orfs = functions.get_orfs(seq)



# 	#-------------------------------Create the Graph-------------------------------------------#
# 	my_graph = functions.get_graph(my_orfs)



# 	#-------------------------------Run Bellman-Ford-------------------------------------------#
# 	source = "Node('source','source',0,0)"
# 	target = "Node('target','target',0," + str(len(seq)+1) + ")"
# 	# Write edges to the fastpath program, and multiply the weight to not lose decimal places
# 	fz.empty_graph()
# 	for e in my_graph.iteredges():
# 		if args.dump: print(e)
# 		ret = fz.add_edge(str(e))

# 	if args.dump: sys.exit()

# 	shortest_path = fz.get_path(source=source, target=target)


	
# 	#-------------------------------Write Output ----------------------------------------------#
# 	file_handling.write_output(id, args, shortest_path, my_graph, my_orfs)

# #--------------------------------------------------------------------------------------------------#
# #                               END                                                                #
# #--------------------------------------------------------------------------------------------------#