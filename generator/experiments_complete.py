# imports 
#import argparse
import time 
import pickle
from io_utils import read_file_bif_format
from LRDP_SOSP_fast import LRDP_SOSP
from INDSEP_fast import INDSEP
from query_processing_fast import Query, QueryProcessor, QueryProcessor_indsep
import os
import sys 

#from SHORTCUT_POTENTIAL import ShortcutPotential
#from graph.potential import PotentialUtil
def read_pickle(filename): 
    objects = []
    with (open(filename, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break
    return objects




if __name__ == '__main__':
    
    #parser = argparse.ArgumentParser(description='Materialization and Re-use in Junction Tree Inference')
    # arguments
    #parser.add_argument('K', help='INDSEP space budget' , type=int)
    #parser.add_argument('k_ls', help='LRDP SOSP space budget' , type=int)
    #parser.add_argument('d', help='dataset')
    #parser.add_argument('outfile', help='output file path')
    #parser.add_argument('n', help='number of queries', type=int,  default = 1000)
    
    #parser.add_argument('probs'Perché ora è un altro tipo di giocatore, dice il procuratore. Che dev'essere un omettino col suo quoziente intellettivo, o giù di lì.
    # , help='probability distriribution of variables being queried', default = "power")
        
    # read the arguments
    #args = parser.parse_args()
        
    #f_write = open(args.outfile,"w") 
    
    #f_write_total_skept = open("outputs/total_skept_" + args.d  + "_" + str(args.K) + "_" + str(args.k_ls) + ".txt","w") 
    
    #now args.K and args.d contains the arguments 
    
    ''' offline computations''' 
    
    #alldatasets = [ "child", "alarm", "insurance"]
                   #, "insurance", "hailfinder", "hepar2", "andes", 
    #alldatasets = [ "alarm", "insurance", "hailfinder", "andes", "hepar2"]
        #"child", "alarm", "insurance", "hailfinder", "andes", "hepar2",
    alldatasets = ["munin", "barley", "mildew" ] 
                   #"pathfinder", "diabetes" ]
    k_indsep_list = [100, 200, 300, 400, 500] 
    
    #k_lrdp_list = []
    
    
    for this_dataset in alldatasets: 
        
        start_time = time.time()
                 
        filepath = "datasets/" + this_dataset + ".bif"
        
        join_tree = read_file_bif_format(filepath)
        
        print( " tree created --- " + str( (time.time() - start_time) ) + " --- seconds \n"  )  
        
        
        for k_indsep in k_indsep_list: 
            
            
            
            outfile = "outputs/total_skept_" + this_dataset + "_" + str(k_indsep)
            outfile_time = "outputs/total_skept_running_time_" + this_dataset + "_" + str(k_indsep)
            
            f_write_total_skept = open(outfile, "w") 
            f_write = open(outfile_time, "w")
        
       
            print("current dataset " + str(this_dataset))
            
            
            start_time = time.time()
            
            '''generate queries'''
            
            # generate queries using a power low distribution for single variable probabilities 
            q = Query(join_tree)
            probabilites, all_queries = q.generate(n = int(40), probs = "power", max_size = 5) 
            q.set_probabilties(probs = probabilites)
            
            # pickle queries 
            #>pickle_path_queries = "/u/50/ciaperm1/unix/Desktop/CODE_FAST_APPROXIMATION_BIN/pickles/queries_" + args.d
            #pickle.dump(all_queries , open(pickle_path_queries,"wb"))
            
            
            print("queries generated -- " + str( (time.time() - start_time) ) + " --- seconds \n"  )  
            
            ''''''''''''''''''''''''''
            '''INDSEP'''
            ''''''''''''''''''''''''''
            
            # create INDSEP 
            
            start_time = time.time()
            
            indsep_budget = k_indsep
            indsep = INDSEP(join_tree , space_budget=indsep_budget)  
            allconnected_components = indsep.construct_INDSEP()
            dataset = this_dataset
           # pickle_path_indsep = "/u/50/ciaperm1/unix/Desktop/CODE_FAST_APPROXIMATION_BIN/pickles/indsep_" + dataset
           # pickle.dump(indsep , open(pickle_path_indsep,"wb"))
            
            print("INDSEP done. Preprocessing time --- " + str( (time.time() - start_time) ) + " --- seconds \n"    )  
            
            
            ''''''''''''''''''''''''''
            '''LRDP SOSP'''
            ''''''''''''''''''''''''''
            # budget for LRDP SOSP 
            
            budget = sum( indsep.weight_sp.values() )   
            
          #  budget =  int(k_ls) 
            
            #print("budget for lrdp sosp " + str(budget))
            
            # LRDP SOSP starting 
            
            start_time = time.time()
            
            # left right dynamic programming 
            lrdp_sosp = LRDP_SOSP(join_tree, space_budget=budget, dataset = this_dataset, epsilon = 1.5)  
            lrdp_sosp.lrdp_single_sp() 
            
            print( "SINGLE LRDPs done --- " + str( (time.time() - start_time) ) + " --- seconds \n" )  
            
            #pickle_path = "/u/50/ciaperm1/unix/Desktop/Code/pickles/lrdp_sosp_" + z dataset
            #pickle.dump(lrdp_sosp , open(pickle_path,"wb"))
            
            start_time = time.time()
            
            lrdp_sosp.reconstruct_solution_lrdp() 
            
            print( "Reconstructed LRDP --- " + str( (time.time() - start_time) ) + " --- seconds \n" )  
            
            
            start_time = time.time()
            
            optimal_benefit = lrdp_sosp.sosp_multiple_sp() 
            
            print( "SOSP executed --- " + str( (time.time() - start_time) ) + " --- seconds \n" )  
            
            
            start_time = time.time()
            
            lrdp_sosp.reconstruct_solution_sosp(0,budget) # assuming root is clique 0
                
            print( "Reconstructed SOSP --- " + str( (time.time() - start_time) ) + " --- seconds \n" )  
            
            
            start_time = time.time()    
            
            lrdp_sosp.compute_all_shortcut_potentials() 
            # lrdp_sosp =  lrdp_sosp.compute_lrdpsosp(budget)
            
            print( "LRDP SOSP done. Preprocessing time --- " + str( (time.time() - start_time) ) + " --- seconds \n"    )  
            
            ''''''''''''''''''''''''''
            ''' Query Processing '''
            ''''''''''''''''''''''''''
            
            qp = QueryProcessor(lrdp_sosp)
            qp.create_map() 
            qp_indsep = QueryProcessor_indsep(indsep.jointree, indsep.all_cliques, indsep) 
            qp_indsep.create_map() 
            
            '''INDSEP'''
            
            total_time_indsep = 0 
            start_time = time.time()
            q_i=0 
            
            #for query in all_queries: 
            while q_i < len(all_queries): 
                
                
                try: 
                    query = all_queries[q_i]    
                    
                    queryvars = set([query[i].id for i in range(len(query))]) 
                    
                    # indsep 
                    # convert query variables 
                    queryvars_indsep = set([indsep.var_mapping_inv[x] for x in queryvars]) 
                    qp_indsep.total_skept = 0
                    output = qp_indsep.query( indsep.root , queryvars_indsep )
                    
                    f_write_total_skept.write( str(qp_indsep.total_skept) + "\n" )
                    
                    f_write.write(  str( (time.time() - start_time) ) + " \n"    )  
                    total_time_indsep += time.time() - start_time  
                    start_time = time.time() 
                
                except: 
                    continue 
                
                q_i+=1
                
                #if q_i % 100 == 0: 
                #print("INDSEP query number" + " " + str(q_i)) 
                    
                
            print("INDSEP total time " + str(total_time_indsep))    
            #f_write.write(  str( (time.time() - start_time) ) + "  seconds \n"    )  
            
            
            f_write.write(  "LRDPSOSP " + str(sum( indsep.weight_sp.values() )) + " \n")  
            f_write_total_skept.write( "LRDPSOSP " + str(sum( indsep.weight_sp.values() )) + " \n"  )
            
            
            '''LRDP SOSP'''
            total_time = 0 
            start_time = time.time()
            q_i=0 
            
            
            #for query in all_queries: 
            while q_i < len(all_queries):
                
                query = all_queries[q_i] 
                
                #print(q_i) 
                
                try: 
            
                    queryvars = set([query[i].id for i in range(len(query))]) 
                            
                    #lrdp-sosp 
                    steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot = qp.extract_steiner_tree(queryvars , 0)
                    steiner_nodes , steiner_children, steiner_parents ,steiner_pivot ,  map_id_sp , total_skept = qp.find_sp(steiner_edges,  queryvars, steiner_nodes, steiner_pivot, steiner_parents, steiner_children)
                    f_write_total_skept.write(str(total_skept) + "\n")
                    
                    # output = qp.query( steiner_nodes, steiner_parents, steiner_children, steiner_pivot,  queryvars, map_id_sp)
                            
                    f_write.write(  str( (time.time() - start_time) ) + " \n"  )  
                    total_time += time.time() - start_time  
                    start_time = time.time() 
            
                except: 
                    continue 
                
                    q_i+=1
                
                #if q_i % 100 == 0: 
                #    print("LRDP SOSP query number" + " " + str(q_i)) 
                    
            
            print("LRDP SOSP total time " + str(total_time))
            
            f_write_total_skept.close()
            f_write.close()