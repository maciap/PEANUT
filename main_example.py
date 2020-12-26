import time 
from io_utils import read_file_bif_format
from LRDP_SOSP import LRDP_SOSP
from query_processing import Query, QueryProcessor, QueryProcessor_naive
import argparse


if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description='Materialization and Re-use in Junction Tree Inference')
    parser.add_argument('dataset', help='dataset')
    parser.add_argument('K', help='LRDP SOSP space budget' , type=int)
    parser.add_argument('epsilon', help='approximation parameter' , type=float)
     
    # read the arguments
    args = parser.parse_args()
        
    # create junction tree 
    filepath = "datasets/" + args.dataset + ".bif" 
    join_tree = read_file_bif_format(filepath)

    start_time = time.time()
  
    '''read queries'''
  
    q = Query(join_tree)
    all_queries_training, all_queries_test = q.generate(n = int(3000), probs = "power", max_size = 6) 
    
    f_write_total_skept = open("outputs/total_skept_" + args.dataset + "_" + str(args.K) + "_" + str(args.epsilon) + ".txt","w") 
  
    f_write_total_skept_absolute = open("outputs/total_skept_absolute_" + args.dataset + "_" + str(args.K) + "_" + str(args.epsilon) + ".txt","w") 
  
    ''' offline computations''' 
                
    budget = args.K 
    
    q.set_probabilties(queries = all_queries_training) 
    lrdp_sosp = LRDP_SOSP(join_tree, space_budget=budget, dataset = args.dataset, epsilon = args.epsilon)  
    lrdp_sosp.lrdp_single_sp() 
    lrdp_sosp.reconstruct_solution_lrdp() 
    optimal_benefit = lrdp_sosp.sosp_multiple_sp() 
    lrdp_sosp.reconstruct_solution_sosp(0,budget)     
    lrdp_sosp.compute_all_shortcut_potentials() 
    
            
    ''''''''''''''''''''''''''
    ''' Query Processing '''
    ''''''''''''''''''''''''''
    # we pass @lrdpsosp for simplicity of initialization 
    qp_naive = QueryProcessor_naive(lrdp_sosp)  
    qp_naive.create_map() 

    qp = QueryProcessor(lrdp_sosp)
    qp.create_map() 


    ''' No Materialization '''
    
    to_skip = set() 
    
    start_time = time.time()
    
    q_i=0 
    
    total_skept_queries = list() 
    
    #for query in all_queries: 
    while q_i < len(all_queries_test): 
            
        query = all_queries_test[q_i] 
       
        queryvars = set([query[i].id for i in range(len(query))]) 
   
        steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot = qp_naive.extract_steiner_tree(queryvars , 0)
                    
        total_skept_this_query = qp_naive.compute_total_skept_jt(queryvars, steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot)
        
        if len(steiner_nodes) == 1 or total_skept_this_query==0: 
            to_skip.add(q_i) 
        
        diameter = qp_naive.compute_diameter(steiner_pivot, steiner_children, steiner_nodes) 
#                            
        total_skept_queries.append(total_skept_this_query)
    
        f_write_total_skept_absolute.write( str(total_skept_this_query) + " " + str(diameter) + " \n" )
      
        q_i+=1
           
                
    '''LRDP SOSP'''
    
    total_time = 0 
    
    start_time = time.time()
    
    q_i=0 
    
    while q_i < len(all_queries_test):
        
        if q_i not in to_skip: 
            
            query = all_queries_test[q_i] 
                                 
            tot_this_query = total_skept_queries[q_i] 
                                     
            queryvars = set([query[i].id for i in range(len(query))]) 
        
            #lrdp-sosp 
            steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot = qp.extract_steiner_tree(queryvars , 0)
   
            steiner_nodes , steiner_children, steiner_parents ,steiner_pivot ,  map_id_sp , total_skept = qp.find_sp(steiner_edges,  queryvars, steiner_nodes, steiner_pivot, steiner_parents, steiner_children)
                                    
            total_done_sp = qp.compute_total_skept_jt( queryvars, steiner_nodes , steiner_parents , steiner_children , steiner_pivot, map_id_sp )
            
            total_skept = tot_this_query - total_done_sp
            
            to_write = total_skept / tot_this_query
    
            f_write_total_skept.write( str(to_write)  + " \n"  )
                                
            total_time += time.time() - start_time  
            start_time = time.time() 
     
        q_i+=1
        
   
    f_write_total_skept.close()
    f_write_total_skept_absolute.close() 
   