#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#import operator 
#import random 
#from collections import defaultdict 
#from INDSEP import compute_single_shortcut_potential 
# from graph.potential import multiply_and_return
# import sys 
import pickle 
import copy 
import itertools
from collections import  defaultdict
import operator 
# import graph.potential
from graph.potential import PotentialUtil , Potential, PotentialEntry
# import numpy as np 
# import copy 
from SHORTCUT_POTENTIAL import ShortcutPotential 
import time 
from anycache import anycache



class LRDP_SOSP():
    """
    INDSEP class (Kanagal) 
    ISSUE: NOT SURE ABOUT HOW TO COMPUTE SHORTCUT POTENTIALS 
    """
    
    def __init__(self, jointree, space_budget, dataset, epsilon = None):
        self.jointree = jointree
        self.K = space_budget  # target disk block size
        self.adj_lst , self.map_nodes  , self.map_nodes_inv = jointree.convert_to_adjancency_list() 
        
        # add dummy clique as parent of the root 
        self.adj_lst[0].append(-1) 
        self.adj_lst[-1].append(0) 
        
        if epsilon == None: 
            self.epsilon = max(1, int(self.K/50))
        else: 
            self.epsilon = epsilon
            
        self.all_cliques = jointree.get_cliques()
        self.initialize()  # dyanmic programming table  
        self.probabilites = dict() # implement function to compute it 
        self.optimal_packing = dict()
        # visited for dfs 
        self.visited = set()  
        self.dfs_labels_map = dict() 
        self.dfs_labels_map_inv = dict() 
        self.dfs_label = 0  
        self.parents = dict() 
        self.children =  defaultdict(set) 
        self.leaves = set() 
        self.visited_leaves = set()
        self.probabilites = dict() 
        self.allpaths = [] 
        self.heights = dict() 
        self.heights[-1] = 0 
        self.dfs(-1) 
        
        self.dataset = dataset 
        self.dict_var_clique = self.create_var_clique_mapping() 
        
        
    def create_var_clique_mapping(self): 
        ''' this function create a mapping between each variable and the 
        clique which contains it closest to the root for efficient query processing
        @returns 
        map_clique_var : dict 
        '''
        
        # cliques sorted in order of increasing heights 
        map_clique_var = dict()
        sorted_cliques = [x[0] for x in sorted(self.heights.items(), key=operator.itemgetter(1))][1:]
        bn_nodes = self.jointree.get_bbn_nodes() 
        for i in range(len(bn_nodes)): 
            this_var = bn_nodes[i]
            for clique in sorted_cliques: 
                if this_var.id in self.all_cliques[clique].node_ids: 
                    map_clique_var[this_var.id] = clique 
                    # go to next bbn node 
                    break

        return map_clique_var
    
    
    
       
    def initialize(self): 
        '''
        intitialize 3d arrays to be used later 

        Returns
        -------
        None.

        '''
        # 2d array 
        
        # create bins 
        
        self.bins = [] 
        
        self.index_bins = dict() 
        
        i = 0 
        
        cnt = 0 
        
        self.bins.append(0) 
        self.index_bins[i] = cnt 
        
        cnt+=1 
        
        while i < self.K + 1:
            
            i += self.epsilon**cnt 
            
            if i < self.K:
                
                i = round(i)
                
                self.bins.append( i ) 
                self.index_bins[i] = cnt
                # increase power 
                cnt+=1 
            
        self.bins.append(self.K) 
        self.index_bins[self.K] = cnt 
        
        print("# bins " + str(len(self.bins)))
        
        
        rows, cols = (len(self.all_cliques), len(self.bins))  
        self.P_final = [[float("-inf") for i in range(cols)] for j in range(rows)] 
        self.I_final = [[0 for i in range(cols)] for j in range(rows)] 
        self.H = [[0 for i in range(cols)] for j in range(rows)] 
        self.H_2 = [[0 for i in range(cols)] for j in range(rows)] 
        
        # 3d array 
        one, two, three = ( len(self.bins) , len(self.all_cliques),   len(self.all_cliques)  )
        self.P =  [[[float("-inf") for i in range(one)] for j in range(two)]   for h in range(three)]
        self.I =  [[[0 for i in range(one)] for j in range(two)]   for h in range(three)]
        self.parent_benefit = [[[float("-inf") for i in range(one)] for j in range(two)]   for h in range(three)]
        self.children_benefit = [[[float("-inf") for i in range(one)] for j in range(two)]   for h in range(three)]
        self.children_scope = [[[set() for i in range(one)] for j in range(two)]   for h in range(three)]
        self.I_children =  [[[0 for i in range(one)] for j in range(two)]   for h in range(three)]
            
        
        
        
    def dfs(self, node):
        '''
        Perform dfs and compute labels, parents, children, leaves and leaf-to-root paths 
        @params : visited = set of visited nodes, initially empty
                   graph: adjancency list of T 
                   node: current node to visit 
        '''
        
        if node not in self.visited:
            self.visited.add(node)

            self.dfs_labels_map[self.dfs_label] = node 
            self.dfs_labels_map_inv[node] = self.dfs_label
            self.dfs_label+=1 
            
            if len(self.adj_lst[node]) == 1 and node != -1:
                
                self.leaves.add(node) 
                # node is a leaf, we add a new path from leaf (node) to root
                this_path = [node] 
                while node!=-1: 
                    this_path.append( self.parents[node] )
                    node = self.parents[node] 
                self.allpaths.append(this_path)     
                    
                    
            for neighbour in self.adj_lst[node]:
                if neighbour not in self.visited: 
                    
                    self.parents[neighbour] = node 
                    self.children[node].add(neighbour) 
                    self.heights[neighbour] = self.heights[node] + 1 
                    
                self.dfs(neighbour)
                
    
    
    
    def compute_all_benefits(self):
        '''
        Compute all benefits of the shortcut potentials in the tree 
       
        @return: benefit_paths which is a dictionary giving the benefit for each node 
        
        '''
   
        benefit_paths = defaultdict(lambda: defaultdict(float))
        
        # one by one we deal with each path 
        for path in self.allpaths: 
                        
            # compute benefits for this single path 
            this_path_benefit = self.compute_single_path(path)
            
            # for each root 
            for root in this_path_benefit.keys(): 
                
             #   print("starting root ")
                # update the overall benefits so far
                for k in this_path_benefit[root].keys(): 

                    # defaultdict initializes to zero if not present already 
                    benefit_paths[root][k] += this_path_benefit[root][k]      
                    
                # print("root " + str(root) + " done")  
                    
                    
        
        return benefit_paths     
            
            
    
    
    def makeCombos(self, arr):
        ''' helper function computing combinations of variables'''
        return (combo for i in range(1, len(arr) + 1) for combo in map(list, itertools.combinations(arr, i)))

    
    
    
    
    def compute_single_path(self, path_nodes): 
        
        ''' compute benefit for a
        single path 
        between c and root
        return all benefit which is a dictionary 
        which maps each 
        benefit (value) to key which is the
        destination node 
        '''
    
        allbenefits = defaultdict(dict)

        for r in range(len(path_nodes)-1): 
            
            # slciing and reversing to get the current nodes in the appropriate order 
            current_path_nodes = path_nodes[:(len(path_nodes)-r)][::-1]
            
            current_root = current_path_nodes[0]
            
            weight = 1
            if r != 0:
                weight *= self.all_cliques[r].get_weight()
                root_sp = self.all_cliques[r].node_ids.intersection(self.all_cliques[self.parents[r]])
            
            if r ==0: 
                root_sp =  self.all_cliques[r].node_ids
        
          
            pr = 0 
        
            all_vars_root = set()
        
            # now accumulating all the query component for the root shortcut potential 
            for node in current_path_nodes[1:]: 
        
                # get variables in current clique 
                vars_current_clique = self.all_cliques[node].nodes
                all_vars_root.update( vars_current_clique )
        
            
            all_combinations = self.makeCombos(all_vars_root)
        
            # now process each combination 
            for comb in all_combinations: 
        
                # for each variable in a single combinations 
                for var in comb:
        
                    # update probability 
                    pr += var.weight * var.probability
        
            # first benefit (for root shortcut potential) 
            allbenefits[current_root][current_path_nodes[1]] =  weight *  pr 
            
            allnodes = []
            
            for i in range(1 , len(current_path_nodes)-1):
                
        
                node = current_path_nodes[i]
                
                allnodes.append(node)
        
                vars_current_clique = self.all_cliques[node].nodes
        
                # find combinations of variables in this clique 
                all_combinations_current_clique = self.makeCombos(vars_current_clique)
        
              
                c_pr = 0 
                for comb in all_combinations_current_clique:
                    for var in comb: 
                        c_pr += var.weight * var.probability
                        # c_pr += var.weight 
        
                # update probability 
                pr -= c_pr 
        
                # update weight 
                weight += self.all_cliques[node].get_weight()  
        
                # CHeck that there is a gain 
                curr = current_path_nodes[i+1] 
                this_int = self.all_cliques[curr].node_ids 
                this_parent_int = self.all_cliques[self.parents[curr]].node_ids 
                this_node_sp = this_int.intersection(this_parent_int)
                union_sp = this_node_sp.union(root_sp)
                
                we = 1 
                for var in union_sp: 
                    we*=var.weight
                    
                other_scope = set() 
                for n in allnodes:
                    other_scope.update(self.all_cliques[n].node_ids) 
                    
                we2 = 1 
                for var in other_scope: 
                    we2*=var.weight    
                    
                allbenefits[current_root][current_path_nodes[i+1]] = 0
                    
                if we2  > we: 
                 
                    allbenefits[current_root][current_path_nodes[i+1]] =   weight *  pr 
                
            
        return allbenefits
    
    
    
    
    
    
    def backward(self, visited, node, r): 
        
        ''' 
        backward move in algorithm for Single optimal shortcut potential  
        @params: 
        visited: set of crurrently visited nodes 
        node: node v to process in the current backward step 
         
        '''

        # we first evaluate combinations of children with the previous 
        # benefit_sum data structure map a key to -inf if key is not present otherwise to a tuple (value, weight1, weight2)
        benefit_sum = defaultdict(lambda: (float("-inf") , 0 , 0))
        
        p_v = self.parents[node] 
        
        # for each weight which is not zero (we can skip zero to avoid problems here 
        # since P will always be -inf for weight zero)
        for c1 in self.bins: 
            
            for c2 in self.bins: 
                
                # note that at the first child this is always -inf 
                comb = self.children_benefit[r][p_v][ self.index_bins[ c1 ] ] +  self.P[r][node][ self.index_bins[ c2 ]  ] 
                
                #weight comb must be given by the combination of (v,p_v) scope and optimal scope at 
                scope_comb = self.children_scope[r][p_v][ self.index_bins[ c1 ] ].union(self.all_cliques[node].node_ids.intersection(self.all_cliques[p_v].node_ids))
                scope_comb_with_root = scope_comb.union(self.all_cliques[r].node_ids.intersection(self.all_cliques[self.parents[r]].node_ids)) 
                
                
                
                weight_comb = 1 
                for var in scope_comb_with_root: 
                    weight_comb = weight_comb * len(self.jointree.get_bbn_node(var).variable.values)  
                
                
                if weight_comb <= self.K:
                    
                    #for c_c in range(weight_comb , self.K+1): 
                    for c_c in self.bins:
                        
                        if c_c >= weight_comb:
                    
                            #right hand side initialized to -inf
                            if comb > benefit_sum[self.index_bins[c_c]][0]: 
    
                                # evaluate sums of benefits 
                                benefit_sum[self.index_bins[c_c]] = (comb, c1, c2) # we add the corresponding tuple 

        
        # now we compare for each weight the benefit associated with each possible case 
        for c in self.bins: 
            
                
            #1) new node is optimal (note that at the first children visited only the last condition matters)
            if self.P[r][node][self.index_bins[c] ] > self.children_benefit[r][p_v][self.index_bins[c]] and self.P[r][node][self.index_bins[c]]  >  benefit_sum[self.index_bins[c]][0] and self.P[r][node][self.index_bins[c]]  >  self.parent_benefit[r][p_v][self.index_bins[c]]: 
                
                # update children benefit becuse new child is alone optimal 
                self.children_benefit[r][p_v][self.index_bins[c]] = self.P[r][node][self.index_bins[c]] 
                
                #update children scope 
                self.children_scope[r][p_v][ self.index_bins[ c ] ] =  self.all_cliques[node].node_ids.intersection(self.all_cliques[p_v].node_ids)
                
                
                # keep also the optimal at the parent level 
                self.P[r][p_v][self.index_bins[c]] = self.P[r][node][self.index_bins[c]] 
                
             
                self.I[r][node][self.index_bins[c]] = 1 
                self.I_children[r][node][self.index_bins[c]] = 1
                
                for ch in self.children[p_v]: 
                    if ch!= node: 
                        self.I[r][ch][self.index_bins[c]] = 0
                        self.I_children[r][ch][self.index_bins[c]] = 0
                        
    
                        
            #2)  parent (untouched by backward procedure) is optimal  
            elif self.parent_benefit[r][p_v][self.index_bins[c]] > self.P[r][node][self.index_bins[c]] and self.parent_benefit[r][p_v][self.index_bins[c]]  >  benefit_sum[self.index_bins[c]][0]  and self.parent_benefit[r][p_v][self.index_bins[c]]  >  self.children_benefit[r][p_v][self.index_bins[c]]:
                
            #   self.I[r][p_v][c] = 1 this is actually not need since in the next backward step we would set this to 1 
                # here we do not need to update the others 
                # children are never set to 1 if parent is optimal since self.parent_benefit stays the same 
            
            
            
                # first we need to set to zero the indicator variable for the children 
                for ch in self.children[p_v]: 
                    self.I[r][ch][self.index_bins[c]] = 0
            
            
              # if parent is optimal we only need to update children benefit because it may later become optimal combined with new children 
                
                # we have three possible cases 
                # 2a) self.P (new node) is alone optimal - update children optimal benefit and indicator for children 
                # 2b) combination with new node is optimal -  update children optimal benefit and indicator for childrens
                # 3c) old children benefit is optimal - do nothing 
                
                if self.P[r][node][self.index_bins[c]] > benefit_sum[self.index_bins[c]][0] and self.P[r][node][self.index_bins[c]] > self.children_benefit[r][p_v][self.index_bins[c]]: 
                    
                    self.children_benefit[r][p_v][self.index_bins[c]] = self.P[r][node][self.index_bins[c]] 
                    
                    #update children scope 
                    self.children_scope[r][p_v][ self.index_bins[ c ] ] =  self.all_cliques[node].node_ids.intersection(self.all_cliques[p_v].node_ids)
                
                    
                    # update indicator variable for children only (not overall)
                    self.I_children[r][node][self.index_bins[c]] = 1
                    for ch in self.children[p_v]: 
                        if ch!= node: 
                            self.I_children[r][ch][self.index_bins[c]] = 0 
                
                elif benefit_sum[self.index_bins[c]][0] > self.P[r][node][self.index_bins[c]] and benefit_sum[self.index_bins[c]][0] > self.children_benefit[r][p_v][self.index_bins[c]]: 
                    v, c1, c2 = benefit_sum[self.index_bins[c]]
                    self.children_benefit[r][p_v][self.index_bins[c]] = v
                    # combination is optimal, in this case we update 
                    #self.children_scope[r][p_v][ self.index_bins[ c ] ].update(self.all_cliques[node].node_ids.intersection(self.all_cliques[p_v].node_ids))
                
                    self.children_scope[r][p_v][ self.index_bins[ c ] ] = self.children_scope[r][p_v][ self.index_bins[ c1 ] ].union(self.all_cliques[node].node_ids.intersection(self.all_cliques[p_v].node_ids))
                    
                    
                    for ch in self.children[p_v]: 
                        if self.I_children[r][ch][self.index_bins[c1]] == 1 or ch == node: 
                            self.I_children[r][ch][self.index_bins[c]] = 1 
                        else: 
                            # make sure that all the others are null
                            self.I_children[r][ch][self.index_bins[c]] = 0
                        
                    
            #3) previous optimal is still optimal - do nothing 
            
            #4) new combination is optimal 
            elif benefit_sum[self.index_bins[c]][0] > self.P[r][node][self.index_bins[c]] and benefit_sum[self.index_bins[c]][0] > self.parent_benefit[r][p_v][self.index_bins[c]] and benefit_sum[self.index_bins[c]][0]  >  self.children_benefit[r][p_v][self.index_bins[c]]:
                
                v, c1, c2 = benefit_sum[self.index_bins[c]] # c1 is the previous children weight contribution , c2 is the new one 
                
                # update children_benefit 
                self.children_benefit[r][p_v][self.index_bins[c]] = v
                #self.children_scope[r][p_v][ self.index_bins[ c ] ].update(self.all_cliques[node].node_ids.intersection(self.all_cliques[p_v].node_ids))
                
                self.children_scope[r][p_v][ self.index_bins[ c ] ] = self.children_scope[r][p_v][ self.index_bins[ c1 ] ].union(self.all_cliques[node].node_ids.intersection(self.all_cliques[p_v].node_ids))
                    
                
                # keep the optimal at the parent level 
                self.P[r][p_v][self.index_bins[c]] = v

                # finally update the indicator variable 
                # we set the otpimal for c plus the new one 
                for ch in self.children[p_v]: 
                    if self.I_children[r][ch][self.index_bins[c1]] == 1 or ch == node: 
                        self.I[r][ch][self.index_bins[c]] = 1 
                        self.I_children[r][ch][self.index_bins[c]] = 1 
                    else: 
                        # make sure that all the others are null
                        self.I[r][ch][self.index_bins[c]] = 0 
                        self.I_children[r][ch][self.index_bins[c]] = 0 
       
                    
         
            
    def reconstruct_solution_lrdp(self): 
        
        '''
        
          in the end self.sep will contain for each root and each weight the set of optimal 
          separators to be used in the 
          procedure
          dict: root / weight to list of edges that are associated with the separators of the optimal shortcut potential 
          for that given weight
          
        '''
                
        self.seps = defaultdict(lambda: defaultdict(set))
        
        for r in range(len(self.all_cliques)): 
            
            #print("I_r")
        
            #print(self.I[r])
                        
            for c in self.bins: 
            
                # optimal at the root (kept for sosp)
                
                self.P_final[r][self.index_bins[c]] = self.P[r][r][self.index_bins[c]] 

                # we need to iterate in a forward fashion 
                to_visit =  copy.deepcopy(self.children[r]) # we initialize the set of nodes to visit with the children of the root
                visited = set() 

                node = r 
                path_nodes = set() 
                
                #for i in range( self.dfs_labels_map_inv[r] + 1 , len(self.all_cliques)): 
                while len(to_visit) > 0: 

                    # get the next node to visit 
                    node = self.dfs_labels_map[ self.dfs_labels_map_inv[node] + 1 ] 
                    
                    
                    if node in to_visit:
                        to_visit = to_visit.difference({node}) 
                        visited.add(node) 
                        
                        if node in self.leaves and self.I[r][node][self.index_bins[c]] == 1:
                            
                            # add optimal separator 
                            self.seps[r][self.index_bins[c]].add(( node , self.parents[node] )) 
                            #to_visit = to_visit.difference({node}) 
                            path_nodes = set() 
    
                        elif self.I[r][node][self.index_bins[c]] == 0 and len(path_nodes)>0 and len(self.children[self.parents[node]].difference(visited)) == 0:
                                                    
                            self.seps[r][self.index_bins[c]].add(( self.parents[node] , self.parents[self.parents[node]] )) 
                            #to_visit = to_visit.difference({node}) 
                            path_nodes = set() 
    
    
                        elif self.I[r][node][self.index_bins[c]] == 1:
                            to_visit.update(self.children[node]) 
                            path_nodes.add(node) 
    
            
        

                    


    def lrdp_single_sp(self): 
        
        '''
        compute optimal single 
        shortcut potential 
        with left-right 
        dynamic programming 
        
        '''
        
        benefits = self.compute_all_benefits() 
        
        print("starting LRDP")
        
        for r in range(len(self.all_cliques)):
                        
            root = self.all_cliques[r]
            
            par_root = self.parents[r]
            
            if par_root == -1: #dummy clique 
                
                root_weight = 1 
                vars_to_multiply_root = set({})
            
            else:    
                
                vars_to_multiply_root = self.all_cliques[ par_root ].node_ids.intersection(root.node_ids)
                
                root_weight = 1 
                for var in vars_to_multiply_root:
                    root_weight = root_weight * len(self.jointree.get_bbn_node(var).variable.values)  
                    
             
            visited = set() 
            
         
            to_visit = copy.deepcopy( self.children[r] ) 
            
            node = r 
                     
            #for i in range( self.dfs_labels_map_inv[r] + 1 , len(self.all_cliques)): 
            while len(to_visit) > 0: 
                                
                # get the next node to visit 
                node = self.dfs_labels_map[ self.dfs_labels_map_inv[node] + 1 ] 
                
              #  print("forward " + str(node))
                                
                if node not in visited: 
                    
                    # update set of visited nodes 
                    visited.add( node ) 
                    
                    # node2 = self.jointree.all_cliques[self.dfs_labels_map[i-1]]
                    par = self.parents[ node  ]
                    
                    vars_to_multiply = self.all_cliques[ par ].node_ids.intersection( self.all_cliques[ node ].node_ids )
                    
                    #if self.jointree.potentials.get(sepset.id):
                    this_sep_set_weight = 1 
                    for var in vars_to_multiply.union(vars_to_multiply_root):
                        this_sep_set_weight = this_sep_set_weight * len(self.jointree.get_bbn_node(var).variable.values)  
                    
                    c_weight = this_sep_set_weight   # could it be get_cost() as well ? 
                    
                 
                    # update with the forward 
                    for c in self.bins: 
                                                
                        if c >= c_weight:
                            
                            self.P[r][node][self.index_bins[c]] = benefits[r][node]
                            # parent benefit is always left unmodified to keep track of the benefit at the parent (above) 
                            self.parent_benefit[r][node][self.index_bins[c]] = benefits[r][node]
                            

                if node in self.leaves: 
                 
                      
                        back_node = node 
                                                
                        # p_v = node 
                        while True:
                                                     
                            # we need to do bacward step 
                            self.backward(visited, back_node, r) # does not return but updating matrix P 
                            
                            # and update the parent
                            back_node = self.parents[back_node] 
                            
                            # when exiting from here we go to the next iteration of the foor lopop 
                            # this can also be replaced with the numbers since we visit nodes in dfs order 
                            if back_node == r or len( self.children[back_node].difference(visited) ) != 0: 
                                break
         
             
                to_visit.update(self.children[node]) 
                to_visit = to_visit.difference({node}) 
    
            
    
    
    
    def sosp_multiple_sp(self): 

        '''
        SOSP algorithm to find optimal K shortcut potentials
        must get as input optimal shortcut potentials for each root-weights
        @return set of optimal K shortcut potentials
        '''        
        
        visited = set()    
        self.D_S =  defaultdict(lambda: defaultdict(set))
        self.w_res_1 = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.w_res_2 = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.w_res_2_ch = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
     
        to_be_processed = [x for x in self.leaves] 
        
        visited.update(set(self.leaves))
        
        # this guarantees that we end when we have visited bottom up all the nodes 
        while len( visited ) < len( self.all_cliques ): 
            
            node = to_be_processed.pop() 
         
            
            visited.add(node) 
            
            # if all che the children have been processed we can move to next node 
            if len(visited.intersection(self.children[ self.parents[node] ])) ==  len(self.children[ self.parents[node] ]): 
                to_be_processed.append(self.parents[node]) 
                
                
            if node not in self.leaves: 
                
                verbose = node == 11
          
                ''' computing H_1 for node'''
                     
                previous_children = set() 
                
                for ch in self.children[node]: 
                    
               
                    dict_optimal_weights_H_1 = dict() 
                    
                    # now we consider the summation 
                    H_1_this_children = [0 for i in self.bins] 
                    
                    for c_1 in self.bins: 
                        
                        for c_2 in self.bins: 
                            
                            # condition in the for loop for c_2 is equivalent to: 
                            if c_1 + c_2 <= self.K: 
                           
                                # we need to find the "c" value that corresponds to the sum 
                                for c_tmp0 in self.bins: 
                                    if c_tmp0 >= c_1 + c_2: 
                                        c_sum = c_tmp0
                                        break 
                               
                                #now the value of c holds the correct value 
                                s = self.H[node][self.index_bins[c_1]] + self.H[ch][self.index_bins[c_2]]
                                
                                
                                if s  > H_1_this_children[self.index_bins[c_sum]]: # new maximum for this children and c_sum 
                                    
                                    H_1_this_children[self.index_bins[c_sum]] = s
                                    
                                    dict_optimal_weights_H_1[c_sum] = (c_1, c_2) 
                                    
                                 
                    
                    for c_c in self.bins: 
                        
                       
                        if  H_1_this_children[self.index_bins[c_c]] > self.H[node][self.index_bins[c_c]]: 
                            
                     
                            self.H[node][self.index_bins[c_c]] = H_1_this_children[self.index_bins[c_c]]
                            
                            c_1, c_2 = dict_optimal_weights_H_1[c_c]
                            
                            
                            self.w_res_1[node][ch][self.index_bins[c_c]] = c_2  
                            
                            
                            for prev_child in previous_children: 
                                
                                    
                                self.w_res_1[node][prev_child][self.index_bins[c_c]] = self.w_res_1[node][prev_child][self.index_bins[c_1]]
                                
                    
                
                    # previous children added to the set to iterate for residual weights            
                    previous_children.add(ch)
                
                ''' computing H_2 for node'''                    
                # for each c we can have a different solution
                for c in self.bins: 
                    
                    visited_D_S = set() 
                    
                    # we consider the successor of the shortcut potential rooted at node 
                    for edge in self.seps[node][self.index_bins[c]]: 
                        
                        path_node = edge[1] 
                        visited_D_S.add( path_node )
                        
                        # bottom node is the children of the new shortcut potentials 
                        self.D_S[node][self.index_bins[c]].add( edge[0] )
                        
                        # add to @visited_D_S all the nodes in the path between parent of bottom node in optimal 
                        # separator and the root (node) 
                        
                        while path_node != node: 
                            
                            path_node = self.parents[path_node]
                            visited_D_S.add( path_node ) 
                            
                    # for each node in visited_D_S if it has children that are not in a path computed above 
                    # such children must be added to self.D_S because it can be the root of a descendant sp 
                    for path_node in visited_D_S: 
                        
                        for children_path_node in self.children[path_node]: 
                            
                            if children_path_node not in visited_D_S: 
                                
                                self.D_S[node][self.index_bins[c]].add( children_path_node )
                                
                    
                    
                    #self.H_2 must be initialized as 
                    self.H_2[node][self.index_bins[c]] = self.P_final[node][self.index_bins[c]]
                   
               
                # then we consider the combination of given values of c and children 
                for c in self.bins: 
                    
                    # now we simply sum the new benefits of the others in left to right fashion as before  
                    previous_children = set() 
                    
                    # this holds the current total for the children of the rooted shortcut potential 
                    H_2_ch = [0 for _ in self.bins]
                    
                    #now we consider l-r all children to try to improve it 
                
                    # for each child 
                    for ch in self.D_S[node][self.index_bins[c]]: 
                        
                        H_2_this_ch = [0 for _ in self.bins]
                        
                        # combinations of weights for the (previous optimal and new) children 
                        for c_1 in self.bins: 
                                              
                            for c_2 in self.bins:    
                                
                                # each time we are updating the optimal in H_2 
                                #children totla weight (new children and H_2_ch) 
                                weight_c = c_1 + c_2
                                
                                # the toal weight of the current solution is @weight_c + @c 
                                if weight_c + c <= self.K: 
                                    # new optimal obtained with the new children s
                                    # if H_2_ch[c_1] + self.P_final[node][c] + self.H[ch][c_2]  >  H_2_ch[weight_c]+ self.P_final[node][c]:
                                    
                                    s = H_2_ch[self.index_bins[c_1]] + self.H[ch][self.index_bins[c_2]] 
                                    
                                    # we need to find the c corresponding to weight_c 
                                    for c_tmp in self.bins: 
                                        if c_tmp >= weight_c: 
                                            c_weight = c_tmp
                                            break
                                        
                                    
                                    if s  >  H_2_ch[self.index_bins[c_weight]]: 

                                        # we keep the optimal children combination
                                        H_2_this_ch[self.index_bins[c_weight]] = s
                                        
                                        
                                        # keep residual weight for children combination 
                                        self.w_res_2_ch[node][ch][self.index_bins[c_weight]] = c_2 
                                        
                                     
                                        for prev_child in previous_children:
                                            self.w_res_2_ch[node][prev_child][self.index_bins[c_weight]] = self.w_res_2_ch[node][prev_child][self.index_bins[c_1]]
                                            
                                       
                                        s2 =  H_2_this_ch[self.index_bins[c_weight]] + self.P_final[node][self.index_bins[c]] 
                                        
                                        # find the "c" corresponding to weight_c + c 
                                        
                                        for c_tmp_2 in self.bins: 
                                            if c_tmp_2 >= weight_c + c: 
                                                c_weight_c_c = c_tmp_2
                                                break 
                                      
                                        if s2 > self.H_2[node][self.index_bins[c_weight_c_c]]: 
                                            
                                            self.H_2[node][self.index_bins[c_weight_c_c]] = s2 
                                       
                                            self.w_res_2[node][node][self.index_bins[c_weight_c_c]] = c # the first key must be the root node 
                                            
                                         
                                            self.w_res_2[node][ch][self.index_bins[c_weight_c_c]] = c_2 

                                            for prev_child in previous_children: 
                                                
                                        
                                                self.w_res_2[node][prev_child][self.index_bins[c_weight_c_c]] = self.w_res_2_ch[node][prev_child][self.index_bins[c_1]]

                        
                        
                        for c_c in self.bins:
                            if H_2_this_ch[self.index_bins[c_c]]>H_2_ch[self.index_bins[c_c]]: 
                                H_2_ch[self.index_bins[c_c]] = H_2_this_ch[self.index_bins[c_c]]
                            
                            
                        
                        # next children             
                        previous_children.add( ch ) 
                              
                
                for c in self.bins:     
                     
                    # last check if we exceed H_2 for the same weight 
                    
                    
                   
                    if self.H_2[node][self.index_bins[c]] > self.H[node][self.index_bins[c]]:
                        
                     
                        # solution includes @node 
                        self.I_final[node][self.index_bins[c]] = 1 
                        
                        # in the end H contains the overall optimal solution ( we are overweiting H if H_2 is larger) 
                        self.H[node][self.index_bins[c]] = self.H_2[node][self.index_bins[c]] 
    
  
   
        # H[0][self.K] holds the optimal benefit for the overall partitioning (assuming the root is 0)     
        return self.H[0][self.index_bins[self.K]]         

    
    
    
    
    

    def reconstruct_solution_sosp(self, node, c): 

        ''' recursively reconstruct the final optimal packing of shortcut potential
        @params: 
            
        node : current node (first function call this is the root 0) 
        current_weight: current weight (first function call this is self.K) 
        '''
        

                    
        if node in self.leaves: 
            
            self.visited_leaves.add(node) 
            
            # making sure it ends if we have visited all the leaves 
            if len(self.leaves.difference(self.visited_leaves)) == 0:  
                return None 


        
                    
        elif self.I_final[node][self.index_bins[c]] == 1: 
            
            # append optimal subtree 
            rooted_weight =  self.w_res_2[node][node][self.index_bins[c]]
            # key is the top separator while all the bottom separators are the value 
    
          
            
            
            
            
            if len(self.seps[node][self.index_bins[rooted_weight]]) > 0: 
                self.optimal_packing[(node , self.parents[node])] = self.seps[node][self.index_bins[rooted_weight]]
            
            
            for ch in self.D_S[node][self.index_bins[rooted_weight]]: 
                self.reconstruct_solution_sosp(ch , self.w_res_2[node][ch][self.index_bins[c]])
                
            
            
        # in this case in is larger H_1
        else:  
            
            for ch in self.children[node]: 
            
                self.reconstruct_solution_sosp(ch , self.w_res_1[node][ch][self.index_bins[c]])

    
                                
                
    def fix_sp(self, x): 
        
        allseps = x.bottom_separators 
        
        sp = ShortcutPotential(self.jointree)
        
        # find top separator 
        allbots = [ ]
        alltops = [ ]
        for sep in allseps: 
            alltops.append(sep[1]) 
            allbots.append(sep[0])
            
        if self.pivot in alltops: 
            sp.top_separator = (0,-1) 
            
        else:
            # find minimum height node 
            top_sep_top = -1
            min_height = float("inf") 
            for bot in allbots: 
                if self.heights[bot] < min_height: 
                    top_sep_top = bot 
                    min_height = self.heights[bot] 
                    
            
            sp.top_separator = (top_sep_top,self.parents[top_sep_top])
            sp.bottom_separators = x.bottom_separators.difference({sp.top_separator})
            sp.potential = Potential()
            sp.scope = self.compute_scope_sp(sp) 
            
            
        return sp 
    
                                    
    def compute_all_shortcut_potentials(self): 
    
        '''
        Compute the shortcut potentials retrieved by our algorithm 
        '''

        self.all_shortcutpotentials = []
        # overall_joint = self.jointree.compute_overall_joint() 

        for edge in self.optimal_packing:
            # new shortcut potential initialization
            sp = ShortcutPotential(self.jointree)
            # top separator 
            sp.top_separator = edge 
            # bottoms separator 
            
           
            
            sp.bottom_separators = self.optimal_packing[edge] 
            # potential 
                
            # sp.potential = self.compute_single_shortcut_potential(self.jointree, set([edge]).union(self.optimal_packing[edge])  ) 
            sp.potential = Potential()
            sp.scope = self.compute_scope_sp(sp) 
            
            self.all_shortcutpotentials.append(  sp   )
            
            #to_be_pickled = self.all_shortcutpotentials
        
        #all_sp_path = "/u/50/ciaperm1/unix/Desktop/CODE_FAST_APPROXIMATION/pickles/all_shortcut_potentials_"  + self.dataset 
        #pickle.dump(to_be_pickled , open(all_sp_path,"wb"))
    
       
        

    
    def compute_scope_sp(self, sp):           
        
        if sp.top_separator!=None:
            allseps = sp.bottom_separators.union({sp.top_separator}) 
        else: 
            allseps = sp.bottom_separators
        
        allvars = set() 
        
        for sep in allseps: 
            
            bottom = sep[0] 
            top = sep[1] 
            
            if top!=-1:
            
                #pot_bottom = self.jointree.potentials[self.all_cliques[bottom].id]
                
                # scope 
                scope_bottom = self.all_cliques[bottom].node_ids
                
               # pot_top = self.all_cliques[top].scope 
                # scope 
                scope_top = self.all_cliques[top].node_ids 
                
                allvars.update( scope_bottom.intersection(scope_top) )
                
            
            
        return allvars     
        
        
  
    def compute_single_shortcut_potential(self, jointree, seps): 
        
        '''
        compute a single shortcut potential from the separators 
        @params
        seps : list of separator potentials for a given shortcut potential to compute (union of top and bottom separators) 
        overall_joint : joint pdf of all the cliques in the tree
        
        '''
        pot = Potential()
    
        return pot 
    

    
    def query(self, steiner_nodes, steiner_parents, steiner_children, queryvars): 
        
        '''  final message passing to
        answer the) current query 
        we start from the lowest nodes and each node send a message to its parent 
        until we have reached the lowest common ancenstor which does not have parents  in the steiner tree 
        '''

        
        # leaves 
        tmp_potentials = copy.deepcopy(self.jointree.potentials)
        
        leaves = set([ x for x in steiner_nodes if x not in steiner_children ])
        
        current_query_vars = set() 
        tobeprocessed = list(leaves)
        visited = set() 
        
        dict_query_vars = defaultdict(set) # this dictionary keeps track of the query variables in the subtree rooted at the key 
            
        while True: 
     
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
          
            if len(visited) == len(steiner_nodes): 
                break 
            
            
             # compute dict_query_vars by merging 
            if node not in leaves:
                for ch in steiner_children[node]: 
                    dict_query_vars[node].update(dict_query_vars[ch])
                 
            dict_query_vars[node].update(  self.all_cliques[node].node_ids.intersection(queryvars)   )
            
            
            current_query_vars = dict_query_vars[node]
            
            
            
            sepset = self.all_cliques[node].get_sep_set(self.all_cliques[steiner_parents[node]])

            if self.jointree.potentials.get(sepset.id): 
                tmp_potentials[self.all_cliques[steiner_parents[node]].id] = PotentialUtil.send_message_online_from_potential( self.jointree,  tmp_potentials[self.all_cliques[node].id] ,  tmp_potentials[sepset.id]   ,tmp_potentials[self.all_cliques[steiner_parents[node]].id], current_query_vars )

                
            else:
                sepset = self.all_cliques[steiner_parents[node]].get_sep_set(self.all_cliques[node])
                tmp_potentials[self.all_cliques[steiner_parents[node]].id] = PotentialUtil.send_message_online_from_potential( self.jointree,  tmp_potentials[self.all_cliques[node].id] , tmp_potentials[sepset.id] ,  tmp_potentials[self.all_cliques[steiner_parents[node]].id], current_query_vars )
      
                
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                tobeprocessed.append(steiner_parents[node]) 
                
        return PotentialUtil.marginalize_for_from_potential(self.jointree, tmp_potentials[self.all_cliques[node].id], [self.jointree.get_bbn_node(x) for x in queryvars])

    
    
    
    
    
    
    def extract_steiner_tree(self, querynodes , pivot): 
    
        '''
        extract Steiner tree 
        param queryvar : set of query variables 
        Returns set of steiner nodes 
        '''
     
        steiner_nodes = set() 
        steiner_edges = set() 
        steiner_parents = dict() 
        steiner_children = defaultdict(list) 
    
        # start with query variables in the steiner tree 
        steiner_nodes.update( {self.dict_var_clique[x] for x in querynodes} )

        # find lowest common ancestor using heights and paths to it 
        heights_queryvars = { k:v for k,v in  self.heights.items() if k in steiner_nodes }

        # highest clique 
        key_min = min( heights_queryvars.keys(), key=(lambda k: heights_queryvars[k]) ) 
        
 
        flag_vars = self.all_cliques[key_min].node_ids.intersection(querynodes)
        steiner_pivot = key_min
    
        # find lowest common ancestor 
        # the stenier tree is found when all the nodes in the dictionary
        # terminates when we have deleted all the nodes but one 
        visited = set() 
        while len(heights_queryvars) > 1: 
            
    
            # find max height 
            key_max = max( heights_queryvars.keys(), key=(lambda k: heights_queryvars[k]) )  
            visited.add(key_max) 
      
            
            steiner_nodes.add( self.parents[key_max] )
    
            steiner_edges.add( (key_max, self.parents[key_max])  )
             # need to check this if self.jointree.parents or self.parents only instead 
    
            steiner_parents[key_max] = self.parents[key_max]
    
            steiner_children[self.parents[key_max]].append(key_max)
    
            # here we are adding the parent so we need to check everytime whether we have this separator as a bottom 
            # separator for a given shortcut potential 
    
            # delete current node with minumum height 
            del heights_queryvars[key_max]
    
            # add parent (in the end all the nodes will have a common parent)
            heights_queryvars[self.parents[key_max]] = self.heights[self.parents[key_max]]
            
            # implement early stopping (if a variable exists in both children and ancestors) 
            if len(flag_vars.difference(self.all_cliques[self.parents[key_max]].node_ids)) == 0:
                                
                # we have also included the top variable we can exit loop  and return 
                # furthermore the highest clique is not achieved hence we remove it from the set of steiner nodes 
                if self.parents[key_max] != key_min:  
                    
                    if key_min in steiner_nodes and len((steiner_nodes.difference({key_min})).difference(visited))==0:
                        steiner_nodes.remove(key_min)
                        steiner_pivot = self.parents[key_max]
                        break 
                    
               # if len(self.children[self.parents[key_max]].intersection(steiner_nodes).difference(visited)) == 0:  
               #     break 
                    
        if len(visited)>0 and self.parents[key_max] != key_min: 
            steiner_pivot = self.parents[key_max] 
            

        return  steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot

    
    
    
    
    @anycache(cachedir='/tmp/anycache.my')
    def compute_lrdpsosp(self , budget): 
    
        start_time = time.time()
    
        # left right dynamic programming 
        self.lrdp_single_sp() 
    
        #pickle_path = "/u/50/ciaperm1/unix/Desktop/Code/pickles/lrdp_sosp_" + dataset
        #pickle.dump(self , open(pickle_path,"wb"))
    
        self.reconstruct_solution_lrdp() 
    
        optimal_benefit = self.sosp_multiple_sp() 
        self.reconstruct_solution_sosp(0,budget) # assuming root is clique 0
    
        print("LRDP SOSP done. Preprocessing time --- " + str( (time.time() - start_time) ) + " --- seconds \n"    )  
    
        self.compute_all_shortcut_potentials() 
        
        return self
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   
        
        
        
        
        
        
        
        
        
        
        
        
        