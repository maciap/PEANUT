import operator 
from collections import defaultdict, Counter, MutableSet 
import numpy as np 
from graph.potential import PotentialUtil, Potential, PotentialEntry 
from SHORTCUT_POTENTIAL import ShortcutPotential
import copy 
import time 
import pickle 
import anycache


class OrderedSet(MutableSet):

    def __init__(self, iterable=None):
        self.end = end = [] 
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)



class QueryProcessor_naive(): 
    
     
    ''' 
    This class is intended to process queries
    without using any materialization. 
    '''
    
    def __init__(self, lrdp_sosp): 
                 
         
        self.jointree = lrdp_sosp.jointree 
        self.lrdp_sosp = lrdp_sosp
        self.all_cliques = lrdp_sosp.all_cliques
        self.heights = lrdp_sosp.heights
        self.parents = lrdp_sosp.parents
        self.children = lrdp_sosp.children
        self.dict_var_clique = self.create_var_clique_mapping()  

        
        
    def create_var_clique_mapping(self): 
      
        '''  create a mapping between each variable and the 
        clique which contains it closest to the root for efficient query processing
        @returns 
        map_clique_var : dict 
        '''
        
        # cliques sorted in order of increasing heights 
        map_clique_var = dict( )
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
    
    
    
    def extract_steiner_tree(self, querynodes , pivot): 
    
        '''
        extract Steiner tree 
        param 
        queryvar : set of query variables 
        pivot : junction tree pivot 
        Returns 
        steiner_nodes: set of steiner nodes 
        steiner_edges : steiner tree edges 
        steiner_parents : steiner tree mapping for parents 
        steiner_children : steiner tree mapping for children 
        steiner_pivot : steiner tree pivot  
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
    
            steiner_parents[key_max] = self.parents[key_max]
    
            steiner_children[self.parents[key_max]].append(key_max)
    
            del heights_queryvars[key_max]
    
            # add parent (in the end all the nodes will have a common parent)
            heights_queryvars[self.parents[key_max]] = self.heights[self.parents[key_max]]
            
            # implement early stopping (if a variable exists in both children and ancestors) 
            if len(flag_vars.difference(self.all_cliques[self.parents[key_max]].node_ids)) == 0:
                                
                if self.parents[key_max] != key_min:  
                    
                    if key_min in steiner_nodes and len((steiner_nodes.difference({key_min})).difference(visited))==0:
                        steiner_nodes.remove(key_min)
                        steiner_pivot = self.parents[key_max]
                        break 
              
        if len(visited)>0 and self.parents[key_max] != key_min: 
            steiner_pivot = self.parents[key_max] 
            

        return  steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot



    def depthOfTree(self, node_i, steiner_children, steiner_nodes): 
        
        if node_i not in steiner_nodes: 
            return 0 
            
        maxdepth = 0 
        for ch_i in steiner_children[node_i]: 
            this_depth = self.depthOfTree(ch_i, steiner_children, steiner_nodes) 
            maxdepth = max(maxdepth, this_depth)
        
        return maxdepth + 1 
        
        

    def compute_diameter(self, node_i , steiner_children, steiner_nodes):
        
        ''' compute the diameter of a function 
        params: 
            steiner_children: children of steiner tree 
            steiner_nodes: nodes in steiner tree '''
        
        max1 = 0 
        max2 = 0 
        
        if node_i not in steiner_nodes:
            return 0 
        
        for ch_i in steiner_children[node_i]: 
            h = self.depthOfTree(ch_i, steiner_children, steiner_nodes)
            
            if h > max1: 
                max2 = max1 
                max1 = h 
                
            elif h > max2:
                max2 = h 
                
        maxChildDia = 0
        for ch_i in steiner_children[node_i]: 
            this_diam = self.compute_diameter(ch_i , steiner_children, steiner_nodes)
            maxChildDia = max(maxChildDia, this_diam)
            
        
        return max( maxChildDia, max1 + max2 + 1 ) 
        

    def create_map(self): 
        
        ''' create map between each bottom separator of a
        shortcut potential and query variables in the 
        subtree below it'''
        
     
        leaves = list(self.lrdp_sosp.leaves)
      
        tobeprocessed = leaves 
        visited = set() 
        
        self.sep_to_node_map = defaultdict(set) 
        
        while True: 
             
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
            
        
            
            if len(visited) == len(self.all_cliques): 
                break 
            
            # update separator to node map with this clique 
            self.sep_to_node_map[node].update(self.all_cliques[node].node_ids)
                        
            if node not in leaves: 
                for ch in self.children[node]: 
                    self.sep_to_node_map[node].update(  self.sep_to_node_map[ch]  )
                    
                    
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(self.children[ self.parents[node] ])) ==  len(self.children[ self.parents[node] ]): 
                tobeprocessed.append(self.parents[node]) 
             
    

    
    
    def compute_total_skept(self, queryvars, steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot):

        '''compute total amount of skept operations
        for a given query'''
        
  
        total_skept_this_query = 0 
        
        leaves = [ x for x in steiner_nodes if x not in steiner_children ] 

        for leaf in leaves:      
            
            #current_query_vars = set() 
            
            node = leaf 
            
            # add leaf cost 
        
            while node!=steiner_pivot: 
                
                current_query_vars = self.sep_to_node_map[node].intersection(queryvars)
                               
                par_scope = self.all_cliques[self.parents[node]].node_ids
                separator_scope = self.all_cliques[node].node_ids.intersection(par_scope)
                vars_to_keep = separator_scope.union(current_query_vars) 
                if len(self.all_cliques[node].node_ids.difference(vars_to_keep)) > 0: 
                    
                    this_message = 1 
                    node_scope = self.all_cliques[node].node_ids 
                    for var in node_scope: 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                        
                    #current_query_vars.update(self.all_cliques[node].node_ids.intersection(queryvars))
                
                    total_skept_this_query += this_message
               
                                
                current_query_vars = self.sep_to_node_map[node].intersection(queryvars)
                
                this_message = 1 
                node_scope = self.all_cliques[node].node_ids 
                par_scope = self.all_cliques[steiner_parents[node]].node_ids 
                int_scope = par_scope.union(node_scope).union(current_query_vars)  
                
                for var in int_scope: 
                    
                    this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
               
                total_skept_this_query += this_message
                
            
        return total_skept_this_query
        
    
    
    
    def compute_total_skept_jt(self, queryvars, steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot):

        '''compute total amount of skept operations
        for a given query'''
        
        
        total_skept_this_query = 0 
        
        leaves = [ x for x in steiner_nodes if x not in steiner_children ] 
        
   
        tobeprocessed = leaves 
        visited = set() 
        dict_query_vars = defaultdict(set)
        
        # we should check whether the queryvars are a subset of the variables in the query 
        f = False 
        min_l = float("inf") 
        for nod in range(len(self.all_cliques)):
            this_scope = self.all_cliques[nod].node_ids 
            if len(this_scope.intersection(queryvars)) == len(queryvars):
                f = True 
                # compute message size 
                this_message = 1 
                node_scope = this_scope
                for var in node_scope: 
                    this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
        
                # update minimum message size 
                if this_message < min_l: 
                    node = nod
                    min_l = this_message 
                    top_scope = self.all_cliques[nod].node_ids 
                    
                    
                if nod in self.parents and self.parents[nod]!=-1: 
                    # also check parent and children 
                    message_parent = 1 
                    parent_scope = self.all_cliques[self.parents[nod]].node_ids  
                    int_scope = node_scope.intersection(parent_scope)
                    
                    if len(int_scope.intersection(queryvars)) == len(queryvars):
            
                        for var in int_scope: 
                            message_parent = message_parent * len(self.jointree.get_bbn_node(var).variable.values) 
                            
                        if message_parent < min_l: 
                            min_l = message_parent 
                            top_scope = int_scope 
                    
            
        if f:     
            if len(queryvars) == len(top_scope): 
                return 0
            else: 
                return  min_l
            
        
        while True: 
             
        
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
        
            
             # compute dict_query_vars by merging 
            if node not in leaves:
                for ch in steiner_children[node]: 
                    dict_query_vars[node].update(dict_query_vars[ch])
                 
            # both if leaf and not 
            # finally add possixble extra query variables that are in this node 
            dict_query_vars[node].update(  self.all_cliques[node].node_ids.intersection(queryvars)   )
            
            current_query_vars = dict_query_vars[node]
            
            if len(visited) == len(steiner_nodes): 
                break 
            
            # we are sending a message to this node to its parent, we need to add clique size with variables as well 
 
            if node in leaves: 
                # check wheter there are redundant variables 
                node_scope = self.all_cliques[node].node_ids  
                par_scope = self.all_cliques[steiner_parents[node]].node_ids
                separator_scope = par_scope.intersection(node_scope) 
                vars_to_keep = current_query_vars.union(separator_scope)
                if len( node_scope.difference(vars_to_keep) ) > 0: 
                    # in this case we want to sum out 
                    # compute message and add to total 
                    this_message = 1 
                    for var in node_scope: 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
                    total_skept_this_query += this_message
                
            else:
                
                for ch in steiner_children[node]:
                    
                    ch_query_vars = dict_query_vars[ch] 
                    ch_separator_vars = self.all_cliques[ch].node_ids.intersection(self.all_cliques[node].node_ids) 
                    ch_vars_to_keep = ch_query_vars.union(ch_separator_vars)
                    ch_node_scope = ch_vars_to_keep.union(self.all_cliques[node].node_ids)
                    
                    this_message = 1 
                    for var in ch_node_scope: 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                        
                    total_skept_this_query += 2 * this_message
                    
                # now we sum out if there are redundant variables only 
                par_scope = self.all_cliques[steiner_parents[node]].node_ids
                separator_scope = par_scope.intersection(self.all_cliques[node].node_ids) 
                vars_to_keep_node = dict_query_vars[node].union(separator_scope)
                
                if len( self.all_cliques[node].node_ids.difference(vars_to_keep_node) ) > 0: 
                    
                    this_message = 1 
                    for var in self.all_cliques[node].node_ids.union(dict_query_vars[node]): 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values) 
                        
                    total_skept_this_query += this_message
                 
                
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                tobeprocessed.append(steiner_parents[node]) 
                
              
        for ch in steiner_children[node]:
                    
            ch_query_vars = dict_query_vars[ch] 
            ch_separator_vars = self.all_cliques[ch].node_ids.intersection(self.all_cliques[node].node_ids) 
            ch_vars_to_keep = ch_query_vars.union(ch_separator_vars)
            ch_node_scope = ch_vars_to_keep.union(self.all_cliques[node].node_ids)
            
            this_message = 1 
            for var in ch_node_scope: 
                this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values) 
                
            total_skept_this_query += 2 * this_message
                    
        
        if len(self.all_cliques[node].node_ids.difference(queryvars)) > 0: 
            this_message = 1 
            for var in self.all_cliques[node].node_ids.union(queryvars): 
                this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
            total_skept_this_query += this_message
            
                
        
        return total_skept_this_query
        
        

class QueryProcessor(): 
    
     
    ''' 
    Process queries with PEANUT  
    '''
    
    def __init__(self, lrdp_sosp): 
                 
         
        self.lrdp_sosp = lrdp_sosp
        self.jointree = lrdp_sosp.jointree 
        self.all_cliques = lrdp_sosp.all_cliques
        self.heights = lrdp_sosp.heights
        self.parents = lrdp_sosp.parents
        self.children = lrdp_sosp.children
        self.all_shortcutpotentials = lrdp_sosp.all_shortcutpotentials
        #self.query = queries
        self.dict_var_clique = self.create_var_clique_mapping()  # move this, heights and all_cliques to junction tree 
        #first we need to find the cliques that are in query 
        # we necessarily need a dictionary to map variables from cliques and viceversa 

        
        
    def create_var_clique_mapping(self): 
      
        ''' this function create a mapping between each variable and the 
        clique which contains it closest to the root for efficient query processing
        @returns 
        map_clique_var : dict 
        '''
        
        # cliques sorted in order of increasing heights 
        map_clique_var = dict( )
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

        

    def create_map(self): 
        
        ''' create map between each bottom separator of a shortcut potential and query variables in the 
        subtree below it'''
        
        # visit all the nodes in a bottom_up fashion 
        #leaves = [ x for x in range(len(self.all_cliques)) if x not in self.children.keys() ] 
        
        leaves = list(self.lrdp_sosp.leaves)
      
        tobeprocessed = leaves 
        visited = set() 
        
        self.sep_to_node_map = defaultdict(set) 
        
        while True: 
             
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
            
        
            
            if len(visited) == len(self.all_cliques): 
                break 
            
            # update separator to node map with this clique 
            self.sep_to_node_map[node].update(self.all_cliques[node].node_ids)
                        
            if node not in leaves: 
                for ch in self.children[node]: 
                    self.sep_to_node_map[node].update(  self.sep_to_node_map[ch]  )
                    
                    
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(self.children[ self.parents[node] ])) ==  len(self.children[ self.parents[node] ]): 
                tobeprocessed.append(self.parents[node]) 
             
                

    
    
    def find_sp(self, steiner_edges,  queryvars, steiner_nodes, steiner_pivot, steiner_parents, steiner_children):
        
        '''manipulates the extracted steiner tree
        in order to check whehter there exists a useful shortcut potential and if there is replace it into the steiner 
        tree so as to obtain a reduced steiner tree 
        @returns 
        steiner_nodes:  set, updated steiner nodes
        steiner_parents: dict, updated steiner_parents
        map_id_sp: dict, map from id to sp object 
        
        
        '''
       
        i = len(self.all_cliques) + 1 
        # we create a correspondence between shortcut
        # potential and identifier 
        map_id_sp = dict() 
        below_nodes_dict = dict() 
        above_nodes_dict = dict() 
        skept_nodes_dict = dict() 
        bool_dict = dict() 
        below_nodes_dict_inv = dict() 
        
        total_skept = 0 
            
        for sp in self.all_shortcutpotentials: 
                        
            # all separators 
            allseps = sp.bottom_separators.union({sp.top_separator}) 
            
            
            inters = steiner_edges.intersection(allseps)
         
            # useful shortcut potentials have at least one edge in the intersection 
            if len( inters )  > 0: 

                vars_inters = set() 
                for t in inters: 
                    vars_inters.update( self.all_cliques[t[0]].node_ids.intersection( self.all_cliques[t[1]].node_ids )) 
                
                toberemoved = set() 
                to_be_marginalized = []
                    
                if sp.top_separator in inters: 
    
                    bottom_seps = inters.difference(set([sp.top_separator]))

                    skept_nodes = OrderedSet([])
                    
                    for sep in bottom_seps:
                        
                        node = sep[0] 
                                        
                        this_sep_skept_nodes = OrderedSet([])
                        
                        while node != sp.top_separator[0]: 
    
                            node = self.parents[node] 
                            
                            this_sep_skept_nodes.add(node) 
                               
                            if len(self.all_cliques[node].node_ids.difference(vars_inters).intersection(queryvars)) > 0: 
    
                                # toberemoved 
                                toberemoved.add( sep )
                                                                    
                                # once we have noticed that we can remoove this separator 
                                this_sep_skept_nodes = set()                                     
                                break 
                  
                        for na in this_sep_skept_nodes:
                            skept_nodes.add(na)
                  
                else:
                    
                    
                    self.skept_nodes_with_pivot = OrderedSet([]) 
                    
                    out = self.dfsUtil_findsp(steiner_pivot, steiner_parents, steiner_children, steiner_pivot, set(), inters, sp, queryvars)
              
                    if out == 0: 
                        
                        # remove them all 
                        toberemoved = inters
                        
                    else:
                        
                        # remove none 
                        toberemoved = set() 
                        skept_nodes = self.skept_nodes_with_pivot
                        
        
                int_seps = inters.difference( toberemoved )     
                
                if sp.top_separator in int_seps: 
                    
                    # check that there is still something useful 
                    if len(int_seps) > 1: 
                        
                        below_nodes = set() 
                        above_nodes = set() 

                         # flatten set of separators 
                        for t in int_seps: 
                            below_nodes.add(t[0])
                            above_nodes.add(t[1])
                            to_be_marginalized.extend( self.all_cliques[t[0]].node_ids.intersection( self.all_cliques[t[1]].node_ids )) 
                        
                        
                        pot =  PotentialUtil.marginalize_for_from_potential(self.jointree,  sp.potential ,  [self.jointree.get_bbn_node(x) for x in to_be_marginalized]) 
                        
                        sp_marginalized = ShortcutPotential( self.jointree ) 
                        sp_marginalized.bottom_separators = int_seps.difference({sp.top_separator})
                        sp_marginalized.top_separator = sp.top_separator 
                        sp_marginalized.potential = pot
                         
                        all_vars_this_sp = set() 
                        for t_1 in int_seps: 
                            all_vars_this_sp.update( self.all_cliques[t_1[0]].node_ids.intersection( self.all_cliques[t_1[1]].node_ids )) 
                        
                        sp_marginalized.scope = all_vars_this_sp
                        
                        # need to check here htan all_vars_this_sp scope is lower than the joint scope of the variables 
                        #that we are skipping 
                        
                        # compute total scope that we skip 
                        scope_skept_nodes = set() 
                        for node_sk in skept_nodes: 
                            scope_skept_nodes.update(self.all_cliques[node_sk].node_ids) 
                        
                        if len(scope_skept_nodes) > len(sp_marginalized.scope): 
                            # mapping id to sp 
                            map_id_sp[i] = sp_marginalized 
                            below_nodes_dict[i] = below_nodes 
                            above_nodes_dict[i] = above_nodes 
                            skept_nodes_dict[i] = skept_nodes
                            bool_dict[i] = True 
                            for node in below_nodes: 
                                if node != sp.top_separator[0]: 
                                    below_nodes_dict_inv[node] = i 
             
                else:

                    
                    if len(int_seps) > 0: 
                        
                        below_nodes = set() 
                        above_nodes = set() 
                        
                        # flatten set of separators 
                        for t in int_seps: 
                            below_nodes.add(t[0])
                            above_nodes.add(t[1])
                            to_be_marginalized.extend( self.all_cliques[t[0]].node_ids.intersection( self.all_cliques[t[1]].node_ids )) 
                        
                        pot =  PotentialUtil.marginalize_for_from_potential(self.jointree,  sp.potential ,  [self.jointree.get_bbn_node(x) for x in to_be_marginalized]) 
        
                        sp_marginalized = ShortcutPotential(  self.jointree  )
                        sp_marginalized.bottom_separators = int_seps
                        sp_marginalized.top_separator = None 
                        sp_marginalized.potential = pot
                  
                        all_vars_this_sp = set() 
                        for t_1 in int_seps: 
                            all_vars_this_sp.update( self.all_cliques[t_1[0]].node_ids.intersection( self.all_cliques[t_1[1]].node_ids )) 
                
                        sp_marginalized.scope = all_vars_this_sp
                            
                        # mapping id to sp 
                        
                        scope_skept_nodes = set() 
                        for node_sk in skept_nodes: 
                            scope_skept_nodes.update(self.all_cliques[node_sk].node_ids) 
                        
                        if len(scope_skept_nodes) > len(sp_marginalized.scope): 
                          
                            
                            map_id_sp[i] = sp_marginalized 
                            below_nodes_dict[i] = below_nodes 
                            for node in below_nodes: 
                                below_nodes_dict_inv[node] = i 
                                
                            
                            skept_nodes_dict[i] = skept_nodes
                            bool_dict[i] = False
                            
            i+=1
            
        # we have the set of leaf for each leave with need to go up to the steiner pivot 
        leaves = [ x for x in steiner_nodes if x not in steiner_children ] 
                
        to_be_skept_set = defaultdict(list) 
        
        children_and_parent = defaultdict(list) 
        
        sp_leaves = defaultdict(list) 
        
        all_belows = defaultdict(list) 
        
        previous_skept_dict = defaultdict(list)
        
        for leaf in leaves: 
            
            # start to process a given leaf 
            node = leaf 
            
            while node != steiner_pivot: 
                
                if below_nodes_dict_inv.get(node)!=None: 
                    
                    this_sp = below_nodes_dict_inv[node]
                    
                    this_sp_skept_nodes = set() 
                    
                    to_skip_candidate = steiner_parents[node]
                    
                    previous_skept = to_skip_candidate
                    
                    while to_skip_candidate in skept_nodes_dict[this_sp]:  
                    
                        this_sp_skept_nodes.add(to_skip_candidate)
                         
                        if steiner_parents.get(to_skip_candidate)!=None: 
                            previous_skept = to_skip_candidate
                            to_skip_candidate = steiner_parents[to_skip_candidate]
                        else: 
                            break
                        
                    to_be_skept_set[leaf].append(this_sp_skept_nodes) 
                                
                    children_and_parent[leaf].append( [node, to_skip_candidate] ) 
                    
                    all_belows[leaf].append(node) 
                    
                    sp_leaves[leaf].append(this_sp) 
                    
                    previous_skept_dict[leaf].append(previous_skept)

                    node = previous_skept
                   
                else: 
                                    
                    node = steiner_parents[node] # at this point we need to check whether this node is a below node 
               
            
        allskept = set() 
        
        allvisited = set() 
                
        for leaf in leaves:  
            
            
            flag_not_to_visit = -2
            
            # start to process a given leaf 
            node = leaf 
            
            this_leaf_cnt_sp = 0 
            
            while node != steiner_pivot: 
                                                
                if node in all_belows[leaf]: 
                    
                    
                    if node!=flag_not_to_visit: 
                       allvisited.add(node)
                    
                    this_children = node 
                    
                    this_sp = sp_leaves[leaf][this_leaf_cnt_sp]
                                                                
                    this_parent = children_and_parent[leaf][this_leaf_cnt_sp][1]
                    
                    if bool_dict[this_sp]: 
                        
                        if node!=flag_not_to_visit: 
                        
                            # remove the children from the old parent 
                            old_parent = steiner_parents[this_children]
                            steiner_children[old_parent].remove(this_children)
                            
                            steiner_parents[this_children] = this_sp
                            
                            steiner_children[this_sp].append(this_children) 
                            
                            steiner_parents[this_sp] = this_parent
                            
                            steiner_children[this_parent].append(this_sp) 
                            
                            node = previous_skept_dict[leaf][this_leaf_cnt_sp]
                            
                            flag_not_to_visit = copy.deepcopy(node) 
                            
                            steiner_nodes.add(this_sp) 
                            
                            allskept.update(  to_be_skept_set[leaf][this_leaf_cnt_sp]  )
                            
                            this_leaf_cnt_sp+=1
                            
                        else: 
                            # in this case the children is another shortcut potential 
                            
                            this_children = sp_leaves[leaf][this_leaf_cnt_sp-1]
                            
                            steiner_parents[this_children] = this_sp
                            
                            steiner_children[this_sp].append(this_children) 
                            
                            steiner_parents[this_sp] = this_parent
                            
                            steiner_children[this_parent].append(this_sp) 
                            
                            node = previous_skept_dict[leaf][this_leaf_cnt_sp]
                            
                            flag_not_to_visit = copy.deepcopy(node) 
                            
                            steiner_nodes.add(this_sp) 
                            
                            allskept.update(  to_be_skept_set[leaf][this_leaf_cnt_sp]  )
                            
                            this_leaf_cnt_sp+=1
                            
                        
                    else: 
                        
                        
                        if node!=flag_not_to_visit: 
                            
                        
                            # remove the children from the old parent 
                            old_parent = steiner_parents[this_children]
                            
                            steiner_children[old_parent].remove(this_children)
                            
                            steiner_parents[this_children] = this_sp
                            
                            steiner_children[this_sp].append(this_children) 
                            
                            steiner_pivot = this_sp 
                            
                            steiner_nodes.add(this_sp) 
                            
                            allskept.update(  to_be_skept_set[leaf][this_leaf_cnt_sp]  )
                        
                            # we have reached the pivot 
                            break 
                        
                        else:
                            
                            this_children = sp_leaves[leaf][this_leaf_cnt_sp-1]
                            
                            steiner_parents[this_children] = this_sp
                            
                            steiner_children[this_sp].append(this_children) 
                            
                            steiner_pivot = this_sp 
                            
                            steiner_nodes.add(this_sp) 
                            
                            allskept.update(  to_be_skept_set[leaf][this_leaf_cnt_sp]  )
                        
                            break 
                            
                else:
            
                    if node!=flag_not_to_visit: 
                        allvisited.add(node)
                        
                    node = steiner_parents[node]
            
        
        to_skip = allskept.difference(allvisited) 
        
        for par in to_skip:
         
            
            steiner_parents.pop(par, None) 
            steiner_children.pop(par, None) 
            steiner_nodes = steiner_nodes.difference({par})
            
            keys_children = set([x for x in steiner_children.keys()]) 
            keys_children = keys_children.difference({par})
            
             # furthermore this node must be removed from the other node children list 
            for k in keys_children: 
                if par in steiner_children[k]: 
                    steiner_children[k].remove(par) 
                    # if no steiner children left for a given node adter we removed 
                    #par we shall drop it 
                if len(steiner_children[k]) == 0: 
                    steiner_children.pop(k, None) 
        
        
        # remove duplicates for multiple paths
        for k in steiner_children.keys():
            steiner_children[k] = list(set(steiner_children[k]))
        
        
        # this is needed also for the case in which i have a "this_parent" of a shortcut 
        keys_children = [x for x in steiner_children.keys()]
        for k in keys_children:
            for v in steiner_children[k]:
                if k!=steiner_parents[v]: 
                    steiner_children[k].remove(v) 
            if len(steiner_children[k]) == 0: 
                steiner_children.pop(k, None) 
                
                                   
        return steiner_nodes , steiner_children, steiner_parents ,steiner_pivot ,  map_id_sp , total_skept
    


    def compute_total_skept_jt(self, queryvars, steiner_nodes , steiner_parents , steiner_children , steiner_pivot, map_id_sp):

        '''compute total amount of skept operations for a given query'''
        
        
        total_skept_this_query = 0 
        
        leaves = [ x for x in steiner_nodes if x not in steiner_children ] 
        
   
        tobeprocessed = leaves 
        visited = set() 
        dict_query_vars = defaultdict(set)
        
        # we should check whether the queryvars are a subset of the variables in the query 
        f = False 
        min_l = float("inf") 
        for nod in range(len(self.all_cliques)):
            this_scope = self.all_cliques[nod].node_ids 
            if len(this_scope.intersection(queryvars)) == len(queryvars):
                f = True 

                # compute message size 
                this_message = 1 
                node_scope = this_scope
                for var in node_scope: 
                    this_message = this_message * len( self.jointree.get_bbn_node(var).variable.values )  
        
                # update minimum message size 
                if this_message < min_l: 
                    node = nod
                    min_l = this_message 
                    top_scope = self.all_cliques[nod].node_ids 
                    
                    
                if nod in self.parents and self.parents[nod]!=-1: 
                    # also check parent and children 
                    message_parent = 1 
                    parent_scope = self.all_cliques[self.parents[nod]].node_ids  
                    int_scope = node_scope.intersection(parent_scope)
                    
                    if len(int_scope.intersection(queryvars)) == len(queryvars):
            
                        for var in int_scope: 
                            message_parent = message_parent * len(self.jointree.get_bbn_node(var).variable.values) 
                            
                        if message_parent < min_l: 
                            min_l = message_parent 
                            top_scope = int_scope 
                    
            
        if f:     
            if len(queryvars) == len(top_scope): 
                return 0
            else: 
                return min_l
            
        
        while True: 
            
            
            
        
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
        
            
            if node < len(self.all_cliques): # otherwise it is a shortcut potential 
                # update query variables 
                dict_query_vars[node].update(  self.all_cliques[node].node_ids.intersection(queryvars)  )
                
            else: 
                # if the node is a shortcut potential 
                dict_query_vars[node].update(  map_id_sp[node].scope.intersection(queryvars)  )
                
                
            # compute dict_query_vars by merging 
            if node not in leaves:
                for ch in steiner_children[node]: 
                    dict_query_vars[node].update(dict_query_vars[ch])
                
            current_query_vars = dict_query_vars[node] 
            
            
            if len(visited) == len(steiner_nodes): 
                break 
            
            # we are sending a message to this node to its parent, we need to add clique size with variables as well 
 
            if node in leaves: 
                # check wheter there are redundant variables 
                if node < len(self.all_cliques):
                    node_scope = self.all_cliques[node].node_ids  
                else:
                    node_scope = map_id_sp[node].scope
                    
                if steiner_parents[node] < len(self.all_cliques): 
                    par_scope = self.all_cliques[steiner_parents[node]].node_ids
                    
                else:
                    par_scope = map_id_sp[steiner_parents[node]].scope
                    
                    
                separator_scope = par_scope.intersection(node_scope) 
                vars_to_keep = current_query_vars.union(separator_scope)
                if len( node_scope.difference(vars_to_keep) ) > 0: 
                    # in this case we want to sum out 
                    # compute message and add to total 
                    this_message = 1 
                    for var in node_scope: 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
                    total_skept_this_query += this_message
                
            else:
                
                for ch in steiner_children[node]:
                    
                    ch_query_vars = dict_query_vars[ch] 
                    
                    if ch < len(self.all_cliques):
                        ch_scope = self.all_cliques[ch].node_ids  
                    else:
                        ch_scope = map_id_sp[ch].scope
                    
                    if node < len(self.all_cliques):
                        node_scope = self.all_cliques[node].node_ids  
                    else:
                        node_scope = map_id_sp[node].scope
                    
                    
                    ch_separator_vars = ch_scope.intersection(node_scope) 
                    ch_vars_to_keep = ch_query_vars.union(ch_separator_vars)
                    ch_node_scope = ch_vars_to_keep.union(node_scope)
                    
                    this_message = 1 
                    for var in ch_node_scope: 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values) 
                        
                    total_skept_this_query += 2 * this_message
                    
                    
                if node < len(self.all_cliques):
                    node_scope = self.all_cliques[node].node_ids  
                else:
                    node_scope = map_id_sp[node].scope
                        
                # now we sum out if there are redundant variables only 
                if steiner_parents[node] < len(self.all_cliques):
                    par_scope = self.all_cliques[steiner_parents[node]].node_ids  
                else:
                    par_scope = map_id_sp[steiner_parents[node]].scope
                    
                separator_scope = par_scope.intersection(node_scope) 
                vars_to_keep_node = dict_query_vars[node].union(separator_scope)
                
                if len( node_scope.difference(vars_to_keep_node) ) > 0: 
                    
                    this_message = 1 
                    for var in node_scope.union(dict_query_vars[node]): 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                        
                    total_skept_this_query += this_message
                 
                
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                tobeprocessed.append(steiner_parents[node]) 
                
                
        
        # when we arrive at this point we are at the root which must also be treated like other nodes,
        # receiving messages froma ll children 
        for ch in steiner_children[node]:
                    
            ch_query_vars = dict_query_vars[ch] 
            
            if ch < len(self.all_cliques):
                ch_scope = self.all_cliques[ch].node_ids  
            else:
                ch_scope = map_id_sp[ch].scope
            
            if node < len(self.all_cliques):
                node_scope = self.all_cliques[node].node_ids  
            else:
                node_scope = map_id_sp[node].scope
    
            
            ch_separator_vars = ch_scope.intersection(node_scope) 
            ch_vars_to_keep = ch_query_vars.union(ch_separator_vars)
            ch_node_scope = ch_vars_to_keep.union(node_scope)
            
            this_message = 1 
            for var in ch_node_scope: 
                this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
            total_skept_this_query += 2 * this_message
                    
        
        if len(node_scope.difference(queryvars)) > 0: 
            this_message = 1 
            for var in node_scope.union(queryvars): 
                this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
            total_skept_this_query += this_message
            
        
        return total_skept_this_query



    
    
    def dfsUtil_findsp(self, node, steiner_parents, steiner_children, steiner_pivot, visited, inters, sp, queryvars): 
    
        if node not in visited: 
       
            visited.add(node) 
        

            if node not in steiner_children: 


                path_node = node 
                
                f = False 
                
                while path_node!= steiner_pivot:
                    
                    par = steiner_parents[path_node] 

                    edge = (path_node, par) 

                    if edge in inters: 
                        
                        sep = edge 
                        f = True

                    path_node = par 


                # when we have finished the while loop 
                if not f: 
                    
                    return 0 

                else: 
                    
                  
                    # we go up and check for query variables 
                    par = sep[0] 

                    while par != steiner_pivot: 
                        
                        
                        par = steiner_parents[par]
                        

                        self.skept_nodes_with_pivot.add(par) 

                        if len(self.all_cliques[par].node_ids.difference(sp.scope).intersection(queryvars)) > 0: 
                            # we have found a query variable 
                            
                            return 0


            else:    
                
                
                # if not leaf 
                for child in steiner_children[node]: 
                    # continue the traversing 

                    if child not in visited: 
                        
                        out = self.dfsUtil_findsp(child, steiner_parents, steiner_children, steiner_pivot, visited, inters, sp, queryvars)
                        
                        if out == 0: 
                            return 0

                
                
    def query(self, steiner_nodes, steiner_parents, steiner_children, steiner_pivot,  queryvars, map_id_sp): 
        
        '''  final message passing to
        answer the) current query 
        we start from the lowest nodes and each node send a message to its parent 
        until we have reached the lowest common ancenstor which does not have parents  in the steiner tree 
        '''
        
        # leaves 
        # to be modified tmp_potentials not 
        tmp_potentials = copy.deepcopy(self.jointree.potentials)
        
        for k, v in map_id_sp.items(): 
            tmp_potentials[k] = v.potential

        
        leaves = [ x for x in steiner_nodes if x not in steiner_children ] 
        
     
        tobeprocessed = leaves 
        visited = set() 
        dict_query_vars = defaultdict(set)
        
        while True: 
             
        
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
            
            
            # we terminate when the last node has been added to the nodes that must be visited 
            # this should work becuase we do not need to arrive to the pivot necessarily
            if len(visited) == len(steiner_nodes): 
                break 
            
           
            if node < len(self.all_cliques): # otherwise it is a shortcut potential 
                # update query variables 
                dict_query_vars[node].update(  self.all_cliques[node].node_ids.intersection(queryvars)  )
                
            else: 
                # if the node is a shortcut potential 
                dict_query_vars[node].update(  map_id_sp[node].scope.intersection(queryvars)  )
                
                
            # compute dict_query_vars by merging 
            if node not in leaves:
                for ch in steiner_children[node]: 
                    dict_query_vars[node].update(dict_query_vars[ch])
                
            current_query_vars = dict_query_vars[node] 
                
            # we need to send a message from node to its parent 
            
            if node >= len(self.all_cliques) and steiner_parents[node] >= len(self.all_cliques): 
                
                # both shortcut potentials 
                tmp_potentials[steiner_parents[node]] = PotentialUtil.send_message_online_no_sep_from_potential(self.jointree , tmp_potentials[node] , tmp_potentials[steiner_parents[node]] , query_vars = current_query_vars)
            
            elif node >= len(self.all_cliques) and steiner_parents[node] < len(self.all_cliques):
            # only sending is a shortcut potential 
                tmp_potentials[self.all_cliques[steiner_parents[node]].id] = PotentialUtil.send_message_online_no_sep_from_potential(self.jointree , tmp_potentials[node] , tmp_potentials[self.all_cliques[steiner_parents[node]].id] , query_vars = current_query_vars)
            
            elif node < len(self.all_cliques) and steiner_parents[node] >= len(self.all_cliques):
                #receriver is a sp 
                    
                tmp_potentials[steiner_parents[node]] = PotentialUtil.send_message_online_no_sep_from_potential(self.jointree  , tmp_potentials[self.all_cliques[node].id] ,  tmp_potentials[steiner_parents[node]] , query_vars = current_query_vars)
            
            # no sp 
            else:

                sepset = self.all_cliques[node].get_sep_set(self.all_cliques[steiner_parents[node]])

                if self.jointree.potentials.get(sepset.id): 
                    tmp_potentials[self.all_cliques[steiner_parents[node]].id] = PotentialUtil.send_message_online_from_potential( self.jointree,  tmp_potentials[self.all_cliques[node].id] ,  tmp_potentials[sepset.id]   ,tmp_potentials[self.all_cliques[steiner_parents[node]].id], current_query_vars )

                else:
                    sepset = self.all_cliques[steiner_parents[node]].get_sep_set(self.all_cliques[node])
                    tmp_potentials[self.all_cliques[steiner_parents[node]].id] = PotentialUtil.send_message_online_from_potential( self.jointree,  tmp_potentials[self.all_cliques[node].id] , tmp_potentials[sepset.id] ,  tmp_potentials[self.all_cliques[steiner_parents[node]].id], current_query_vars )

                    
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                tobeprocessed.append(steiner_parents[node]) 
                
            
            # when we arrive here we have to marginalize with respect to the query variables only 
        if steiner_pivot < len(self.all_cliques): 
            return PotentialUtil.marginalize_for_from_potential(self.jointree, tmp_potentials[self.all_cliques[node].id], [self.jointree.get_bbn_node(x) for x in queryvars])
        else:
            return PotentialUtil.marginalize_for_from_potential(self.jointree, tmp_potentials[node], [self.jointree.get_bbn_node(x) for x in queryvars])
            







class QueryProcessor_indsep(): 

    ''' 
    Process queries with INDSEP
    '''
    
    def __init__(self, jointree, all_cliques, indsep):
        
        
        self.jointree = jointree 
        self.all_cliques = all_cliques
        self.indsep = indsep
        self.heights = indsep.heights
        # self.dict_var_clique = self.create_var_clique_mapping()  
        self.node2pot = dict() 
        self.total_skept = 0 
    

        
    def answer_query(self, steiner_nodes, steiner_parents, steiner_children, queryvars): 
    
        ''' 
        subroutine for answer_query_base_case
        
        final message passing to
        answer the) current query at the clique level 
        we start from the lowest nodes and each node send a message to its parent 
        until we have reached the lowest common ancenstor which does not have parents  in the steiner tree 
        '''
        
  
        leaves = [ x for x in steiner_nodes if  len(steiner_nodes.intersection(self.indsep.children[x])) == 0 ] 
        
        tobeprocessed = leaves 
        visited = set() 
            
        dict_query_vars = defaultdict(set)
        
        while True: 
                                     
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
            
            #add import cost 
            node_scope = self.all_cliques[node].node_ids  
            this_message = 1 
            for var in node_scope: 
                this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
        
            self.total_skippable += this_message
        
             # compute dict_query_vars by merging 
            if node not in leaves:
                for ch in steiner_children[node]: 
                    dict_query_vars[node].update(dict_query_vars[ch])
                 
            # finally add possixble extra query variables that are in this node 
            dict_query_vars[node].update(  self.all_cliques[node].node_ids.intersection(queryvars)   )
            
            current_query_vars = dict_query_vars[node]
            
            if len(visited) == len(steiner_nodes): 
                break 
            
            # we are sending a message to this node to its parent, we need to add clique size with variables as well 
            if node in leaves: 
                # check wheter there are redundant variables 
                node_scope = self.all_cliques[node].node_ids  
                par_scope = self.all_cliques[steiner_parents[node]].node_ids
                separator_scope = par_scope.intersection(node_scope) 
                vars_to_keep = dict_query_vars[node].union(separator_scope)
                if len( node_scope.difference(vars_to_keep) ) > 0: 
                    
                    
                    this_message = 1 
                    for var in node_scope: 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
                    self.total_skippable += this_message
                
            else:
                
                for ch in steiner_children[node]:
                    
                    ch_query_vars = dict_query_vars[ch] 
                    ch_separator_vars = self.all_cliques[ch].node_ids.intersection(self.all_cliques[node].node_ids) 
                    ch_vars_to_keep = ch_query_vars.union(ch_separator_vars)
                    ch_node_scope = ch_vars_to_keep.union(self.all_cliques[node].node_ids)
                    
                    this_message = 1 
                    for var in ch_node_scope: 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                        
                    self.total_skippable += 2 * this_message
                    
                # now we sum out if there are redundant variables only 
                par_scope = self.all_cliques[steiner_parents[node]].node_ids
                separator_scope = par_scope.intersection(self.all_cliques[node].node_ids) 
                vars_to_keep_node = dict_query_vars[node].union(separator_scope)
                
                if len( self.all_cliques[node].node_ids.difference(vars_to_keep_node) ) > 0: 
                    
                    this_message = 1 
                    for var in self.all_cliques[node].node_ids.union( dict_query_vars[node] ): 
                        this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                        
                    self.total_skippable += this_message
                 
                
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                tobeprocessed.append(steiner_parents[node]) 
                
                
        
        # when we arrive at this point we are at the root which must also be treated like other nodes,
        # receiving messages froma ll children 
        for ch in steiner_children[node]:
                    
            ch_query_vars = dict_query_vars[ch] 
            ch_separator_vars = self.all_cliques[ch].node_ids.intersection(self.all_cliques[node].node_ids) 
            ch_vars_to_keep = ch_query_vars.union(ch_separator_vars)
            ch_node_scope = ch_vars_to_keep.union(self.all_cliques[node].node_ids)
            
            this_message = 1 
            for var in ch_node_scope: 
                this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
            self.total_skippable += 2 * this_message
                    
        
        if len(self.all_cliques[node].node_ids.difference(queryvars)) > 0: 
            
            this_message = 1 
            for var in self.all_cliques[node].node_ids.union(dict_query_vars[node]): 
                this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
            self.total_skippable += this_message
            
        

                
        
    def answer_query_gm(self,  steiner_nodes, steiner_parents, steiner_children, queryvars): 
        
        tmp_potentials = dict()

        leaves = set([ x for x in steiner_nodes if  len(steiner_nodes.intersection(steiner_children[x])) == 0 ])

        #leaves_tuples = set([ x for x in steiner_tuples if  len(steiner_nodes.intersection(steiner_children[x[0]]))== 0 ])

        queryvars_mapped = set([self.indsep.var_mapping[x] for x in queryvars])
        # intitialize leaf potentials 
   
        for leaf in leaves: 
            tmp_potentials[leaf] = self.node2pot[leaf]

        # we start by processing the leaves 
        tobeprocessed = [x for x in leaves]
        visited = set() 
        #current_query_vars_gm = set()

        dict_query_vars = defaultdict(set)
        
        while True: 

            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 

            # we terminate when the last node has been added to the nodes that must be visited 
            if len(visited) == len(steiner_nodes):
                break 

          
            this_node_vars = set() 
            for lst_ch in self.indsep.range_list[node]: 
                this_node_vars.update(set([x for x in range( lst_ch[0], lst_ch[1] + 1 )]))
            
            #this_node_vars = set([x for x in range(self.indsep.range_list[node][0],  self.indsep.range_list[node][1]+1)]).union(self.indsep.addList[node])

            this_node_vars_mapped = set([self.indsep.var_mapping[x] for x in this_node_vars])
            
            dict_query_vars[node].update(this_node_vars_mapped.intersection(queryvars_mapped))
            
            # compute dict_query_vars by merging 
            if node not in leaves:
                for ch in steiner_children[node]: 
                    dict_query_vars[node].update(dict_query_vars[ch])
                
            current_query_vars_gm = dict_query_vars[node] 

            #get sepset potential 
            for sep_pot_tuple in self.indsep.all_separator_potentials[node]: 
                if sep_pot_tuple[1] == steiner_parents[node]:
                    sep_pot = sep_pot_tuple[2] 
                    break


            tmp_potentials[steiner_parents[node]] = PotentialUtil.send_message_online_from_potential( self.jointree,  tmp_potentials[node] ,    sep_pot   , self.node2pot[steiner_parents[node]] , current_query_vars_gm )
          

            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                tobeprocessed.append(steiner_parents[node]) 


        # when we arrive here we have to marginalize with respect to the query variables only 
        return PotentialUtil.marginalize_for_from_potential(self.jointree, tmp_potentials[node], [self.jointree.get_bbn_node(x) for x in queryvars_mapped])



                
                
    def compute_max_skippable_gm(self,  steiner_nodes, steiner_parents, steiner_children, queryvars): 
        
        tmp_potentials = dict()
        
        leaves = set([ x for x in steiner_nodes if  len(steiner_nodes.intersection(steiner_children[x])) == 0 ])

        #leaves_tuples = set([ x for x in steiner_tuples if  len(steiner_nodes.intersection(steiner_children[x[0]]))== 0 ])

        queryvars_mapped = set([self.indsep.var_mapping[x] for x in queryvars])
        # intitialize leaf potentials 
   
        for leaf in leaves: 
            tmp_potentials[leaf] = Potential()

        dict_query_vars = defaultdict(set)
        
        tobeprocessed = list([ x for x in steiner_nodes if  len(steiner_nodes.intersection(steiner_children[x])) == 0 ]) 
        
        visited = set() 
            
        dict_query_vars = defaultdict(set)
        
        min_size = float("inf") 
        # no message passing case, al the query variables are contained into this steiner index node 
        f = False 
        for node in steiner_nodes:
        
            if len(queryvars_mapped.difference(self.dict_skippable[node])) == 0:  
               f = True
                     
               if len(self.dict_skippable[node]) < min_size :
                    
                   min_size = len(self.dict_skippable[node]) 
                   min_node = node 
                    
                    
        if f: 
            
            if min_size > len(queryvars_mapped): 
                
                this_message = 1 
                for var in self.dict_skippable[min_node]: 
                    this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                
                self.total_skippable += this_message
                

        else: 
            
            
            while True: 
                                         
        
                # we get the first element 
                node = tobeprocessed.pop(0)
                visited.add(node) 
               
                  
                this_node_vars_mapped = set([x for x in self.dict_skippable[node] ])
                
                dict_query_vars[node].update(this_node_vars_mapped.intersection(queryvars_mapped))
                            
                # compute dict_query_vars by merging 
                if node not in leaves:
                    for ch in steiner_children[node]: 
                        dict_query_vars[node].update(dict_query_vars[ch])
                    
                current_query_vars_gm = dict_query_vars[node]
                
                if len(visited) == len(steiner_nodes): 
                    break       
                
                
                if node in leaves: 
                    # check wheter there are redundant variables 
                    node_scope = self.dict_skippable[node]
                    par_scope = self.dict_skippable[steiner_parents[node]]
                    separator_scope = par_scope.intersection(node_scope) 
                    vars_to_keep = dict_query_vars[node].union(separator_scope)
                    if len( node_scope.difference(vars_to_keep) ) > 0: 
                        # in this case we want to sum out 
                        # compute message and add to total 
                        this_message = 1 
                        for var in node_scope: 
                            this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                    
                        self.total_skippable += this_message
                    
                else:
                    
                    for ch in steiner_children[node]:
                        
                        ch_query_vars = dict_query_vars[ch] 
                        ch_separator_vars = self.dict_skippable[ch].intersection(self.dict_skippable[node]) 
                        ch_vars_to_keep = ch_query_vars.union(ch_separator_vars)
                        ch_node_scope = ch_vars_to_keep.union(self.dict_skippable[node])
                        
                        this_message = 1 
                        for var in ch_node_scope: 
                            this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                            
                        self.total_skippable += 2 * this_message
                        
                    # now we sum out if there are redundant variables only 
                    par_scope =  self.dict_skippable[steiner_parents[node]]
                    separator_scope = par_scope.intersection(self.dict_skippable[node]) 
                    vars_to_keep_node = dict_query_vars[node].union(separator_scope)
                    
                    if len( self.dict_skippable[node].difference(vars_to_keep_node) ) > 0: 
                        
                        this_message = 1 
                        for var in self.dict_skippable[node].union(dict_query_vars[node]): 
                            this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                            
                        self.total_skippable += this_message
                 
                
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
                if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                    tobeprocessed.append(steiner_parents[node]) 
                                        
            
            for ch in steiner_children[node]:
                    
                ch_query_vars = dict_query_vars[ch] 
                ch_separator_vars = self.dict_skippable[ch].intersection(self.dict_skippable[node]) 
                ch_vars_to_keep = ch_query_vars.union(ch_separator_vars)
                ch_node_scope = ch_vars_to_keep.union(self.dict_skippable[node])
                
                this_message = 1 
                for var in ch_node_scope: 
                    this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                    
                self.total_skippable += 2 * this_message
                        
        
            if len(self.dict_skippable[node].difference(dict_query_vars[node])) > 0: 
                this_message = 1 
                for var in self.dict_skippable[node].union(dict_query_vars[node]): 
                    this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  
                    
                self.total_skippable += this_message
            
            
       



    def get_steiner_tree_children(self, inode, steiner_nodes):
        
        #steiner_nodes = set([steiner_nodes]) 
        steiner_edges = set() 
        steiner_parents = dict() 
        steiner_seps = dict()
        steiner_children = defaultdict(list)
 
        # find lowest common ancestor using heights and paths to it 
        heights_querynodes = { k:v for k,v in  self.indsep.heights_current_level.items() if k in steiner_nodes }

    
        while len(heights_querynodes) > 1: 
    
            # find max height key
            key_max = max( heights_querynodes.keys(), key=(lambda k: heights_querynodes[k]) )  
    
            steiner_nodes.add( self.indsep.parents_current_level[key_max] )
    
            steiner_edges.add( (key_max, self.indsep.parents_current_level[key_max])  )
             # need to check this if self.jointree.parents or self.parents only instead 
    
            steiner_parents[key_max] = self.indsep.parents_current_level[key_max]
    
            steiner_children[self.indsep.parents_current_level[key_max]].append(key_max)
            
            for sep_pot in self.indsep.all_separator_potentials[key_max]:
                if sep_pot[1] == self.indsep.parents_current_level[key_max]: 
                    steiner_seps[key_max] = sep_pot
                
                    break 
            

            del heights_querynodes[key_max]
    
            heights_querynodes[self.indsep.parents_current_level[key_max]] = self.indsep.heights_current_level[self.indsep.parents_current_level[key_max]]
    
        
        return  steiner_nodes , steiner_parents ,steiner_children ,  steiner_seps 



    def compute_max_skippable_base_case(self,  inode,  query_vars):
                             
                             
        '''subroutine to add a clique 
        and associated separators for the base case'''
        
        steiner_edges = set() 
        steiner_nodes = set() 
        steiner_parents = dict() 
        steiner_children = defaultdict(list)
        allvars = set() # for message passing in the end 
        
        for current_var in query_vars: 
            var_mapped = self.indsep.var_mapping[current_var]
            allvars.add(var_mapped)
            # we need to find the clique which contains this variable 
            for clique in self.indsep.id_to_connected_components_cliques[inode]: 
                if var_mapped in self.all_cliques[clique].node_ids:
                    steiner_nodes.add(clique) 
                    break
     
        heights_querynodes = { k:v for k,v in  self.heights.items() if k in steiner_nodes }

        # highest clique 
        key_min = min( heights_querynodes.keys(), key=(lambda k: heights_querynodes[k]) ) 
        
        
        flag_vars = self.all_cliques[key_min].node_ids.intersection(allvars)
        steiner_pivot = key_min

        # find lowest common ancestor 
        visited = set() 
        while len(heights_querynodes) > 1: 

            # find max height key
            key_max = max( heights_querynodes.keys(), key=(lambda k: heights_querynodes[k]) )  
            
            visited.add(key_max) 

            steiner_nodes.add( self.indsep.parents[key_max] )

            steiner_edges.add( (key_max, self.indsep.parents[key_max])  )
             # need to check this if self.jointree.parents or self.parents only instead 

            steiner_parents[key_max] = self.indsep.parents[key_max]

            steiner_children[self.indsep.parents[key_max]].append(key_max)

            # delete current node with minumum height 
            del heights_querynodes[key_max]

            # add parent (in the end all the nodes will have a common parent)
            heights_querynodes[self.indsep.parents[key_max]] = self.heights[self.indsep.parents[key_max]]

            if len(flag_vars.difference(self.all_cliques[self.indsep.parents[key_max]].node_ids)) == 0:
                            
            # we have also included the top variable we can exit loop  and return 
            # furthermore the highest clique is not achieved hence we remove it from the set of steiner nodes 
                if self.indsep.parents[key_max] != key_min:  
                    
                    if key_min in steiner_nodes and len((steiner_nodes.difference({key_min})).difference(visited))==0:
                        steiner_nodes.remove(key_min)
                        steiner_pivot = self.indsep.parents[key_max]
                        break 
         
        if len(visited)>0 and self.indsep.parents[key_max] != key_min: 
            steiner_pivot = self.indsep.parents[key_max] 
 
        
        self.answer_query(steiner_nodes, steiner_parents, steiner_children, allvars)
        
    
                     


    def create_map(self): 
        
        ''' create map between each bottom separator of a shortcut potential and query variables in the 
        subtree below it'''
        
        # visit all the nodes in a bottom_up fashion 
        # leaves = [ x for x in self.all_cliques if x not in self.children ] 
        
        leaves = list(self.indsep.leaves)
        
        tobeprocessed = leaves 
        visited = set() 
        
        self.sep_to_node_map = defaultdict(set) 
        
        while True: 
             
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
            
            if len(visited) == len(self.all_cliques): 
                break 
            
            # update separator to node map with this clique 
            self.sep_to_node_map[node].update(self.all_cliques[node].node_ids)
                        
            if node not in leaves: 
                for ch in self.indsep.children[node]: 
                    self.sep_to_node_map[node].update(  self.sep_to_node_map[ch]  )
                    
                    
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(self.indsep.children[ self.indsep.parents[node] ])) ==  len(self.indsep.children[ self.indsep.parents[node] ]): 
                tobeprocessed.append(self.indsep.parents[node]) 
                    
                    
                
    
    
    def extract_steiner_tree_sp( self, querynodes ): 
    
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
      
            
            steiner_nodes.add( self.indsep.parents[key_max] )
    
            steiner_edges.add( (key_max, self.indsep.parents[key_max])  )
             # need to check this if self.jointree.parents or self.parents only instead 
    
            steiner_parents[key_max] = self.indsep.parents[key_max]
    
            steiner_children[self.indsep.parents[key_max]].append(key_max)
    
            # here we are adding the parent so we need to check everytime whether we have this separator as a bottom 
            # separator for a given shortcut potential 
    
            # delete current node with minumum height 
            del heights_queryvars[key_max]
    
            # add parent (in the end all the nodes will have a common parent)
            heights_queryvars[self.indsep.parents[key_max]] = self.heights[self.indsep.parents[key_max]]
            
            # implement early stopping (if a variable exists in both children and ancestors) 
            if len(flag_vars.difference(self.indsep.all_cliques[self.indsep.parents[key_max]].node_ids)) == 0:
                                
                # we have also included the top variable we can exit loop  and return 
                # furthermore the highest clique is not achieved hence we remove it from the set of steiner nodes 
                if self.indsep.parents[key_max] != key_min:  
                    
                    if key_min in steiner_nodes and len((steiner_nodes.difference({key_min})).difference(visited))==0:
                        steiner_nodes.remove(key_min)
                        steiner_pivot = self.indsep.parents[key_max]
                        break 
                    
               # if len(self.children[self.parents[key_max]].intersection(steiner_nodes).difference(visited)) == 0:  
               #     break 
                    
        if len(visited)>0 and self.indsep.parents[key_max] != key_min: 
            steiner_pivot = self.indsep.parents[key_max] 
            

        return  steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot

   
                            

    def query(self, inode, query_vars):
        
        
        
        
        ''' 
        this function implements the recursive 
        query processing algorithm of Kanagal (Algorithm I)
        @inode: input index node 
        @params: query_vars: variables which are queried for (in index renaming)
        '''
        
        # base case: if we have reached the leaf partitions we simply need to add those to the final graphical model 
        if self.indsep.index_children[inode] == None:
            #  we do the message passing to create a potential with query and neighbor variables 
            #for clique in self.indsep.id_to_connected_component[node]: 
            
                        
           # print("entering base case for node " + str(inode))
            self.compute_max_skippable_base_case( inode ,  query_vars  )
            pot = Potential() 
            self.node2pot[inode] = pot 

            return ( inode, pot ) 
            # quitting the recursion 
                             
        
        # convert query vars 
        query_vars_mapped = set([self.indsep.var_mapping[x] for x in query_vars])
        

        # we first search the query variables in the children of inode  
        found = dict()
    
        for var in query_vars: 

            # note there might be multiple children this variable belongs to
            for ch in self.indsep.index_children[inode]: 
                
                #for ch in self.indsep.index_children[inode]; :
                range_list_ch = [] 
                for lst in self.indsep.range_list[ch]: 
                    range_list_ch.extend( [x for x in range(lst[0], lst[1] + 1)] )
                    
            
                if var in range_list_ch or var in self.indsep.addList[ch]: 
                    found[var] = ch 


        # now we need to check wheter all the variables are assigned to the same children 
        if len( set( found.values()) ) == 1: 
            
                   
                             
            #c = found[var]  # all the values are the same 
            c = list( found.values() )[0]
            
            answer_sep = None
            # now we check whether all the variables are in the same separator 

            for sep in self.indsep.all_separator_potentials[c]: 
                sep_pot = sep[2]
                
             
                clique_1 = sep[3] 
                clique_2 = sep[4] 
                    
                    
                vars_se_pot = self.all_cliques[clique_1].node_ids.intersection(self.all_cliques[clique_2].node_ids) 
                    
                
                query_vars_set = set( query_vars_mapped ) 

                inters = vars_se_pot.intersection(query_vars_set)
                if len(inters) == len(query_vars_set):
                    answer_sep = sep_pot 
                    break 

            if answer_sep != None:
                
                
                # query nodes 
                query_bbn_nodes = [self.jointree.get_bbn_node(x) for x in query_vars_mapped]

                # we add a node 
                # sep_pot = PotentialUtil.marginalize_for_from_potential(self.jointree, answer_sep, query_bbn_nodes)
                sep_pot = Potential()
                # we add this node with the marginalized potential 
                self.node2pot[c] = sep_pot
                
                # all the query variables are in the same subpartition , we can update also the parent with the same potential
                self.node2pot[inode] = sep_pot
                
                return (c, sep_pot)

            else: 

                # recurse in c 
                            
                node, pot = self.query(c, query_vars)
                self.node2pot[inode] = pot
                return (node, pot)

                
        # otherwise we create the steiner tree         
        else: 
        
                        
            # initialize steiner_nodes with those that contain variables 
            steiner_nodes = set() 
            for var in found.keys(): 
                
          
                steiner_nodes.add(found[var])
       

            # Steiner tree for the children of @inode 
            steiner_nodes , steiner_parents ,steiner_children ,  steiner_seps = self.get_steiner_tree_children(inode, steiner_nodes) 
            gm = []
            
            #steiner_seps_dest = set([x[1] for x in steiner_seps])
            
               
            # visit nodes 
            for node in steiner_nodes:
                
                neighbors_mapped = set() 
                
                if node in steiner_parents: 
                    
            
                    sep_tup = steiner_seps[node] # this is the parent separator 
                    
                
                    sep = sep_tup[2] # first two elements of the tuples are id and destination node 
                    
                    
                    clique_1 = sep_tup[3] 
                    clique_2 = sep_tup[4] 
                    
                    
                    this_sep_var = self.all_cliques[clique_1].node_ids.intersection(self.all_cliques[clique_2].node_ids) 
                    
                    neighbors_mapped.update(this_sep_var)
                    
                
                # if also is not in steiner_parents 
                    
                # now we add the children separator as well 
                for ch in steiner_children[node]: 
                    
                    spot = steiner_seps[ch]
                        
              
                    clique_1 = spot[3] 
                    clique_2 = spot[4] 
        
                    this_sep_var = self.all_cliques[clique_1].node_ids.intersection(self.all_cliques[clique_2].node_ids) 
                    
                    neighbors_mapped.update(this_sep_var)
                                        
                            
                
                
                neighbors = {self.indsep.var_mapping_inv[x] for x in neighbors_mapped}
                
                range_list_set = set() 
                for lst in self.indsep.range_list[node]: 
                    range_list_set.update( set([x for x in range(lst[0], lst[1] + 1)]) )

              
                I_V =  range_list_set.union(set(self.indsep.addList[node])).intersection(set(query_vars))

                query_and_neighbors = I_V.union(neighbors)
                
                self.dict_skippable[node] = set([self.indsep.var_mapping[x] for x in query_and_neighbors])

                f = False
                min_sp_l = float("inf") 
                for sh_pot in self.indsep.all_shortcut_potentials[node]: 
                 
                    #this_sh_pot_vars = [x[0] for x in sh_pot.potential.entries[0].get_entry_keys()]
                    this_sh_pot_vars = sh_pot.scope
                   
                    if len(neighbors_mapped.difference(this_sh_pot_vars)) == 0: 
                        f = True
                        # when break sh_pot is the shortcut potential that is going to be used 
                        
                        this_message = 1 
                        for var in this_sh_pot_vars: 
                            this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  

                        this_sh_pot_vars_len = this_message
                        if this_sh_pot_vars_len < min_sp_l: 
                            min_sp_l = this_sh_pot_vars_len
                            sh_pot_vars = sh_pot.scope
                            
                            
    
                # if there is null intersection between I_V and query_vars 
                if len(I_V) == 0 and f:
             
                    skept_cliques = self.indsep.id_to_connected_components_cliques[node] 
                  
                    this_neighbours = copy.deepcopy(neighbors_mapped)
                    this_skept_cliques = list(copy.deepcopy(skept_cliques)) 
                
                    
                    self.node2pot[node] = Potential()
                    gm.append((node, sh_pot)) 
                        
                       
                   
                else: 
                                            
                    gm.append(self.query(node, query_and_neighbors))


            
            self.compute_max_skippable_gm(  steiner_nodes, steiner_parents, steiner_children, query_vars )
            pot = Potential() 
            self.node2pot[inode] = pot 
            return (inode, pot)
        
        
        
        
    def query_nosp(self, inode, query_vars):
    
    
        ''' 
        this function implements the recursive 
        query processing algorithm of Kanagal (Algorithm I)
        @inode: input index node 
        @params: query_vars: variables which are queried for (in index renaming)
        '''
        if self.indsep.index_children[inode] == None:
            
            self.compute_max_skippable_base_case( inode ,  query_vars  )
            pot = Potential() 
            self.node2pot[inode] = pot 

        return ( inode, pot ) 
            # quitting the recursion 
                             
        # convert query vars 
        query_vars_mapped = set([self.indsep.var_mapping[x] for x in query_vars])
        
        # we first search the query variables in the children of inode  
        found = dict()
    
        for var in query_vars: 

            # note there might be multiple children this variable belongs to
            for ch in self.indsep.index_children[inode]: 
                
                #for ch in self.indsep.index_children[inode]; :
                range_list_ch = [] 
                for lst in self.indsep.range_list[ch]: 
                    range_list_ch.extend( [x for x in range(lst[0], lst[1] + 1)] )
                    
                if var in range_list_ch or var in self.indsep.addList[ch]: 
                    found[var] = ch 


        # now we need to check wheter all the variables are assigned to the same children 
        if len( set( found.values()) ) == 1: 
            
                   
                             
            #c = found[var]  # all the values are the same 
            c = list( found.values() )[0]
            
            answer_sep = None
            # now we check whether all the variables are in the same separator 

            for sep in self.indsep.all_separator_potentials[c]: 
                sep_pot = sep[2]
                
             
                clique_1 = sep[3] 
                clique_2 = sep[4] 
                    
                    
                vars_se_pot = self.all_cliques[clique_1].node_ids.intersection(self.all_cliques[clique_2].node_ids) 
         
                query_vars_set = set( query_vars_mapped ) 

                inters = vars_se_pot.intersection(query_vars_set)
                if len(inters) == len(query_vars_set):
                    answer_sep = sep_pot 
                    break 

            if answer_sep != None:
                
                
                # query nodes 
                query_bbn_nodes = [self.jointree.get_bbn_node(x) for x in query_vars_mapped]
                
                sep_pot = Potential()

                self.node2pot[c] = sep_pot
                
                self.node2pot[inode] = sep_pot
                
                return (c, sep_pot)

            else: 

                # recurse in c 
                            
                node, pot = self.query_nosp(c, query_vars)
                self.node2pot[inode] = pot
                return (node, pot)

                
        # otherwise we create the steiner tree         
        else: 
        
                        
            # initialize steiner_nodes with those that contain variables 
            steiner_nodes = set() 
            for var in found.keys(): 
          
                steiner_nodes.add(found[var])
       

            # Steiner tree for the children of @inode 
            steiner_nodes , steiner_parents ,steiner_children ,  steiner_seps = self.get_steiner_tree_children(inode, steiner_nodes) 
            gm = []
            
            #steiner_seps_dest = set([x[1] for x in steiner_seps])
            
               
            # visit nodes 
            for node in steiner_nodes:
                
                neighbors_mapped = set() 
                
                if node in steiner_parents: 
                    
                
                    sep_tup = steiner_seps[node] # this is the parent separator 
                    
                
                    sep = sep_tup[2] # first two elements of the tuples are id and destination node 
                    
                    
                    clique_1 = sep_tup[3] 
                    clique_2 = sep_tup[4] 
                    
                    
                    this_sep_var = self.all_cliques[clique_1].node_ids.intersection(self.all_cliques[clique_2].node_ids) 
                    
                    neighbors_mapped.update(this_sep_var)
                    
                
                # if also is not in steiner_parents 
                    
                # now we add the children separator as well 
                for ch in steiner_children[node]: 
                    
                    spot = steiner_seps[ch]
                        
              
                    clique_1 = spot[3] 
                    clique_2 = spot[4] 
        
                    this_sep_var = self.all_cliques[clique_1].node_ids.intersection(self.all_cliques[clique_2].node_ids) 
                    
                    neighbors_mapped.update(this_sep_var)
                                        
                            
                
                
                neighbors = {self.indsep.var_mapping_inv[x] for x in neighbors_mapped}
                
                range_list_set = set() 
                for lst in self.indsep.range_list[node]: 
                    range_list_set.update( set([x for x in range(lst[0], lst[1] + 1)]) )

                # variables in this index node 
                # we must use mapping between query variables with variable renaming 
                
      
                    
                I_V =  range_list_set.union(set(self.indsep.addList[node])).intersection(set(query_vars))

                query_and_neighbors = I_V.union(neighbors)
                
                self.dict_skippable[node] = set([self.indsep.var_mapping[x] for x in query_and_neighbors])

                f = False
                min_sp_l = float("inf") 
                for sh_pot in self.indsep.all_shortcut_potentials[node]: 
                 
                    #this_sh_pot_vars = [x[0] for x in sh_pot.potential.entries[0].get_entry_keys()]
                    this_sh_pot_vars = sh_pot.scope
                   
                    if len(neighbors_mapped.difference(this_sh_pot_vars)) == 0: 
                        f = True
                        # when break sh_pot is the shortcut potential that is going to be used 
                        
                        this_message = 1 
                        for var in this_sh_pot_vars: 
                            this_message = this_message * len(self.jointree.get_bbn_node(var).variable.values)  

                        this_sh_pot_vars_len = this_message
                        if this_sh_pot_vars_len < min_sp_l: 
                            min_sp_l = this_sh_pot_vars_len
                            sh_pot_vars = sh_pot.scope
                            
              
                                            
                gm.append(self.query_nosp(node, query_and_neighbors))

            
            
            #print("entering final gm with steiner nodes : " + str(steiner_nodes))
            self.compute_max_skippable_gm(  steiner_nodes, steiner_parents, steiner_children, query_vars )
            pot = Potential() 
            self.node2pot[inode] = pot 
            return (inode, pot)

    
       
        
    
    
        
class Query(): 
    
    def __init__(self, jointree):
        self.jointree = jointree 

      
    def generate(self, n, probs = None, sizes = None, max_size = 5): 
        ''' generate random queries 
        @ n: query 
        @ probs: probability of each bbn node 
        @ size of each query 
        '''
        
        node_list = self.jointree.get_bbn_nodes()
        
        all_queries_training=[] 
        all_queries_test=[] 
        
        # if sizes not specified we want to randomly sample sizes 
        if sizes == None: 
            # we need to find all the variables in the junction tree 
            
            allsizes = [x for x in range(1, max_size)] 
            sizes = np.random.choice( allsizes, n )
            
        
        if probs == "power": 
            
            
            # exponent hyperparameter fixed to 5 
            
            probs = np.random.power(4 , size = len(node_list))
            tot_prob = np.sum(probs) 
            probs = [x/tot_prob for x in probs]
            
            for i in range(n): 
                all_queries_training.append(np.random.choice(node_list, size = sizes[i], p = probs, replace = False))
            
            
            for i in range(2, 6): 
                for j in range(500): 
                    all_queries_test.append(np.random.choice(node_list, size = i, p = probs , replace = False))
                
        
        
        if probs == "power2": 
            
            # exponent hyperparameter fixed to 5 
            
            probs = np.random.power(2 , size = len(node_list))
            tot_prob = np.sum(probs) 
            probs = [x/tot_prob for x in probs]
            
            for i in range(n): 
                all_queries_training.append(np.random.choice(node_list, size = sizes[i], p = probs, replace = False))
            
            
            for i in range(2, 6): 
                for j in range(500): 
                    all_queries_test.append(np.random.choice(node_list, size = i, p = probs , replace = False))
            
        
        if probs == "uniform": 
         
            for i in range(n): 
                all_queries_training.append(np.random.choice(node_list, size = sizes[i], replace = False))
            
            
            for i in range(2, 6): 
                for j in range(500): 
                    all_queries_test.append(np.random.choice(node_list, size = i, replace = False))
            
                    
  
        return all_queries_training, all_queries_test
   

    def set_probabilties(self, queries):
        '''
        compute probabilities from generated queries
        self.probabilities is a dictionary mapping each variable to a probability of being queried for
        @params: 
        queries: queries to compute probabilities of 
        '''
        
        allvars = [] 
        for query in queries: 
            for var in query: 
                allvars.append(var) 
    
        
        # counts 
        cnt_vars = Counter(allvars) 
        
        # total 
        tot = sum(cnt_vars.values())
        
            # estimated probabilites 
        node_list = self.jointree.get_bbn_nodes()
        for i in range(len(node_list)): 
            
            var = node_list[i] 
            
            if var not in cnt_vars.keys(): 
                var.probability=0 
                
            
            var.probability = cnt_vars[var] / tot 
    



    
        