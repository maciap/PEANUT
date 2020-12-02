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

          

class QueryProcessor(): 
    
     
    ''' 
    Process queries with our approach 
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

        
    
    
    def find_sp(self, steiner_edges,  queryvars, steiner_nodes, steiner_pivot, steiner_parents, steiner_children):
        
        '''manipulates the extracted steiner tree
        in order to check whehter there exists a useful shortcut potential
        @returns 
        steiner_nodes:  set, updated steiner nodes
        steiner_parents: dict, updated steiner_parents
        map_id_sp: dict, map from id to sp object 
        '''
       
        i = len(self.all_cliques) + 1 # we start the count of the cliques from here 
        
        # we create a correspondence between shortcut potential and identifier 
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
    
                # we distinguish the case in which the top separator is in the intersection 
                # or not (which means that the steiner tree pivot is inside the nodes shortcutted by @sp ) 
                
                 # we need to look for query variables 
                
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
    
                            node = self.parents[node]  # queryvars must contain the variables in the query 
    
                            this_sep_skept_nodes.add(node) 
        
                            #this_scope = set([x for x in sp.potential.entry_dicts[0].keys()]) 
                       
                           
                            if len(self.all_cliques[node].node_ids.difference(vars_inters).intersection(queryvars)) > 0: 
    
                                # toberemoved 
                                toberemoved.add( sep )
                                                                    
                                # once we have noticed that we can remoove this separator 
                                this_sep_skept_nodes = set()                                     
                                
                                break 
                            
                        for na in this_sep_skept_nodes:
                            skept_nodes.add(na)
        
                # this is the case in which the pivot is in the cliques that are skept thanks to the sp 
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
                        
                        # we can use this sp 
        
                # we remove those that are not useful due to the presence of query variables 
                int_seps = inters.difference( toberemoved )     
                
         
                
                if sp.top_separator in inters: 
                    
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
                        #sp_marginalized.potential.compute_stride(self.jointree) 
                         
                        sp_marginalized.scope = self.lrdp_sosp.compute_scope_sp(sp)
                            
                        
                        # mapping id to sp 
                        map_id_sp[i] = sp_marginalized 
                        below_nodes_dict[i] = below_nodes 
                        above_nodes_dict[i] = above_nodes 
                        skept_nodes_dict[i] = skept_nodes
                        bool_dict[i] = True 
                        for node in below_nodes: 
                            if node != sp.top_separator[0]: 
                                below_nodes_dict_inv[node] = i 
                    
                    # modify steiner tree 
                    
                    # set as parent 
                    #print("skept nodes") 
                    #print(skept_nodes)
                
                    
                    #for node in below_nodes:
                    #    if node!=sp.top_separator[0]: 
                    #        steiner_parents[node] = i 
                    #        steiner_children[i].append(node)
                            
                            # now node has "i" as only parent and cannot be the children of other nodes as well 
                            # drop from parents - no need 
                            # update steiner children dictionary 
                            
                     #       dict_iter = copy.deepcopy(steiner_children)
                     #       for u in dict_iter: 
                     #           if node in steiner_children[u] and u != i: 
                     #               steiner_children[u].remove(node) 
                     #               if len(steiner_children[u]) == 0: 
                     #                   eqp = steiner_children.pop(u, None) 
                     #                  eqp = steiner_parents.pop(u, None) 
                     #                  steiner_nodes = steiner_nodes.difference({u}) 
                                        
                                    
                        

                    # finally set the parent of the sp 
                    #steiner_parents[i] = sp.top_separator[1]
                    #steiner_children[sp.top_separator[1]].append(i) 
                    # steiner children of top separator should not contain those that we skip thanks to the shortcut potential 
                    #for ch in steiner_children[sp.top_separator[1]]: 
                    #    if ch in skept_nodes: 
                    #        steiner_children[sp.top_separator[1]].remove(ch) 
                            
                            
                    #steiner_children[sp.top_separator[1]] = [x for x in steiner_children[sp.top_separator[1]] if x not in skept_nodes]
                
                    # and to steiner nodes 
                    #steiner_nodes.add(i)
                    
                    # remove nodes 
                    #steiner_nodes = steiner_nodes.difference(skept_nodes) 
                    
                    # remove redundant keys 
                    #for k in skept_nodes: 
                    #    eqp = steiner_parents.pop(k, None) 
                    #    eqp = steiner_children.pop(k, None) 
                    
                        
                        for node in skept_nodes: 
                            total_skept += len( self.jointree.potentials[ self.all_cliques[node].id ].entry_dicts )                
                        
    
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
                        #sp_marginalized.potential.compute_stride( self.jointree )
                        sp.scope = self.lrdp_sosp.compute_scope_sp(sp)
                            
                        
                        # mapping id to sp 
                        map_id_sp[i] = sp_marginalized 
                        below_nodes_dict[i] = below_nodes 
                        for node in below_nodes: 
                            below_nodes_dict_inv[node] = i 
                            
                        
                        skept_nodes_dict[i] = skept_nodes
                        bool_dict[i] = False
                    # set as parent 
                    #for node in below_nodes:
                    #    steiner_parents[node] = i 
                    #    steiner_children[i].append(node)
                        
                     #   dict_iter = copy.deepcopy(steiner_children)
                     #   for u in dict_iter: 
                     #       if node in steiner_children[u] and u != i: 
                     #           steiner_children[u].remove(node) 
                     #           if len(steiner_children[u]) == 0: 
                     #               eqp = steiner_children.pop(u, None)
                     #               eqp = steiner_parents.pop(u, None) 
                     #               steiner_nodes = steiner_nodes.difference({u}) 
                                
                                    
                        
                    # in this case the sp has no steiner parent 
                    # the shortcut potential needs an identifier , to treat like it is a clique 
                    #steiner_nodes.add(i)
                    
                    # in this case we also need to update the pivot 
                        #steiner_pivot = i 
                    
                    # remove nodes 
                    #steiner_nodes = steiner_nodes.difference(skept_nodes) 
                    
                    #for k in skept_nodes: 
                    #    eqp = steiner_parents.pop(k, None) 
                    #    eqp = steiner_children.pop(k, None) 
                        
                    
                    
                        for node in skept_nodes: 
                            total_skept += len( self.jointree.potentials[ self.all_cliques[node].id ].entry_dicts )                
                            
                    
            i+=1
            
    
   
        
        # we have the set of leaf for each leave with need to go up to the steiner pivot 
        leaves = [ x for x in steiner_nodes if x not in steiner_children ] 
        
        leaves_set = set(leaves)
        
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
                    
              
                                        
                    #to_skip_candidate = node 
                    
                    this_sp_skept_nodes = set() 
                    
                    to_skip_candidate = steiner_parents[node]
                    
                    while to_skip_candidate in skept_nodes_dict[this_sp]:  
                        
                    
                        
                        this_sp_skept_nodes.add(to_skip_candidate)
                        
                        previous_skept = to_skip_candidate
                        
                        if steiner_parents.get(to_skip_candidate)!=None: 
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
                        
                            # we have reached the pivot 
                            break 
                            
                        
 
                            
                else:
            
                    if node!=flag_not_to_visit: 
                        allvisited.add(node)
                        
                    node = steiner_parents[node]
                    
                        
                        
    
    

        
        ''' new general approach ''' 
        
        
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
    
        
    
        #now we change the steiner tree in view of the shortcut potentials 
        #for sp_i in map_id_sp: 
            
        #    steiner_nodes.add(sp_i)
            
        #    not_to_skip = set() 
            
            # this is the case in which we have
            # also the parent (not the pivot inside) 
        #    if bool_dict[sp_i]: 
                
        #        print("above nodes dict")
        #        print(above_nodes_dict)
            
                # first think about skept nodes 
        #        for node in skept_nodes_dict[sp_i]: 
                    
        #            print("node") 
        #            print(node) 
                    
        #            this_nodes = copy.deepcopy(steiner_children[node])
                    
        #            print("this-nodes") 
        #            print(this_nodes) 
                    
        #            print("steiner children")
        #            print(steiner_children)
                    
        #            for ch in this_nodes: 
                        
        #               print("ch") 
        #               print(ch) 
                        
        #              if (ch in skept_nodes_dict[sp_i] or ch in below_nodes_dict[sp_i]) and ch not in not_to_skip:
        #                  steiner_children[node].remove(ch) 
                            
                 
                
                
         #           print("steiner children before not to skip") 
         #           print(steiner_children)
         #           # if this skept_node does not have children anymore we need to drop it 
         #          if len(steiner_children[node]) == 0: 
         #              eqp = steiner_children.pop(node, None)
         #               eqp = steiner_parents.pop(node, None) 
         #               steiner_nodes = steiner_nodes.difference({node}) 
         #           else:
         #               print("not to skip before") 
         #               print(not_to_skip)
         #               not_to_skip.add(node)
         #               print("not to skip") 
         #              print(not_to_skip)
            
                # update below nodes parents and sp pot children
         #      for bel_node in below_nodes_dict[sp_i]: 
         #           if bel_node not in skept_nodes_dict[sp_i]: 
         #               steiner_parents[bel_node] = sp_i
         #               steiner_children[sp_i].append(bel_node) 
                
                # update shortcut potential children 
                # for ch_sp_pot in below_nodes_dict[sp_i]: 
                
                # now we set the parent of the sp 
         #       target_node = map_id_sp[sp_i].top_separator[1]
                
                # search target node in the skept nodes of the other sp potentials 
                
         #       flag_found = False  
         #       for sp_j in map_id_sp: 
         #           for skept_node in skept_nodes_dict[sp_j]: 
         #               if skept_node == target_node: 
         #                   # we have found a node matching the target one 
         #                   if len(set(skept_nodes_dict[sp_i]).intersection(set(below_nodes_dict[sp_j]))) > 0: 
         #                       steiner_parents[sp_i] = sp_j 
         #                       steiner_children[sp_j].append(sp_i) 
         #                       for ch in steiner_children[sp_j]: 
         #                           if ch in skept_nodes_dict[sp_i] and ch not in not_to_skip: 
          #                              steiner_children[sp_j].remove(ch)
          #                              
          #                      flag_found = True 
          #                      break 
                
          #      if not flag_found: 
           #         steiner_parents[sp_i] =  target_node
           #         steiner_children[target_node].append(sp_i)
                    
                    # steiner children of @target_node should not contain those that 
                    # are skept 
         ##           for ch_target in steiner_children[target_node]: 
          #              if ch_target in skept_nodes_dict[sp_i] and ch_target not in not_to_skip:
          #                  steiner_children[target_node].remove(ch_target)
    
            
            
           # else:           
                
                # first think about skept nodes 
         #       for node in skept_nodes_dict[sp_i]: 
                    
         #           this_nodes = copy.deepcopy(steiner_children[node])
                    
        #            print("this-nodes") 
       #             print(this_nodes) 
                    
        #           print("steiner children")
        #            print(steiner_children)
                    
        #            for ch in this_nodes: 
                        
        #                print("ch") 
        #                print(ch) 
                        
                        
          #              if (ch in skept_nodes_dict[sp_i] or ch in below_nodes_dict[sp_i]) and ch not in not_to_skip:
         #              
          #                  steiner_children[node].remove(ch)
                    
                    # if this skept_node does not have children anymore we need to drop it 
          #          if len(steiner_children[node]) == 0: 
           #             eqp = steiner_children.pop(node, None)
           #             eqp = steiner_parents.pop(node, None) 
           #             steiner_nodes = steiner_nodes.difference({node})
           #         else:
            #            print("not to skip before") 
            #            print(not_to_skip)
            #            not_to_skip.add(node)
            #            print("not to skip") 
            #            print(not_to_skip)
                        
                
                # update below nodes parents and sp pot children
           #     for bel_node in below_nodes_dict[sp_i]: 
           #         if bel_node not in skept_nodes_dict[sp_i]:
                        
           #             flag_f = False
           #             for sk in skept_nodes_dict[sp_i]: 
           #                 if bel_node in self.children[sk]: 
            #                    flag_f = True 
            #                    break
                                
            #            if flag_f: 
             #               steiner_parents[bel_node] = sp_i
              #              steiner_children[sp_i].append(bel_node) 
                    
                
        
        
        # make sure 
      #  for n in steiner_children: 
      #      for ch_n in steiner_children[n]: 
      ##          if ch_n not in steiner_parents: 
      #             steiner_children[n].remove(ch_n) 
      #          else:
      #              if steiner_parents[ch_n]!=n: 
      #                  steiner_children[n].remove(ch_n) 
        
      #  steiner_children_it = copy.deepcopy(list(steiner_children.keys()))
      #  for n in steiner_children_it: 
      #      if len(steiner_children[n]) == 0: 
      #          steiner_children.pop(n, None)
            

       
    # return steiner_nodes , steiner_children, steiner_parents ,steiner_pivot ,  map_id_sp , total_skept
    


    
    def dfsUtil_findsp(self, node, steiner_parents, steiner_children, steiner_pivot, visited, inters, sp, queryvars): 
    
        if node not in visited: 
       
            visited.add(node) 
        

            if node not in steiner_children: 

                #we have found a leaf and we start to go up and we take the highest separator in inters 

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
        
       
        
      #  current_query_vars = set() 
        
        # we keep a FIFO queue , each time a node has received all message from its children
        # we use a standard list for it which is initialized with the leaves 
    
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
                
                scope = set([x for x in map_id_sp[node].potential.entry_dicts[0].keys()])
                dict_query_vars[node].update(  scope.intersection(queryvars)  )
                
           
                
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
                tmp_potentials[self.all_cliques[steiner_parents[node]].id] = PotentialUtil.send_message_online_no_sep_from_potential(self.jointree , tmp_potentials[node] , tmp_potentials[self.all_cliques[steiner_parents[node]].id] , query_vars = current_query_vars )
            
            elif node < len(self.all_cliques) and steiner_parents[node] >= len(self.all_cliques):
                #receriver is a sp 
                    
                tmp_potentials[steiner_parents[node]] = PotentialUtil.send_message_online_no_sep_from_potential(self.jointree  , tmp_potentials[self.all_cliques[node].id] ,  tmp_potentials[steiner_parents[node]] , query_vars = current_query_vars )
            
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
        
        # leaves 
        # to be modified tmp_potentials not 
        tmp_potentials = copy.deepcopy(self.jointree.potentials)
        
        leaves = [ x for x in steiner_nodes if  len(steiner_nodes.intersection(self.indsep.children[x])) == 0 ] 
        
        #current_query_vars = set() 
        
        # we keep a FIFO queue , each time a node has received all message from its children
        # we use a standard list for it which is initialized with the leaves 
    
        tobeprocessed = leaves 
        visited = set() 
        dict_query_vars = defaultdict(set)
        
        
        while True: 
            
            # we get the first element 
            node = tobeprocessed.pop(0)
            visited.add(node) 
            
            # update query variables 
            # current_query_vars.update(  self.all_cliques[node].node_ids.intersection(queryvars)  )
            
            
           
            # we terminate when the last node has been added to the nodes that must be visited 
            if len(visited) == len(steiner_nodes):
                break 

            
            dict_query_vars[node].update(  self.all_cliques[node].node_ids.intersection(queryvars)  )
            
            # compute dict_query_vars by merging 
            if node not in leaves:
                for ch in steiner_children[node]: 
                    dict_query_vars[node].update(dict_query_vars[ch])
                
            
            current_query_vars = dict_query_vars[node] 
             
            
            # we need to send a message from node to its parent - query variables 
            
            x = tmp_potentials[self.all_cliques[node].id]
            y = tmp_potentials[self.all_cliques[self.indsep.parents[node]].id]
            
            sepset = self.all_cliques[node].get_sep_set(self.all_cliques[self.indsep.parents[node]])
            
            if self.jointree.potentials.get(sepset.id): 
                s = tmp_potentials[sepset.id]
                tmp_potentials[self.all_cliques[self.indsep.parents[node]].id] = PotentialUtil.send_message_online_from_potential(self.jointree, x , s , y , current_query_vars)
                
                                                                                                                    
            else:
                
                sepset = self.all_cliques[self.indsep.parents[node]].get_sep_set(self.all_cliques[node])
                s = tmp_potentials[sepset.id]                                                                                                           
                tmp_potentials[self.all_cliques[self.indsep.parents[node]].id] = PotentialUtil.send_message_online_from_potential(self.jointree, x , s , y , current_query_vars)
                
                    
            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                tobeprocessed.append(steiner_parents[node]) 
                
            # when we arrive here we have to marginalize with respect to the query variables only 
    
        return PotentialUtil.marginalize_for_from_potential(self.jointree, tmp_potentials[self.all_cliques[node].id], [self.jointree.get_bbn_node(x) for x in queryvars])

        
        
                
        
    def answer_query_gm(self,  steiner_nodes, steiner_parents, steiner_children, queryvars): 
        
        tmp_potentials = dict()
        
        #print("steiner nodes") 
        #print(steiner_nodes)
        
        #print("steiner parents") 
        #print(steiner_parents) 
        
        #print("steiner children") 
        #print(steiner_children) 

        leaves = set([ x for x in steiner_nodes if  len(steiner_nodes.intersection(steiner_children[x])) == 0 ])

        #leaves_tuples = set([ x for x in steiner_tuples if  len(steiner_nodes.intersection(steiner_children[x[0]]))== 0 ])

        #print("leaves") 
        #print(leaves) 
        
        queryvars_mapped = set([self.indsep.var_mapping[x] for x in queryvars])

        # intitialize leaf potentials 
        
   
        for node in steiner_nodes: 
            tmp_potentials[node] = self.node2pot[node]
            #print("leaf scope") 
            #print(leaf) 
            #print([x for x in tmp_potentials[leaf].entry_dicts[0].keys()])



        
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
                
                a = set([x for x in tmp_potentials[node].entry_dicts[0].keys()]) 
                queryvars_mapped = queryvars_mapped.intersection(a)
                
                break 

            # we need to send a message from node to its parent 
            # sepset_pot = self.gm_edges[(node, self.gm_parents[node])]  

            # query variables so far 
            # query_vars_current_node =  set([x[0] for x in s.entries[0].get_entry_keys()]).intersection(self.all_query_vars)


            this_node_vars = set() 
            for lst_ch in self.indsep.range_list[node]: 
                this_node_vars.update(set([x for x in range( lst_ch[0], lst_ch[1] + 1 )]))


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


            # get receiver potential 
            #for tup in gm: 
            #    if tup[0] == steiner_parents[node]: 
            #        y_pot = tup[1]
            #        break 
            
            
            tmp_potentials[steiner_parents[node]] = PotentialUtil.send_message_online_from_potential( self.jointree,  tmp_potentials[node] ,    sep_pot   , tmp_potentials[steiner_parents[node]] , current_query_vars_gm )
          

            # if the parent of this node has received all the messages from the children that are also in the steiner tree we can 
            if len(visited.intersection(steiner_children[ steiner_parents[node] ])) ==  len(steiner_children[ steiner_parents[node] ]): 
                tobeprocessed.append(steiner_parents[node]) 


        # when we arrive here we have to marginalize with respect to the query variables only 
        return PotentialUtil.marginalize_for_from_potential(self.jointree, tmp_potentials[node], [self.jointree.get_bbn_node(x) for x in queryvars_mapped])





    def get_steiner_tree_children(self, inode, steiner_nodes):
        
        #steiner_nodes = set([steiner_nodes]) 
        steiner_edges = set() 
        steiner_parents = dict() 
        steiner_seps = set() 
        steiner_children = defaultdict(list) 
        
      

        # find lowest common ancestor using heights and paths to it 
        heights_querynodes = { k:v for k,v in  self.indsep.heights_current_level.items() if k in steiner_nodes }
    
        # highest clique 
        key_min = min( heights_querynodes.keys(), key=(lambda k: heights_querynodes[k]) ) 
    
        # find lowest common ancestor 
        # the stenier tree is found when all the nodes in the dictionary
        # terminates when we have deleted all the nodes but one 
    
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
                    steiner_seps.add(sep_pot)
                    break 
            
       
            # here we are adding the parent so we need to check everytime whether we have this separator as a bottom 
            # separator for a given shortcut potential 
    
            # delete current node with minumum height 
            del heights_querynodes[key_max]
    
            # add parent (in the end all the nodes will have a common parent)
            heights_querynodes[self.indsep.parents_current_level[key_max]] = self.indsep.heights_current_level[self.indsep.parents_current_level[key_max]]
    

        return  steiner_nodes , steiner_parents ,steiner_children ,  steiner_seps 



    def answer_query_base_case(self,  inode,  query_vars):
                             
                             
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
                            
        
        # case in which all the query variables are in the same clique (no message passing needed)
        if len(steiner_nodes) == 1: 
            single_clique = steiner_nodes.pop() 
            return self.jointree.potentials[self.all_cliques[single_clique].id]
                             

        # other steiner nodes 
        # find lowest common ancestor using heights and paths to it 
        heights_querynodes = { k:v for k,v in  self.heights.items() if k in steiner_nodes }

        # highest clique 
        key_min = min( heights_querynodes.keys(), key=(lambda k: heights_querynodes[k]) ) 

        # find lowest common ancestor 
        # the stenier tree is found when all the nodes in the dictionary
        # terminates when we have deleted all the nodes but one 

        while len(heights_querynodes) > 1: 

            # find max height key
            key_max = max( heights_querynodes.keys(), key=(lambda k: heights_querynodes[k]) )  

            steiner_nodes.add( self.indsep.parents[key_max] )

            steiner_edges.add( (key_max, self.indsep.parents[key_max])  )
             # need to check this if self.jointree.parents or self.parents only instead 

            steiner_parents[key_max] = self.indsep.parents[key_max]

            steiner_children[self.indsep.parents[key_max]].append(key_max)

            # here we are adding the parent so we need to check everytime whether we have this separator as a bottom 
            # separator for a given shortcut potential 

            # delete current node with minumum height 
            del heights_querynodes[key_max]

            # add parent (in the end all the nodes will have a common parent)
            heights_querynodes[self.indsep.parents[key_max]] = self.heights[self.indsep.parents[key_max]]


        # run the message passing to compute potential 
        pot = self.answer_query(steiner_nodes, steiner_parents, steiner_children, allvars)
        
        return pot 
                             
                             
   
                            

    def query(self, inode, query_vars):
        
        
        ''' 
        
        this function implements the recursive 
        query processing algorithm of Kanagal (Algorithm I)
        @inode: input index node 
        @params: query_vars: variables which are queried for (in index renaming)
        
        '''
        
        
        
        #print("index parents") 
        #print(self.indsep.index_parents)
        
        #print("index parent this node ") 
        #print(self.indsep.index_parents[inode])
        
        
  
                                     
        # base case: if we have reached the leaf partitions we simply need to add those to the final graphical model 
        if self.indsep.index_children[inode] == None:
            #  we do the message passing to create a potential with query and neighbor variables 
            #for clique in self.indsep.id_to_connected_component[node]: 
            
            
            pot = self.answer_query_base_case( inode ,  query_vars  )
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
                
             
                
                range_list_ch = [] 
                for lst in self.indsep.range_list[ch]: 
                    range_list_ch.extend( [x for x in range(lst[0], lst[1] + 1)] )
                
         
                
                if var in range_list_ch or var in self.indsep.addList[ch]: 
                    found[var] = ch 


        
        # now we need to check wheter all the variables are assigned to the same children 
        if len(set(found.values())) == 1: 
                             
            c = list( found.values() )[0] # all the values are the same 

            answer_sep = None
            # now we check whether all the variables are in the same separator 

            for sep in self.indsep.all_separator_potentials[c]: 
                sep_pot = sep[2]
                vars_se_pot = set([x for x in sep_pot.entry_dicts[0].keys()])
                query_vars_set = set( query_vars_mapped ) 

                inters = vars_se_pot.intersection(query_vars_set)
                if len(inters) == len(query_vars_set):
                    answer_sep = sep_pot 
                    break 

            if answer_sep != None:

                # query nodes 
                query_bbn_nodes = [self.jointree.get_bbn_node(x) for x in query_vars_mapped]

                # we add a node 
                sep_pot = PotentialUtil.marginalize_for_from_potential(self.jointree, answer_sep, query_bbn_nodes)

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

            # visit nodes 
            for node in steiner_nodes:

                # find query variables in this node and neighbours 
                # neig_nodes = self.indsep.all_adjacency_lists[node] 
                neighbors_mapped = set() 
                for sep_tup in self.indsep.all_separator_potentials[node]: 
               
                    
                    steiner_seps_nodes_to_check = [x[0] for x in steiner_seps]
                    
                    if sep_tup[0] in steiner_seps_nodes_to_check: 

                        sep = sep_tup[2] # first two elements of the tuples are id and destination node 
                        this_sep_var = set( [x for x in sep.entry_dicts[0].keys()] )
                  
                        neighbors_mapped.update(this_sep_var)


               
                neighbors = {self.indsep.var_mapping_inv[x] for x in neighbors_mapped}
                
               
                
                
                range_list_set = set() 
                for lst in self.indsep.range_list[node]: 
                    range_list_set.update( set([x for x in range(lst[0], lst[1] + 1)]) )
                    
            
                # variables in this index node 
                # we must use mapping between query variables with variable renaming 
                I_V =  range_list_set.union(set(self.indsep.addList[node])).intersection(set(query_vars))

                query_and_neighbors = I_V.union(neighbors)


           

                f = False
                for sh_pot in self.indsep.all_shortcut_potentials[node]: 
                
                    
                    this_sh_pot_vars = [x for x in sh_pot.potential.entry_dicts[0].keys()]
                    if len(neighbors_mapped.difference(this_sh_pot_vars)) == 0: 
                        f = True
                        # when break sh_pot is the shortcut potential that is going to be used 
                        break 

                # if there is null intersection between I_V and query_vars 
                if len(I_V) == 0 and f:
                    # in this case we can use the shortcut potential to answer the query 
                    # add node - index node id and shortcut potential 
                    
                    self.node2pot[node] = sh_pot.potential
                    gm.append((node, sh_pot)) 
                    
                    skept_nodes = self.indsep.id_to_connected_components_cliques[node] 
                        
                    for clique_node in skept_nodes: 
                        self.total_skept += len( self.jointree.potentials[ self.indsep.all_cliques[clique_node].id ].entry_dicts )
                    
                    
                    # using shortcut potential and quitting the recursion going to the next steiner node 

                else: 
                    
                    gm.append(self.query(node, query_and_neighbors))


            # do the message passing (in the end, it shoud return 9 - root and associated potential)
            
            
            
            pot = self.answer_query_gm(  steiner_nodes, steiner_parents, steiner_children, query_vars )
            self.node2pot[inode] = pot 
            return (inode, pot)
        
           
            






    
class Query(): 
    
    def __init__(self, jointree):
        #self.nodes = nodes 
        self.jointree = jointree 

      
    def generate(self, n, probs = None, sizes = None, max_size = 5): 
        ''' generate random queries 
        @ n query 
        @ probability of each node 
        @ size of each query 
        '''
        
        node_list = self.jointree.get_bbn_nodes()
        all_queries=[] 
        
        # if sizes not specified we want to randomly sample sizes 
        if sizes == None: 
            # we need to find all the variables in the junction tree 
            
            allsizes = [x for x in range(1, max_size)] 
            sizes = np.random.choice( allsizes, n )
            
        
        if probs == "power": 
            
            # exponent hyperparameter fixed to 5 
            
            probs = np.random.power(5 , size = len(node_list))
            tot_prob = np.sum(probs) 
            probs = [x/tot_prob for x in probs]
            
            for i in range(n): 
                all_queries.append(np.random.choice(node_list, size = sizes[i], p = probs))
        
        if probs == "uniform": 
            
            probs = np.random.uniform(0,1, size = len(node_list)) 
            tot_prob = np.sum(probs) 
            probs = [x/tot_prob for x in probs]
            
            for i in range(n): 
                all_queries.append(np.random.choice(node_list, size = sizes[i], p = probs))
            
        
        
        #if probs == "random": 
            
        #    s = np.random.normal(size = 1000)
        #    probs = [x / sum(s) for x in s] # normalize to 0 1 and to sum 1 
            
        #for i in range(n):
        #    all_queries.append( list( np.random.choice( node_list , size = sizes[i],  p = probs )  ) )
            # if probs not specified (probs = None) it assumes equal probabiltiies 
        
        return probs, all_queries
   

    def set_probabilties(self, queries = None, probs = None):
        '''
        compute probabilities from generated queries
        self.probabilities is a dictionary mapping each variable to a probability of being queried for
        @params: 
        queries: queries to compute probabilities of 
        '''
        
        if probs == None:
            
            allvars = [] 
            for query in queries: 
                
            
                break 
                
                for var in query: 
                    allvars.append(var) 
            
            # counts 
            cnt_vars = Counter(allvars) 
            
            # total 
            tot = len(allvars)
            
                # estimated probabilites 
            for var in cnt_vars: 
                var.probability = cnt_vars[var] / tot 
    
        else: 
            
            node_list = self.jointree.get_bbn_nodes()
            for i in range(len(node_list)): 
                
                node = node_list[i] 
                node.probability = probs[i] 
            
            
            
    
        