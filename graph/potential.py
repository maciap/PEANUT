import itertools
from collections import defaultdict

class Potential(object):
    """
    Potential.
    """

    def __init__(self):
        """
        Ctor.
        """
        self.entries = []

    def add_entry(self, entry):
        """
        Adds a PotentialEntry.

        :param entry: PotentialEntry.
        :return: This potential.
        """
        self.entries.append(entry)
        return self

    def get_matching_entries(self, entry):
        """
        Gets all potential entries matching the specified entry.

        :param entry: PotentialEntry.
        :return: Array of matching potential entries.
        """
        return [e for e in self.entries if e.matches(entry)]

    @staticmethod
    def to_dict(potentials):
        """
        Converts potential to dictionary for easy validation.

        :param potentials: Potential.
        :return: Dictionary representation. Keys are entries and values are probabilities.
        """
        return {tup[0]: tup[1] for tup in [pe.get_kv() for p in potentials for pe in p.entries]}

    def __str__(self):
        return str.join('\n', [entry.__str__() for entry in self.entries])
    
    
    
    def create_dict_potential(self):
    
        ''' convert potetial to dictionary structure giving the index of each possible combination of 
        variables and values
        for efficient multiplication 
        to call for each potential when its entries have been added'''
        
        self.d = defaultdict(list) 
    
        for i in range(len(self.entries)): 
    
            this_tuple = tuple( sorted(self.entries[i].entries.items()) )
    
            # update dictionary with this entry   
            for r in range(1, len(this_tuple)+1):
                for k in itertools.combinations(this_tuple, r):
                     self.d[k].append(i)
    
    
    
    
        
    
    def compute_stride(self, join_tree): 
        
        '''procedure to compute the stride of all 
        the variables in a given potential 
        '''
        
        self.strides = defaultdict(int) 
        
        # we compute the strides of each variable 
        # for the first we need to find the first variable to change value 
        
        first_entry = self.entries[0].entries
        second_entry = self.entries[1].entries
            
        for k in first_entry:
            if first_entry[k] != second_entry[k]: 
                break 
        
        
        # now k is the key that change first 
        var = join_tree.get_bbn_node(k).variable
        card = var.card 
        
        # set stride 
        self.strides[k] = 1
        
        # now we need to compute the stride of all the remaining variables 
        previous_entry = first_entry
        
        next_pos = self.strides[k] * card 
        
        
        while len(self.strides) < len(first_entry): 
            
            next_entry = self.entries[next_pos].entries
     
        
            # now find out which entry has changed 
            for k in next_entry:
                if next_entry[k] != previous_entry[k] and k not in self.strides: 
                    break 
                    
            # now k contains the variable which has changed             
            self.strides[k] = next_pos
                    
            # recompute next pos        
            var = join_tree.get_bbn_node(k).variable
            card = var.card 
    
            
            previous_entry = next_entry
            next_pos = self.strides[k] * card 
            
            
        
            
            

class PotentialEntry(object):
    """
    Potential entry.
    """

    def __init__(self):
        """
        Ctor.
        """
        self.entries = dict()
        self.value = 1.0

    def add(self, k, v):
        """
        Adds a node id and its value.

        :param k: Node id.
        :param v: Value.
        :return: This potential entry.
        """
        self.entries[k] = v
        return self

    def matches(self, that):
        """
        Checks if this potential entry matches the specified one. A match is determined with all the keys
        and their associated values in the potential entry passed in matches this one.

        :param that: PotentialEntry.
        :return:
        """
        for k, v in that.entries.items():
            if k not in self.entries or v != self.entries[k]:
                return False
        return True

    def duplicate(self):
        """
        Duplicates this entry.

        :return: PotentialEntry.
        """
        entry = PotentialEntry()
        for k, v in self.entries.items():
            entry.add(k, v)
        entry.value = self.value
        return entry

    def get_entry_keys(self):
        """
        Gets entry keys sorted.

        :return: List of tuples. First tuple is id of variable and second tuple is value of variable.
        """
        return sorted([(k, v) for k, v in self.entries.items()], key=lambda tup: tup[0])

    def get_kv(self):
        """
        Gets key-value pair that may be used for storage in dictionary.

        :return: Key-value pair.
        """
        return '|'.join(list(map(lambda tup: '{}={}'.format(tup[0], tup[1]), self.get_entry_keys()))), self.value

    def __str__(self):
        arr = [(k, v) for k, v in self.entries.items()]
        arr = sorted(arr, key=lambda tup: tup[0])
        arr = ['{}={}'.format(tup[0], tup[1]) for tup in arr]
        s = str.join(',', arr)
        return '{}|{:.5f}'.format(s, self.value)


class PotentialUtil(object):
    
    """
    Potential util.
    """

    @staticmethod
    def pass_single_message(join_tree, x, s, y):
        """
        Single message pass from x -- s -- y (from x to s to y).

        :param join_tree: Join tree.
        :param x: Clique.
        :param s: Separation-set.
        :param y: Clique.
        """
        
        old_sep_set_potential = join_tree.potentials[s.id]
        new_sep_set_potential = PotentialUtil.marginalize_for(join_tree, x, s.nodes)
        join_tree.potentials[s.id] = new_sep_set_potential
        y_potential = join_tree.potentials[y.id]

        ratio = PotentialUtil.divide(new_sep_set_potential, old_sep_set_potential)
        PotentialUtil.multiply(y_potential, ratio)
        
        
    def mulpot(join_tree, pota, potb):

        '''implement multiplication of potentials 
        @params 
        pota, potb : potential
        @returns 
        potprod : potential 
        '''
    
        potprod = Potential() 
        
        if isinstance(potb, (int, float)):
            potprod = pota 
            for i in range(len(potprod.entries)): 
                potprod.entries[i].values *= potb
        
    
        elif isinstance(pota, (int, float)):
            potprod = potb 
            for i in range(len(potprod.entries)): 
                potprod.entries[i].values *= pota 
    
        else:
            
            
            common_vars = list( set(pota.entries[0].entries).intersection(set(potb.entries[0].entries)) )
                    
            nodes = [join_tree.get_bbn_node(x).variable.values for x in common_vars]
    
            for vals in itertools.product(*nodes): 
                
                lst_pota = [] 
                lst_potb = []
    
                intersection_dict = {k:v for k,v in zip(common_vars , vals)}
                
                
                for i in range(len( pota.entries ) ): 
                    
                    if intersection_dict.items() <= pota.entries[i].entries.items() : 
                        lst_pota.append(pota.entries[i])
    
                
                for j in range(len( potb.entries ) ):
                    
                    if intersection_dict.items() <= potb.entries[j].entries.items(): 
                        lst_potb.append(potb.entries[j])
                
                for entry_a in lst_pota:
    
                    for entry_b in lst_potb: 
    
                        # add new entry 
                        e = PotentialEntry() 
                        e.entries = {**entry_a.entries, **entry_b.entries}
                        e.value = entry_a.value * entry_b.value
                        potprod.add_entry(e) 
    
        
        return potprod

    @staticmethod
    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
    
    def divpot(join_tree, pota, potb):
        
        '''implement division of potentials 
        @params 
        pota, potb : potential
        @returns 
        potdiv : potential 
        '''
    
        potdiv = Potential()
    
        if isinstance(potb, (int, float)) and not PotentialUtil.isclose(potb,0):
            potdiv = pota 
            for i in range(len(potdiv.entries)): 
                potdiv.entries[i].values /= potb
    
    
        elif isinstance(pota, (int, float)) and not PotentialUtil.isclose(pota,0):
            potdiv = potb 
            for i in range(len(potdiv.entries)): 
                potdiv.entries[i].values /= pota 
    
        else:
            
            common_vars = list( set(pota.entries[0].entries).intersection(set(potb.entries[0].entries)) )
                    
            nodes = [join_tree.get_bbn_node(x).variable.values for x in common_vars]
    
            for vals in itertools.product(*nodes): 
                
                lst_pota = [] 
                lst_potb = []
    
                intersection_dict = {k:v for k,v in zip(common_vars , vals)}
                
                for i in range(len( pota.entries ) ): 
                    
                    if intersection_dict.items() <= pota.entries[i].entries.items() : 
                        lst_pota.append(pota.entries[i])
    
                
                for j in range(len( potb.entries ) ):
                    
                    if intersection_dict.items() <= potb.entries[j].entries.items(): 
                        lst_potb.append(potb.entries[j])
                
                for entry_a in lst_pota:
    
                    for entry_b in lst_potb: 
    
                        # add new entry 
                        e = PotentialEntry() 
                        e.entries = {**entry_a.entries, **entry_b.entries}
                        
                        if PotentialUtil.isclose(entry_b.value , 0):
                            if PotentialUtil.isclose(entry_a.value , 0):
    
                                e.value = 0 
                                potdiv.add_entry(e) 
    
                            else: 
    
                                e.value = float("inf") 
                                potdiv.add_entry(e) 
    
                        else: 
    
                            e.value = entry_a.value / entry_b.value
                            potdiv.add_entry(e) 


        return  potdiv          
    
    
    
    @staticmethod
    def marginalize_for_from_potential(join_tree, potx, nodes):

        """
        Marginalizes the specified clique's potential over the specified nodes.
        :param join_tree: Join tree.
        :param clique_potential: Clique potential 
        :param nodes: List of BBN nodes.
        :return: Potential.
        """
        potential  = Potential() 
        
        collector_sums = defaultdict(float) 
        
        vars_to_keep = [i.id for i in nodes]
    
        
        for i in range(len( potx.entries ) ): 
            
            k = {key: potx.entries[i].entries[key] for key in vars_to_keep}
                    
            collector_sums[tuple(k.items())] += potx.entries[i].value
            
        # now fill the potential 
        for k,v in collector_sums.items(): 
            
            # add new entry 
            e = PotentialEntry() 
            e.entries = dict(k) 
            e.value = v
            potential.add_entry(e) 
            

        return potential
            
    
    #@staticmethod
    #def marginalize_for_from_potential(join_tree, clique_potential, nodes):
          

    #    potential = PotentialUtil.get_potential_from_nodes(nodes)

    #    for entry in potential.entries:
    #        matched_entries = clique_potential.get_matching_entries(entry)
    #        entry.value = sum([matched_entry.value for matched_entry in matched_entries])

     #   return potential
            
    
    @staticmethod   
    def send_message_online(join_tree , x , s , y, query_vars = None): 
        
        '''
        Send a message from a node to its parent in the calibrated junction tree 
        @param: 
        join tree
        x : sending clique 
        s : separator between x and y 
        y:  receiver clique 
        query_vars: query variables s
        
        @return: 
        potential 
        '''
        
        sep_set_pot = join_tree.potentials[s.id] 
        
        if query_vars != None: 
            querynodesvars = [join_tree.get_bbn_node(x) for x in query_vars]
            nodes_to_carry = list( set(s.nodes).union(set(querynodesvars)) )
        else:
            nodes_to_carry = s.nodes
            
        # nodes_to_carry = [join_tree.get_bbn_node(x) for x in nodes_to_carry]
        
        x_marginaled_pot = PotentialUtil.marginalize_for( join_tree, x , nodes_to_carry )
        
        y_pot = join_tree.potentials[y.id]
        
        prod_pot = PotentialUtil.mulpot(join_tree,  x_marginaled_pot ,  y_pot  )
        
        return PotentialUtil.divpot(join_tree, prod_pot , sep_set_pot)
    
    
    @staticmethod   
    def send_message_online_from_potential(join_tree , x , s , y, query_vars = None): 
        
        '''
        Send a message from a node to its parent in the calibrated junction tree 
        @param: 
        join tree
        x : sending clique 
        s : separator between x and y 
        y:  receiver clique 
        query_vars: query variables s
        
        @return: 
        potential 
        '''
        
        snodes = [x[0] for x in s.entries[0].get_entry_keys()]
        if query_vars != None: 
            nodes_to_carry = list( set(snodes).union(set(query_vars)) )
        else:
            nodes_to_carry = snodes
            
        nodes_to_carry = [join_tree.get_bbn_node(x) for x in nodes_to_carry]
        
        x_marginaled_pot = PotentialUtil.marginalize_for_from_potential( join_tree, x , nodes_to_carry )
      
        prod_pot = PotentialUtil.mulpot(join_tree,  x_marginaled_pot ,  y  )
        
        return PotentialUtil.divpot(join_tree, prod_pot , s)
    
    
    @staticmethod   
    def send_message_online_no_sep(join_tree , x ,  y, query_vars = None): 
        
        '''
        Send a message from a node to its parent in the calibrated junction tree 
        In this case the separator potential is not given as input but computed by marginalization 
        (useful e.g. when using shortcut potentials) 
        @param: 
        join tree
        x : sending clique 
        y:  receiver clique 
        query_vars: query variables s
        
        @return: 
        potential 
        '''
        
        
        nodes_inters = set(x.nodes).interesection(set(y.nodes))
        
        # marginalize from the smallest potential 
        if len(join_tree.potentials[x.id].entries) < len(join_tree.potentials[y.id].entries): 
            sep_set_pot = PotentialUtil.marginalize_for_from_potential( join_tree, x , nodes_inters )
        else:
            sep_set_pot = PotentialUtil.marginalize_for_from_potential( join_tree, y , nodes_inters )
                    
        
        if query_vars != None: 
            nodes_to_carry = list( set(nodes_inters).union(set(query_vars)) )
        else:
            nodes_to_carry = nodes_inters
        
        x_marginaled_pot = PotentialUtil.marginalize_for_from_potential( join_tree, x , nodes_to_carry )
        
        y_pot = join_tree.potentials[y.id]
        
        prod_pot = PotentialUtil.mulpot(join_tree,  x_marginaled_pot ,  y_pot  )
        
        return PotentialUtil.divpot(join_tree, prod_pot , sep_set_pot)
    
    
    
    @staticmethod   
    def send_message_online_no_sep_from_potential(join_tree , potx ,  poty, query_vars = None): 
        
        '''
        Send a message from a node to its parent in the calibrated junction tree 
        In this case the separator potential is not given as input but computed by marginalization 
        (useful e.g. when using shortcut potentials) 
        @param: 
        join tree
        x : sending clique potential 
        y:  receiver clique potential 
        query_vars: query variables 
        
        @return: 
        potential 
        '''
        
        
        scope_potx = [x[0] for x in potx.entries[0].get_entry_keys()]
        
        scope_poty = [x[0] for x in poty.entries[0].get_entry_keys()]
        
        nodes_inters = set(scope_potx).intersection(set(scope_poty))
        
        bbn_nodes_inters = [join_tree.get_bbn_node(x) for x in nodes_inters]
        
        # marginalize from the smallest potential 
        if len(potx.entries) < len(poty.entries): 
            sep_set_pot = PotentialUtil.marginalize_for_from_potential( join_tree, potx , bbn_nodes_inters )
        else:
            sep_set_pot = PotentialUtil.marginalize_for_from_potential( join_tree, poty , bbn_nodes_inters )
                    
        
        if query_vars != None: 
            nodes_to_carry = list( set(nodes_inters).union(set(query_vars)) )
        else:
            nodes_to_carry = nodes_inters
        
        nodes_to_carry = [join_tree.get_bbn_node(x) for x in nodes_to_carry]
        
        x_marginaled_pot = PotentialUtil.marginalize_for_from_potential( join_tree, potx , nodes_to_carry )
        
        prod_pot = PotentialUtil.mulpot(join_tree,  x_marginaled_pot ,  poty  )
        
        return PotentialUtil.divpot(join_tree, prod_pot , sep_set_pot)
    
    
    @staticmethod   
    def send_message_online_no_marginalization(join_tree , x , s , y): 
        
        '''
        Send a message from a node to its parent in the calibrated junction tree 
        without marginalization 
        can be used to compute the overall joint as well (as product of clique potentials to separator potentials)
        @param: 
        join tree
        x : sending clique 
        s : separator between x and y 
        y:  receiver clique 
        query_vars: query variables s
        
        @return: 
        potential 
        '''
        
        sep_set_pot = join_tree.potentials[s.id] 
        
        x_pot = join_tree.potentials[x.id]
        
        y_pot = join_tree.potentials[y.id]
        
        prod_pot = PotentialUtil.mulpot(join_tree,  x_pot ,  y_pot  )
        
        return PotentialUtil.divpot(join_tree, prod_pot , sep_set_pot)
    
    
    #@staticmethod
    #def compute_single_shortcut_potential(jointree, seps): 
        
   
        # i need to get the joint for the entire tree and then marginalize for the 
        # variables that are in the shortcut potentials 
        
   #     allvars = set() 
        
   #     for sep in seps: 
            
            # if they are both different from -1, otherwise we dont need to add anything to @allvars
   #         if sep[0]!= -1 and sep[1]!=-1 : 
            
                # get interesection of the variables in the bottom and top 
   #             u = jointree.get_cliques()[sep[0]].nodes
   #             v = jointree.get_cliques()[sep[1]].nodes
                
                # compute intersection 
   #             uv = set(u).intersection(set(v))
                
                # update set of variable 
   #             allvars.update(uv)
                
        # compute potential 
        #pot = PotentialUtil.marginalize_for_from_potential( jointree, overall_joint, list( allvars ) )
                
    #    steiner_nodes , steiner_edges, steiner_parents , steiner_children , steiner_pivot =  self.extract_steiner_tree( allvars , 0 )
    
    #    sp.potential = self.query(steiner_nodes, steiner_parents, steiner_children, allvars)
        
    
    #   return pot 


    
    '''
    def send_message_online_no_early_marginalization(self, join_tree , x , s , y) : 
        
        Send a single message to answer out-of-clique queries
        
        :param join_tree: Join tree.
        :param x: Clique sender .
        :param s: Separation-set.
        :param y: Clique receiver .
        
        
        sep_set_potential = join_tree.potentials[s.id]
        x_potential = join_tree.potentials[x.id]
        y_potential = join_tree.potentials[y.id]
        product = PotentialUtil.multiply_online( x_potential ,  y_potential )
        
        return  PotentialUtil.divide(product, sep_set_potential)
        
    
    @staticmethod
    def send_message_online(self, join_tree , x , s , y): 
        
        Send a single message to answer out-of-clique queries (with early marginalization)
        
        :param join_tree: Join tree.
        :param x: Clique sender .
        :param s: Separation-set.
        :param y: Clique receiver .
        
        
        sep_set_potential = join_tree.potentials[s.id]
        x_marg_potential = PotentialUtil.marginalize_for(join_tree, x, s.nodes)
        y_potential = join_tree.potentials[y.id]
        product = PotentialUtil.multiply_online( x_marg_potential ,  y_potential )
        
        return  PotentialUtil.divide(product, sep_set_potential)
    
    
    @staticmethod
    def multiply_online(bigger, smaller):
        
        """
        Multiplies two potentials. Order matters.
        :param bigger: Bigger potential.
        :param smaller: Smaller potential.
        """
        
        # adjust so that the one with the larger amount of entries we swap 
        # swap if needed 
        if len(smaller.entries) > len(bigger.entries): 
            smaller , bigger = bigger , smaller 
            
        for entry in smaller.entries:
            for e in bigger.get_matching_entries(entry):
                d = e.value * entry.value
                e.value = d
                
        # bigger should store the value of the product s
        return bigger    
        
        @staticmethod
    def marginalize_for_from_potential(join_tree, clique_potential, nodes):
        """
        Marginalizes the specified clique's potential over the specified nodes.

        :param join_tree: Join tree.
        :param clique_potential: Clique potential 
        :param nodes: List of BBN nodes.
        :return: Potential.
        """
        potential = PotentialUtil.get_potential_from_nodes(nodes)
    #     clique_potential = join_tree.potentials.get(clique.id)

        for entry in potential.entries:
            matched_entries = clique_potential.get_matching_entries(entry)
            entry.value = sum([matched_entry.value for matched_entry in matched_entries])

        return potential
    '''
    
    @staticmethod
    def marginalize_for(join_tree, clique, nodes):

        potential  = Potential() 
        
        collector_sums = defaultdict(float) 
        
        vars_to_keep = [i.id for i in nodes]
    
        potx = join_tree.potentials.get(clique.id) 
        
        for i in range(len( potx.entries ) ): 
            
            k = {key: potx.entries[i].entries[key] for key in vars_to_keep}
                    
            collector_sums[tuple(k.items())] += potx.entries[i].value
            
        # now fill the potential 
        for k,v in collector_sums.items(): 
            
            # add new entry 
            e = PotentialEntry() 
            e.entries = dict(k) 
            e.value = v
            potential.add_entry(e) 
        

        return potential
    
    
    

    @staticmethod
    def normalize(potential):
        """
        Normalizes the potential (make sure they sum to 1.0).

        :param potential: Potential.
        :return: Potential.
        """
        total = sum([entry.value for entry in potential.entries])

        if total != 0.0:
            for entry in potential.entries:
                d = entry.value / total
                entry.value = d

        return potential

    @staticmethod
    def divide(numerator, denominator):
        """
        Divides two potentials.

        :param numerator: Potential.
        :param denominator: Potential.
        :return: Potential.
        """
        potential = Potential()
        for i, entry in enumerate(numerator.entries):
            e = denominator.entries[i]
            d = 0.0 \
                if PotentialUtil.is_zero(entry.value) or PotentialUtil.is_zero(e.value) \
                else (entry.value / e.value)
            new_entry = entry.duplicate()
            new_entry.value = d
            potential.add_entry(new_entry)
        return potential

    @staticmethod
    def is_zero(d):
        """
        Checks if the specified value is 0.0.

        :param d: Value.
        :return: A boolean indicating if the value is zero.
        """
        return 0.0 == d

    @staticmethod
    def multiply(bigger, smaller):
        """
        Multiplies two potentials. Order matters.

        :param bigger: Bigger potential.
        :param smaller: Smaller potential.
        """
        for entry in smaller.entries:
            for e in bigger.get_matching_entries(entry):
                d = e.value * entry.value
                e.value = d
                
                
    @staticmethod
    def multiply_and_return(bigger, smaller):
        """
        Multiplies two potentials. Order matters.

        :param bigger: Bigger potential.
        :param smaller: Smaller potential.
        """
        for entry in smaller.entries:
            for e in bigger.get_matching_entries(entry):
                d = e.value * entry.value
                e.value = d
                
        return bigger 

    @staticmethod
    def get_potential(node, parents):
        """
        Gets the potential associated with the specified node and its parents.

        :param node: BBN node.
        :param parents: Parents of the BBN node (that themselves are also BBN nodes).
        :return: Potential.
        """
        potential = PotentialUtil.get_potential_from_nodes(PotentialUtil.merge(node, parents))
        for i in range(len(potential.entries)):
            prob = node.probs[i]
            potential.entries[i].value = prob
        return potential

    @staticmethod
    def get_potential_from_nodes(nodes):
        """
        Gets a potential from a list of BBN nodes.

        :param nodes: Array of BBN nodes.
        :return: Potential.
        """
        lists = [node.variable.values for node in nodes]
        cartesian = PotentialUtil.get_cartesian_product(lists)

        potential = Potential()
        
        for values in cartesian:
                        
            entry = PotentialEntry()
            for i in range(len(nodes)):
                
                entry.add(nodes[i].id, values[i])
            potential.add_entry(entry)
            
        return potential

    @staticmethod
    def get_cartesian_product(lists):
        """
        Gets the cartesian product of a list of lists of values. For example, if the list is

        * [ ['on', 'off'], ['on', 'off'] ]

        then the result will be a list of the following

        * [ 'on', 'on']
        * [ 'on', 'off' ]
        * [ 'off', 'on' ]
        * [ 'off', 'off' ]

        :param lists: List of list of values.
        :return: Cartesian product of values.
        """
        return [i for i in itertools.product(*lists)]

    @staticmethod
    def merge(node, parents):
        """
        Merges the nodes into one array.

        :param node: BBN node.
        :param parents: BBN parent nodes.
        :return: Array of BBN nodes.
        """
        nodes = [parent for parent in parents]
        nodes.append(node)
        return nodes
