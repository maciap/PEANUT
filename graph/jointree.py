from collections import defaultdict
from copy import deepcopy
from enum import Enum
# from collections import defaultdict
# import queue 
from graph.edge import JtEdge
from graph.graph import Ug
from graph.node import SepSet, Clique, BbnNode
from graph.potential import Potential, PotentialEntry, PotentialUtil
from graph.variable import Variable



class JoinTree(Ug):
    """
    Join tree.
    """

    def __init__(self):
        # self.all_cliques = self.get_cliques() 
        # self.adjacency_list , self.map_nodes, self.map_nodes_inv = self.convert_to_adjancency_list()
        """
        Ctor.
        """
        Ug.__init__(self)
        self.potentials = dict()
        self.evidences = dict()
        self.listener = None
        self.parent_info = defaultdict(set)
        self.pivot = 0 # can change later 
        # self.adjacency_list , self.map_nodes, self.map_nodes_inv = self.convert_to_adjancency_list()
        # self.all_cliques = self.get_cliques() 
        # self.__all_nodes__ = None

    def __deepcopy__(self, memodict={}):
        nodes = deepcopy(self.nodes, memodict)
        edges = deepcopy(self.edges, memodict)
        edge_map = deepcopy(self.edge_map, memodict)
        neighbors = deepcopy(self.neighbors, memodict)
        potentials = deepcopy(self.potentials, memodict)
        evidences = deepcopy(self.evidences, memodict)
        parent_info = deepcopy(self.parent_info, memodict)

        jt = JoinTree()
        jt.nodes = nodes
        jt.edges = edges
        jt.edge_map = edge_map
        jt.neighbors = neighbors
        jt.potentials = potentials
        jt.evidences = evidences
        jt.parent_info = parent_info
        return jt

    def get_posteriors(self):
        """
        Gets the posterior for all nodes.

        :return: Map. Keys are node names; values are map of node values to posterior probabilities.
        """
        bbn_nodes = self.get_bbn_nodes()

        posteriors = {}

        for bbn_node in bbn_nodes:
            potential = self.get_bbn_potential(bbn_node)

            m = {}
            for potential_entry in potential.entries:
                k = ''.join([f'{y}' for _, y in potential_entry.entries.items()])
                m[k] = potential_entry.value

            name = bbn_node.variable.name
            posteriors[name] = m

        return posteriors
    
    
    def convert_to_adjancency_list(self): 
            
        '''
        Convert junction tree to standard adjacency list for easier manipulation 
        Rerturn dict: integer id to adjancemcy list 
        and dicts mapping simple id to string "1-2-3" and vice versa 
        '''

        map_nodes = dict() 
        adj_lst = defaultdict(set) 
    
        for i in range(len( self.get_cliques() )):     
            map_nodes[i] = self.get_cliques()[i].id
        map_nodes_inv = {v:k for k,v in map_nodes.items()}
    
        
        for sepset in self.get_sep_sets(): 
            
            pair = []
            
            for i in range(len( self.get_cliques() )):   
                
                l = len(  self.get_cliques()[i].id.split("-")  )
                
                #for j in range(i+1, len( join_tree.get_cliques() )): 
                if sepset.id.split("-")[:l] == self.get_cliques()[i].id.split("-") or sepset.id.split("-")[-l:] == self.get_cliques()[i].id.split("-"):
                    pair.append(i) 
                    
                    if len(pair) == 2:
                        adj_lst[pair[0]].add(pair[1]) 
                        adj_lst[pair[1]].add(pair[0]) 
                        # next sepset 
                        break 
                        
            
        # convert to list         
        adj_lst_lst = defaultdict(list)
    
        for k , v in adj_lst.items(): 
            adj_lst_lst[k] = list(v)

        return adj_lst_lst , map_nodes , map_nodes_inv    


    # function to determine level of  
    # each node starting from x using BFS  
    def get_parents_and_children(self, graph, x): 
        
        '''
        performs a simple bfs traversal of the tree 
        
        params: 
        @graph: adjacency list 
        @x: pivot (root) 
        
        returns:
        level: height level in the tree 
        parent: dict giving the parent of key 
        children: dictionary giving the children of the key 
        
        '''
          
        # array to store level of each node  
        level = dict()   
        parent = dict()
        children = defaultdict(list) 
        
      
        # create a queue  
        q = list()
        
        #set of visited nodes 
        visited = set() 
      
        # enqueue element x  
        q.append(x) 
      
        # initialize level of source  
        # node to 0  
        level[x] = 0
      
        # marked it as visited  
        visited.add(x)
        
        # do until queue is empty  
        while len(q) > 0:  
      
            # get the first element of queue  
            x = q[0] 
            del q[0] 
      
            # traverse neighbors of node x 
            for b in range(len(graph[x])): 
                  
                # b is neighbor of node x  
      
                # if b is not marked already  
                if b not in visited:  
      
                    # enqueue b in queue  
                    q.append(b)  
      
                    # level of b is level of x + 1  
                    level[b] = level[x] + 1
                    parent[b] = x
                    children[x].append(b) 
      
                    # mark b  
                    visited.add(b) 
      
        return level  , parent , children
    
    
    
    def compute_heights_parents(self, pivot = 0): 
        
        tr = self.convert_to_adjancency_list()[0]
        heights , parent = self.Levels(tr, len(sorted(list(tr.keys()))), pivot) 
        
        return heights, parent
    
    
    def extract_steiner_tree(self , querynodes, pivot): 
        
        '''
        extract Steiner tree 
        
        param queryvar : set of query variables 
        Returns set of steiner ndoes 
        '''
        # heights, parents = self.compute_heights_parents(pivot)        
        tr, map_nodes = self.convert_to_adjancency_list() # this we have to change all this should be accessible from the self.jointree object 
        
        V = len( list(map_nodes.keys()) )
        
        heights, parents = self.Levels(tr, V, pivot) 
        
        steiner_nodes = set() 
        
        # start with query variables in the steiner tree 
        steiner_nodes.update(querynodes)
        
        # find lowest common ancestor using heights and paths to it 
        heights_queryvars = { k:v for k,v in  heights.items() if k in  querynodes }
        
        # find lowest common ancestor 
        # the stenier tree is found when all the nodes in the dictionary
        # terminates when we have deleted all the nodes but one 
        
        while len(heights_queryvars) > 1: 
            
            # find minimum height 
            key_min = min( heights_queryvars.keys(), key=(lambda k: heights_queryvars[k]) )  
    
            steiner_nodes.add( parents[key_min] )
            
            # delete current node with minumum height 
            del heights_queryvars[key_min]
            
            # add parent (in the end all the nodes will have a common parent)
            heights_queryvars[parents[key_min]] = heights[parents[key_min]]
            
        return  steiner_nodes


    def get_bbn_potential(self, node):
        """
        Gets the potential associated with the specified BBN node.

        :param node: BBN node.
        :return: Potential.
        """
        clique = node.metadata['parent.clique']
        potential = PotentialUtil.normalize(PotentialUtil.marginalize_for(self, clique, [node]))
        return potential

    def unmark_cliques(self):
        """
        Unmarks the cliques.
        """
        for clique in self.get_cliques():
            clique.unmark()
            

    def update_bbn_cpts(self, cpts):
        """
        Updates the CPTs of the BBN nodes.

        :param cpts: Dictionary of CPTs. Keys are ids of BBN node and values are new CPTs.
        :return: None
        """
        bbn_nodes = {node.id: node for clique in self.get_cliques() for node in clique.nodes}
        for idx, cpt in cpts.items():
            if idx in bbn_nodes:
                bbn_nodes[idx].probs = cpt
                bbn_nodes[idx].potential = None

    def get_bbn_node_and_parents(self):
        """
        Gets a map of nodes and its parents.

        :return: Map. Keys are node ID and values are list of nodes.
        """
        bbn_nodes = {node.id: node for clique in self.get_cliques() for node in clique.nodes}
        result = {node: [pa for pa_id, pa in bbn_nodes.items() if pa_id in self.parent_info[node_id]]
                  for node_id, node in bbn_nodes.items()}
        return result

    def __get_bbn_nodes__(self):
        """
        Gets all BBN nodes (cached).

        :return: Dictionary of BBN nodes.
        """
        # if self.__all_nodes__ is None:
        #     self.__all_nodes__ = {node.id: node for clique in self.get_cliques() for node in clique.nodes}
        # return self.__all_nodes__
        result = {node.id: node for clique in self.get_cliques() for node in clique.nodes}
        return result

    def get_bbn_nodes(self):
        """
        Gets all the BBN nodes in this junction tree.

        :return: List of BBN nodes.
        """
        return list(self.__get_bbn_nodes__().values())

    def get_bbn_node(self, id):
        """
        Gets the BBN node associated with the specified id.

        :param id: Node id.
        :return: BBN node or None if no such node exists.
        """
        bbn_nodes = self.__get_bbn_nodes__()
        if id in bbn_nodes:
            return bbn_nodes[id]
        return None

    def get_bbn_node_by_name(self, name):
        """
        Gets the BBN node associated with the specified name.

        :param name: Node name.
        :return: BBN node or None if no such node exists.
        """
        bbn_nodes = {node.variable.name: node for clique in self.get_cliques() for node in clique.nodes}
        if name in bbn_nodes:
            return bbn_nodes[name]
        return None

    def find_cliques_with_node_and_parents(self, id):
        """
        Finds all cliques in this junction tree having the specified node and its parents.

        :param id: Node id.
        :return: Array of cliques.
        """
        ids = self.__get_parent_ids__(id)
        ids.append(id)

        set1 = set(ids)
        result = [clique for clique in self.get_cliques() if clique.get_node_ids().issuperset(set1)]
        return result

    def add_potential(self, clique, potential):
        """
        Adds a potential associated with the specified clique.

        :param clique: Clique.
        :param potential: Potential.
        :return: This join tree.
        """
        self.potentials[clique.id] = potential
        return self
    
    def add_potential_sepset(self, sepset, potential):
        """
        Adds a potential associated with the specified separator .
        For all the subtrees which have intersection = 1 we need to manually add each 
        variable unfortunately 

        :param clique: Separator.
        :param potential: Potential.
        :return: This join tree.
        """
        self.potentials[sepset.id] = potential
        return self
    
    
    def add_potential_single_var_sepset(self):
        """
        Adds all the potentials associated with separator potential
        :return: This join tree.
        """
        
        list_cliques = list(self.get_cliques()) 
        for i in range(len(list_cliques)-1):
            for j in range(i+1, len(list_cliques)): 
                
                clique1 = list_cliques[i]
                clique2 = list_cliques[j]
                bool_int, scope1, scope2, int_scope = clique1.intersects(clique2)
                
                if bool_int and len(int_scope) == 1: 
                    
                    sofint = clique1.get_sep_set(clique2)
                    name = sofint.get_sid()
                    nod = self.get_bbn_node_by_name(name)
                    potential = self.get_bbn_potential(nod)
                    self.potentials[sofint.id] = potential
        
        return self
    
    
    def get_cliques(self):
        """
        Gets all the cliques in this junction tree.

        :return: Array of cliques.
        """
        return [clique for clique in self.get_nodes() if not isinstance(clique, SepSet)]

    def get_sep_sets(self):
        """
        Gets all the separation sets in this junction tree.

        :return: Array of separation sets.
        """
        return [sep_set for sep_set in self.get_nodes() if isinstance(sep_set, SepSet)]

    def add_edge(self, edge):
        """
        Adds an JtEdge.

        :param edge: JtEdge.
        :return: This join tree.
        """
        
        
        if not isinstance(edge, JtEdge):
            return self

        sep_set = edge.sep_set
        lhs = edge.i
        rhs = edge.j

        if self.__shouldadd__(edge):
            self.add_node(sep_set)
            self.add_node(lhs)
            self.add_node(rhs)

            self.edge_map[lhs.id].add(sep_set.id)
            self.edge_map[rhs.id].add(sep_set.id)
            self.neighbors[lhs.id].add(sep_set.id)
            self.neighbors[rhs.id].add(sep_set.id)

            self.edge_map[sep_set.id].add(lhs.id)
            self.edge_map[sep_set.id].add(rhs.id)
            self.neighbors[sep_set.id].add(lhs.id)
            self.neighbors[sep_set.id].add(rhs.id)

            self.edges[edge.key] = edge

        return self
    
    
    def compute_overall_joint(self):    

        ''' compute the overall joint probability distribution associated with the junction tree 
        It can be obtained as the ratio of clique potentials to separator potentials 
        This it is then 
        used to compute shortcut
        potentials through marginalization '''
        
        all_cliques = self.get_cliques() 
        all_sepsets = self.get_sep_sets() 
        
        #initialize numerator with first clique potential 
        prod_pot_num = self.potentials[all_cliques[0].id] 
        for cli in all_cliques[1:]:
            prod_pot_num = PotentialUtil.mulpot( prod_pot_num , self.potentials[cli.id] ) 
            
        #initialize denominator with first clique potential 
        prod_pot_den = self.potentials[ all_sepsets[0].id ] 
        for s in all_sepsets[1:]:
            prod_pot_den = PotentialUtil.mulpot( prod_pot_den  , self.potentials[s.id] ) 


        # finally divide 
        return PotentialUtil.divpot(prod_pot_num , prod_pot_den)
    
    

    def get_flattened_edges(self):
        """
        Gets all the edges "flattened" out. Since separation-sets are really hyper-edges, this method breaks
        separation-sets into two edges.

        :return: Array of edges.
        """
        edges = []
        for edge in self.edges.values():
            edges.append(edge.get_lhs_edge())
            edges.append(edge.get_rhs_edge())
        return edges

    def set_listener(self, listener):
        """
        Sets the listener.

        :param listener: JoinTreeListener.
        """
        self.listener = listener

    def get_evidence(self, node, value):
        """
        Gets the evidence associated with the specified BBN node and value.

        :param node: BBN node.
        :param value: Value.
        :return: Potential (the evidence).
        """
        if node.id not in self.evidences:
            self.evidences[node.id] = dict()

        if value not in self.evidences[node.id]:
            entry = PotentialEntry()
            entry.add(node.id, value)
            entry.value = 1.0

            potential = Potential()
            potential.add_entry(entry)

            self.evidences[node.id][value] = potential

        result = self.evidences[node.id][value]
        return result

    def get_change_type(self, evidences):
        """
        Gets the change type associated with the specified list of evidences.

        :param evidences: List of evidences.
        :return: ChangeType.
        """
        changes = []
        for evidence in evidences:
            node = evidence.node
            potentials = self.evidences[node.id]
            change = evidence.compare(potentials)
            changes.append(change)

        count = len([change_type for change_type in changes if ChangeType.RETRACTION == change_type])
        if count > 0:
            return ChangeType.RETRACTION

        count = len([change_type for change_type in changes if ChangeType.UPDATE == change_type])
        if count > 0:
            return ChangeType.UPDATE

        return ChangeType.NONE

    def get_unobserved_evidence(self, node):
        """
        Gets the unobserved evidences associated with the specified node.

        :param node: BBN node.
        :return: Evidence.
        """
        evidence = Evidence(node, EvidenceType.UNOBSERVE)
        for value in node.variable.values:
            evidence.add_value(value, 1.0)
        return evidence

    def unobserve(self, nodes):
        """
        Unobserves a list of nodes.

        :param nodes: List of nodes.
        :return: This join tree.
        """
        evidences = [self.get_unobserved_evidence(node) for node in nodes]
        self.update_evidences(evidences)
        return self

    def unobserve_all(self):
        """
        Unobserves all BBN nodes.

        :return: This join tree.
        """
        self.unobserve(self.get_bbn_nodes())
        return self

    def update_evidences(self, evidences):
        """
        Updates this join tree with the list of specified evidence.

        :param evidences: List of evidences.
        :return: This join tree.
        """
        for evidence in evidences:
            evidence.validate()
        change = self.get_change_type(evidences)
        for evidence in evidences:
            node = evidence.node
            potentials = self.evidences[node.id]

            for k, v in evidence.values.items():
                potential = potentials[k]
                potential.entries[0].value = v
        self.__notify_listener__(change)
        return self

    def set_observation(self, evidence):
        """
        Sets a single observation.

        :param evidence: Evidence.
        :return: This join tree.
        """
        potentials = self.evidences[evidence.node.id]

        pvalues = []
        for v, potential in potentials.items():
            entry = potential.entries[0]
            p = entry.value
            if 1.0 == p:
                pvalues.append(v)

        cvalues = []
        for v, likelihood in evidence.values.items():
            if 1.0 == likelihood:
                cvalues.append(v)

        if 1 == len(pvalues):
            last_value = pvalues[0]
            curr_value = cvalues[0]
            if last_value == curr_value:
                self.unobserve([evidence.node])
            else:
                self.update_evidences([evidence])
        else:
            self.update_evidences([evidence])

        return self

    @staticmethod
    def to_dict(jt):
        """
        Converts a junction tree to a serializable dictionary.

        :param jt: Junction tree.
        :return: Dictionary.
        """

        def nodes_to_dict(nodes):
            d = {}
            for n in nodes:
                if isinstance(n, SepSet):
                    d[n.id] = {
                        'left': n.left.id,
                        'right': n.right.id,
                        'type': 'sepset'
                    }
                elif isinstance(n, Clique):
                    d[n.id] = {
                        'node_ids': list(n.node_ids),
                        'type': 'clique'
                    }
            return d

        def edges_to_dict(edges):
            return [e.sep_set.id for e in edges]

        bbn_nodes = {n.id: n.to_dict() for n in jt.get_bbn_nodes()}
        jt_nodes = nodes_to_dict(jt.get_nodes())
        jt_edges = edges_to_dict(jt.get_edges())

        return {
            'bbn_nodes': bbn_nodes,
            'jt': {
                'nodes': jt_nodes,
                'edges': jt_edges,
                'parent_info': jt.parent_info
            }
        }

    @staticmethod
    def from_dict(d):
        """
        Converts a dictionary to a junction tree.

        :param d: Dictionary.
        :return: Junction tree.
        """

        def get_variable(d):
            return Variable(d['id'], d['name'], d['values'])

        def get_bbn_node(d):
            return BbnNode(get_variable(d['variable']), d['probs'])

        def get_clique(d, bbn_nodes):
            return Clique([bbn_nodes[idx] if idx in bbn_nodes else bbn_nodes[str(idx)] for idx in d['node_ids']])

        def get_sep_set(lhs_clique, rhs_clique):
            _, lhs, rhs, intersection = lhs_clique.intersects(rhs_clique)
            return SepSet(lhs_clique, rhs_clique, lhs, rhs, intersection)

        bbn_nodes = {k: get_bbn_node(n) for k, n in d['bbn_nodes'].items()}

        cliques = [get_clique(clique, bbn_nodes)
                   for k, clique in d['jt']['nodes'].items() if clique['type'] == 'clique']
        cliques = {c.id: c for c in cliques}

        sepsets = [get_sep_set(cliques[s['left']], cliques[s['right']])
                   for k, s in d['jt']['nodes'].items() if s['type'] == 'sepset']
        sepsets = {s.id: s for s in sepsets}

        edges = [JtEdge(sepsets[e]) for e in d['jt']['edges']]

        jt = JoinTree()
        if len(edges) > 0:
            for e in edges:
                jt.add_edge(e)
        else:
            jt.nodes = cliques

        jt.parent_info = {int(k): v for k, v in d['jt']['parent_info'].items()}
        return jt
    
   

    def __shouldadd__(self, edge):
        """
        Checks if the specified edge should be added.

        :param edge: Edge.
        :return: A boolean indicating if the specified edge should be added.
        """
        lhs = edge.i
        rhs = edge.j

        if lhs.id == rhs.id:
            return False

        if not PathDetector(self, lhs.id, rhs.id).exists():
            return True

        return False

    def __get_parent_ids__(self, id):
        """
        Gets the parent ids of the specified node id.

        :param id: Node id.
        :return: Array of parent ids.
        """
        node = self.get_bbn_node(id)
        if 'parents' in node.metadata:
            return [n.id for n in node.metadata['parents']]
        return []

    def __notify_listener__(self, change):
        """
        Notifies the JoinTreeListener, if any.

        :param change: ChangeType.
        """
        if self.listener is None:
            return
        if ChangeType.RETRACTION == change:
            self.listener.evidence_retracted(self)
        elif ChangeType.UPDATE == change:
            self.listener.evidence_updated(self)


class PathDetector(object):
    """
    Detects path between two nodes.
    """

    def __init__(self, graph, start, stop):
        """
        Ctor.

        :param graph: Join tree.
        :param start: Start node id.
        :param stop: Stop node id.
        """
        self.graph = graph
        self.start = start
        self.stop = stop
        self.seen = set()

    def exists(self):
        """
        Checks if a path exists.

        :return: True if a path exists, otherwise, false.
        """
        if self.start == self.stop:
            return True
        else:
            return self.__find__(self.start)

    def __find__(self, i):
        """
        Checks if a path exists from the specified node to the stop node.

        :param i: Node id.
        :return: True if a path exists, otherwise, false.
        """
        neighbors = set()

        try:
            neighbors = self.graph.get_neighbors(i)
        except KeyError:
            pass

        if self.stop in neighbors:
            return True

        self.seen.add(i)
        for neighbor in neighbors:
            if neighbor not in self.seen and self.__find__(neighbor):
                return True

        return False


class JoinTreeListener(object):
    """
    Interface like class used for listening to a join tree.
    """

    def evidence_retracted(self, join_tree):
        """
        Evidence is retracted.

        :param join_tree: Join tree.
        """
        pass

    def evidence_updated(self, join_tree):
        """
        Evidence is updated.

        :param join_tree: Join tree.
        """
        pass


class EvidenceType(Enum):
    """
    Evidence type.
    """
    VIRTUAL = 1
    FINDING = 2
    OBSERVATION = 3
    UNOBSERVE = 4


class ChangeType(Enum):
    """
    Change type.
    """
    NONE = 1
    UPDATE = 2
    RETRACTION = 3


class EvidenceBuilder(object):
    """
    Evidence builder.
    """

    def __init__(self):
        """
        Ctor.
        """
        self.values = dict()
        self.node = None
        self.type = EvidenceType.OBSERVATION

    def with_node(self, node):
        """
        Adds a BBN node.

        :param node: BBN node.
        :return: Builder.
        """
        self.node = node
        return self

    def with_type(self, type):
        """
        Adds the EvidenceType.

        :param type: EvidenceType.
        :return: Builder.
        """
        self.type = type
        return self

    def with_evidence(self, val, likelihood):
        """
        Adds evidence.

        :param val: Value.
        :param likelihood: Likelihood.
        :return: Builder.
        """
        self.values[val] = likelihood
        return self

    def build(self):
        """
        Builds an evidence.

        :return: Evidence.
        """
        evidence = Evidence(self.node, self.type)
        for k, v in self.values.items():
            evidence.add_value(k, v)
        return evidence


class Evidence(object):
    """
    Evidence.
    """

    def __init__(self, node, type):
        """
        Ctor.

        :param node: BBN node.
        :param type: EvidenceType.
        """
        self.node = node
        self.type = type
        self.values = dict()

    def add_value(self, value, likelihood):
        """
        Adds a value.

        :param value: Value.
        :param likelihood: Likelihood.
        :return: This evidence.
        """
        self.values[value] = likelihood
        return self

    def compare(self, potentials):
        """
        Compares this evidence with previous ones.

        :param potentials: Map of potentials.
        :return: The ChangeType from the comparison.
        """
        that = self.__convert__(potentials)

        that_unobserve = self.__is_unobserved__(that)
        this_unobserve = self.__is_unobserved__(self.values)

        if that_unobserve and this_unobserve:
            return ChangeType.NONE

        that_observe = self.__is_observed__(that)
        this_observe = self.__is_observed__(self.values)

        if that_observe and this_observe:
            s1 = self.__get_observed_value__(that)
            s2 = self.__get_observed_value__(self.values)

            if s1 == s2:
                return ChangeType.NONE
            else:
                return ChangeType.RETRACTION

        return ChangeType.RETRACTION

    @staticmethod
    def __convert__(potentials):
        """
        Converts potentials to a map (dict).

        :param potentials: Potentials.
        :return: Dict where keys are BBN node values and values are likelihoods.
        """
        m = dict()
        for k, v in potentials.items():
            m[k] = v.entries[0].value
        return m

    @staticmethod
    def __is_unobserved__(values):
        """
        Checks if the values represent an unobserved evidence. If all likelihoods are 1.0, then this
        map of values represent unobserved evidence.

        :param values: Map of values, where keys are values and values are likelihoods.
        :return: A boolean indicating if the values represent unobserved evidence.
        """
        count = 0
        for k, v in values.items():
            count += v
        return count == len(values)

    @staticmethod
    def __is_observed__(values):
        """
        Checks if the values represent an observed evidence. If all likelihoods are 0 with exactly
        one of them being 1, then this map of values represent observed evidence.

        :param values: Map of values, where keys are values and values are likelihoods.
        :return: A boolean indicating if the values represent observed evidence.
        """
        one = 0
        zero = 0

        for k, v in values.items():
            if 1.0 == v:
                one += 1
            else:
                zero += 1

        return 1 == one and len(values) - 1 == zero

    @staticmethod
    def __get_observed_value__(values):
        """
        Gets the value that is observed (the value whose likelihood is 1.0).

        :param values: Map of values, where keys are values and values are likelihoods.
        :return: Observed value.
        """
        strs = [k for k in values if 1.0 == values[k]]
        return strs[0]

    @staticmethod
    def __normalize_value_between_zero_one__(val):
        """
        Normalizes the specified value to the range [0.0, 1.0]. If the specified value is less than 0.0 then
        0.0 is returned; if the specified value is greater than 1.0 then 1.0 is returned.

        :param val: Value.
        :return: A value in the range [0.0, 1.0].
        """
        if 0.0 <= val <= 1.0:
            return val
        elif val < 0.0:
            return 0.0
        else:
            return 1.0

    @staticmethod
    def __normalize_value_zero_or_one__(val):
        """
        Normalizes the specified value to either 0.0 or 1.0 (and nothing else). If the specified value is anything
        greater than 0.0, then a 1.0 is returned, else a 0.0 is returned.

        :param val: Value.
        :return: 0.0 or 1.0.
        """
        if val > 0.0:
            return 1.0
        else:
            return 0.0

    def validate(self):
        """
        Validates this evidence.

        * virtual evidence: each likelihood must be in the range [0, 1].
        * finding evidence: all likelihoods must be exactly 1.0 or 0.0.
        * observation evidence: exactly one likelihood is 1.0 and all others must be 0.0.
        """
        for value in self.node.variable.values:
            if value not in self.values:
                self.values[value] = 0.0

        if EvidenceType.VIRTUAL == self.type:
            for value in self.node.variable.values:
                self.values[value] = Evidence.__normalize_value_between_zero_one__(self.values[value])
        elif EvidenceType.FINDING == self.type:
            for value in self.node.variable.values:
                self.values[value] = Evidence.__normalize_value_zero_or_one__(self.values[value])

            count = sum([x for x in self.values.values()])
            if 0.0 == count:
                for value in self.node.variable.values:
                    self.values[value] = 1.0
        elif EvidenceType.OBSERVATION == self.type:
            tuples = []
            for k, v in self.values.items():
                pair = (k, v)
                tuples.append(pair)
            tuples = sorted(tuples, key=lambda x: (x[1]), reverse=True)
            key = tuples[0][0]
            for value in self.node.variable.values:
                if key == value:
                    self.values[value] = 1.0
                else:
                    self.values[value] = 0.0
        elif EvidenceType.UNOBSERVE == self.type:
            for value in self.node.variable.values:
                self.values[value] = 1.0
