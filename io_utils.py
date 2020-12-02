# imports 
import argparse

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import warnings
from collections import defaultdict

from generator.bbngenerator import generate_singly_bbn, generate_multi_bbn, convert_for_exact_inference
from generator.bbngenerator import convert_for_drawing
from pptc.inferencecontroller import InferenceController


from graph.dag import Bbn
from graph.edge import Edge, EdgeType
from graph.jointree import EvidenceBuilder
from graph.node import BbnNode
from graph.variable import Variable
from pptc.inferencecontroller import InferenceController


def read_file_bif_format(filepath): 
    
    ''' read file in bif format and create bayesian network
    and build junction tree
    @ params: 
    filepath : path to input file in bif format 
    @ returns 
    join_tree: junction tree 
    
    '''
    
    f = open(filepath, "r")
    
    var_to_values = dict() 
    probabilities = defaultdict(list)
    dependencies = dict()
    allnames = []
    
    # iterate all rows of the file 
    with open(filepath, 'r') as f:
        
        lines = f.readlines()
        
        i = 0 
        
        while i < len(lines): 
            
            # next line 
            line = lines[i]
            
            tokens = line.split(" ") 
            
            if tokens[0]=='variable':
            
                # new variable
                name = tokens[1] 
                allnames.append(name)
                
                # next line contains values 
                i+=1
                line = lines[i] 
                
                # get all the values assuming they are always expressed in the same way 
                this_var_values = line.split("{")[-1][1:-3].split(", ")
                
                var_to_values[name] = this_var_values
                
                # next line is closed curly brackets - we skip it  
                i+=2
                    
                
            elif tokens[0]=='probability':
                
                name = tokens[2] 
                
                dependencies[name] = [token.replace(',', '') for token in tokens[4:-2]] 
                            
                
                #  get probabilities 
                i+=1
                line = lines[i] 
                                
                while line!="}\n":
                    
                    tokens = line.split(" ")
                 
            
                                        
                    for token in tokens:   
                        
                        if sum([char.isdigit() for char in token[:-1]]) > 0 and "." in token:
                        #and not any(c.isalpha() for c in token[:-1]) : 
                            try: 
                                probabilities[name].append(float(token[:-1]))
                            except: 
                                probabilities[name].append(float(token[:-2]))
                                
                    i+=1
                    line = lines[i]
                            
                # finished the while loop we go to the next          
                i+=1     
            
            
            else: 
                # only next line (extra lines e.g. first)
                i+=1 
            
            
        # now create bbn 
        bbn = Bbn() 
        
        name_to_node = dict() 
       
        #add nodes 
        for i in range(len(allnames)): 
            
            name = allnames[i] 
            
            node = BbnNode(Variable(i, name, var_to_values[name]), probabilities[name])
            
            bbn.add_node(node) 
            
            # for edges 
            name_to_node[name] = node 
            
        #add edges 
        for i in range(len(allnames)): 
            
            name = allnames[i] 
            
            this_dependecies = dependencies[name] 
            
            for dependency in this_dependecies: 
            
                bbn.add_edge(Edge(name_to_node[dependency], name_to_node[name], EdgeType.DIRECTED))
        
            
        # convert bbn to junction tree 
        join_tree = InferenceController.apply(bbn)
                
        
        return join_tree



def generate_random_bbn(n): 
    
    ''' generate random bayesian network and return junction tree 
    @params: 
    n: number of nodes in the bbn network '''
    
    g, p = generate_singly_bbn(parser.n , max_iter=10)
    s_bbn = convert_for_exact_inference(g, p)
    bbn = convert_for_drawing(s_bbn)
    
    join_tree = InferenceController.apply(bbn)
        
    return join_tree
    
    

    


    