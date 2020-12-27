# PEANUT
WORKLOAD AWARE MATERIALIZATION IN JUNCTION TREE INFERENCE 

Contact author: martino.ciaperoni@aalto.fi 

## Code 

### Cython

To use the code, first run 'python setup.py build_ext --inplace'. command builds the .c files created by Cython. Alternatively, without running the mentioned command, it is possible to directly execute the Python code.

### Execution Example 

**Run**:
main_example.py asia 1000 1.0 <br/>
with **arguments**: <br/>
dataset:asia (Bayesian network - bif format) <br/>
K: 1000 (target budget) <br/>
epsilon: 1.0 approximation level <br/>


After the execution, the outputs folder contain the average cost savings obtained with the exact optimal shortcut potentials for a set of generate queries. 


