# PEANUT
WORKLOAD AWARE MATERIALIZATION IN JUNCTION TREE INFERENCE 

Contact author: martino.ciaperoni@aalto.fi 

## Code 

### Cython

To use the code, first run 'python setup.py build_ext --inplace' from the folder 'span_cores/'. This command builds the .c files created by Cython. Alternatively, without running the mentioned command, it is possible to directly execute the Python code.

### Execution Example 

**Run**:
example.py <br/>
with **arguments**: <br/>
dataset:asia (Bayesian network - bif format) <br/>
K: 1000 (budget) <br/>
epsilon: 1.0 approximation level <br/>

After the execution, the outputs folder contain the average cost savings obtained with the exact optimal shortcut potentialsfor a set of generate queries. 


