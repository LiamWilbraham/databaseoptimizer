# Database Optimizer

A python implementation of an algorithm for obtaining an optimally diverse subset of molecules from a large database.[1]

## Installing the module
Simply clone the repository and add the path to `databaseoptimizer` to your `.bashrc` or `.bash_profile`.

## Using the module
```python
from databaseoptimizer import DatabaseOptimizer

optimizer = DatabaseOptimizer(smiles_list, desired_library_size)
optimizer.optimize()
```
where `smiles_list` is a list of smiles that form the larger database and `desired_library_size` is the required size of the optimized molecule subset. Once `optimize()` is complete, the optimized molecule subset is written to a file `optimized-library.csv` and is also returned directly as a list `optimizer.optimized_library` (i.e. an attribute of `DatabaseOptimizer`).

## Requirements
* RDKit
* numpy

### RDKit
RDKit can be installed via anaconda (recommended)
`conda install -c rdkit rdkit`

## References
[1] John D. Holliday  Sonia S. Ranade  Peter Willett, **1995**, doi.org/10.1002/qsar.19950140602
