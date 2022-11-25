# VirtualLeadOpt
See demo.ipynb as an example. Formal documentation is under construction.

Parameters to initialize LeadOptimizer:
* lead_candidate: A SMILES string of the lead candidate
* receptor: A pdb or pdbqt file that contains the structure of the receptor protein
* dock_config: A txt file that contains docking configuration
* opt_target: A set of constraints that defines the target of the optimization
* iter_num: Number of iterations
* beam_width: Beam width in the searching algorithm
* exhaustiveness: Defines the exhaustiveness of the searching algorithm (between 0 and 1)
* epsilon: Defines the preference to give retry chances to the suboptimal solutions
* pass_line: Deprecated
* checkpoint: Checkpoint filename
* mmpdb: The mmpdb file used by mmpdb transform. See more on: https://github.com/rdkit/mmpdb
* max_variable_size: Argument used in mmpdb transform
* prediction_workers: Number of threads making prediction on molecular properties
* dock_method: Docking method
* dock_cpu: Number of cpu cores used in docking procedure.
