# VASP_Tool
Post-processing of VASP output for electrochemical simulation

For HER&wf.py, ase package is needed and numpy >1.17 is necessary.

Execute wf&HER.py with `python3 HER&wf.py z 1e-3`. direction x,y,z and cutoff vacuum charge density can be specified.

Execute xdat2pos.py with `python3 xdat2pos.py 1 2 ...`. MD step number can be specified, the POSCAR will be generated under sub folder with INCAR POTCAR and KPOINTS files are copied.
