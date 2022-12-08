

## DeepProLigand
DeepProLigand is a deep learning bioinformatics pipeline developed by BML lab at University of Missouri, Columbia for modeling protein-ligand interaction with cryo-EM data in 2021 Ligand Model Challenge.

## Setup Environment
1. The required packages are found in ``environment.yml`` 

``conda env create -f environment.yml``

2. Once the installation is completed, activate the conda environment

``conda activate deep-pro-ligand``

3. The pipeline makes use of DeepTracer, UCSF Chimera, and PyRosetta. Please run DeepTracer through the webpage and install USCF Chimera and PyRosetta through their websites:

[DeepTracer](https://deeptracer.uw.edu/home), this link will take users to: ```https://deeptracer.uw.edu/home```

[UCSF Chimera Download](https://www.cgl.ucsf.edu/chimera/download.html), this link will take users to: ```https://www.cgl.ucsf.edu/chimera/download.html```

[PyRosetta Download](https://www.pyrosetta.org/downloads), this link will take users to: ```https://www.pyrosetta.org/downloads```

## Run
1. Using the protein density map as an input, please run DeepTracer to generate the atomic backbone structure of the protein.

2. The output of DeepTracer and reference structure (deposited into PDB) for the density map needs to be superimposed using UCSF Chimera's matchmaker function. The matchmaker function can be enabled from UCSF Chimera application as : Tools > Structure Comparision > Matchmaker. Once the structures are superimposed, they need to be saved in a single PDB file. 
Please find the superimposed structures for all three targets in ``'/data/'`` directory of this repository.

3. The next steps requires PyRosetta to be installed into the system. To generate the average protein structure run as shown below for each density map:

       python3 avg_7770.py

``7770`` is the density map's name. For, 30210 density map, run ``python3 avg_30210.py``. Similarly, for 22898 density map, run ``python3 avg_22898.py``.
The outputs are currently saved in 
``'/data/averaged'``, output directory can be change from the program.

4. To generate the refined protein structure run the following:

       python3 fastrelax.py

The outputs of fastrelax program are currently saved in the below directory:
``'/data/refined/'``

