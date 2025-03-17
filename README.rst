QCxMS2 Plot
========================

This simple Python script is a tool to visualize the spectrum and fragmentation pathways of QCxMS2. 

### Installation

You can install the project in an existing virtual environment (provided for example by the package managers `conda` or `mamba` (see also [here](https://github.com/conda-forge/miniforge) and [here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html))).
With `mamba`, a matching Python environment can be set up and activated as follows:

```
git clone https://github.com/gorges97/plotQCxMS2/
cd plotQCxMS2
mamba env create -f environment.yml
mamba activate plotQCxMS2
```

### Usage

To visualize the spectrum and fragmentation pathways of a QCxMS2 output file, run in the directory of a finished QCxMS2 calculation:

.. code::
        
        python main.py  

This program will generate a plot of the spectrum and the fragmentation pathways of the QCxMS2 output file.
