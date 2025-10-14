# Use and installation

The package works without installation. Just download the package folder and run the script *jupyter lab* as described in the documentation.<br>

You will need to have Python and Jupyter installed, along with the following modules:
  - jupyterlab
  - ipywidgets
  - numpy
  - scipy
  - pandas
  - python
  - matplotlib

#### Recommended: Use conda environment!

As an alternative, it is recommended to run the  Jupyter notebook within a **Conda environment**.
This ensures that the correct Python installation and module dependencies are in place in a separate environment.

0. Download and extract the repository into a folder of your choice. Ideally, this should be on the same drive where Python is installed.
1. Install Miniconda (or Anaconda) correctly by _strictly_ following the instructions at https://conda.io/projects/conda/en/latest/user-guide/install/index.html.
2. Open a command-line prompt. On Windows, open either "Anaconda Prompt (Miniconda3)” or “Anaconda Powershell Prompt (Miniconda3)" and create the **Giove** environment by typing:

```
conda env create -f <path-to>/environment.yml
```

3. Activate the conda environment:
```
conda activate SAS-MoCa
```

4. Open _jupyter lab_ :
```
jupyter lab
```


For more information about managing conda environments go to https://conda.io/projects/conda/en/latest/user-guide/tasks/index.html.
