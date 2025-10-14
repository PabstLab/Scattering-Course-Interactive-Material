## What you need for interactive sessions and data-analysis templates

It is necessary to have Python installed on your laptop to run the scripts for the interactive sessions. <br>
You will need to have Python and Jupyter installed, along with the following modules:
  - jupyterlab
  - ipywidgets
  - numpy
  - scipy
  - pandas
  - python
  - matplotlib

#### Recommended: Use conda environment!

As an alternative it is recommended to work on an isolated “box”: a **Conda environment**. 
This will ensure to have all the right Python and Jupyter installation and module dependencies in a separated environment.

0. Download and extract the repository into a folder of your choice. Ideally, this should be on the same drive where Python is installed.
1. Install Miniconda (or Anaconda) correctly by _strictly_ following the instructions at https://conda.io/projects/conda/en/latest/user-guide/install/index.html.
2. Open a command-line prompt. On Windows, open either "Anaconda Prompt (Miniconda3)” or “Anaconda Powershell Prompt (Miniconda3)" and create the **Giove** environment by typing:

```
conda env create -f <path-to>/environment.yml
```

3. Activate the conda environment:
```
conda activate Giove
```

4. Open _jupyter lab_ :
```
jupyter lab
```

5. When it is no longer needed, type
```
conda deactivate
```

For more information about managing conda environments go to https://conda.io/projects/conda/en/latest/user-guide/tasks/index.html.
