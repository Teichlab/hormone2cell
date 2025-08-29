# hormone2cell
A toolkit for predicting hormone producing and receiving strength in single cell datasets.


## Installation

```bash

# Create and activate conda environment
conda create -n hormone2cell_env python=3.10 -y
conda activate hormone2cell_env

# Clone repository
git clone https://github.com/Teichlab/hormone2cell.git # Alternatively, download the package via Code → Download ZIP and install it locally. A pip source of the package will be available soon.
cd hormone2cell

# Install in editable mode
pip install -e .
```


## Run example jupyternotebook using jupyter notebook
```bash
# Install Jupyter kernel and Launch JupyterLab
pip install ipykernel
python -m ipykernel install --user --name hormone2cell_env --display-name "hormone2cell_env"
jupyter-lab --no-browser --ip=* --port 8990

```

Then Open Tutorial_pancreas.ipynb and select the hormone2cell_env kernel to run the tutorial. **Further explanations and detailed usage of the tool is provided in this notebook.**


## Citation
