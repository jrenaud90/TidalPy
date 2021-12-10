import os
import subprocess
import warnings

INSTALL_AS_EDITABLE = True
RAISE_CONDA_ERROR = True
cwd = os.getcwd()
conda_file = 'conda_requirements.txt'
pip_file = 'requirements.txt'

# Install conda requirements
try:
    print('Attempting to install requirements using Conda.')
    subprocess.run(['conda', 'install', '--file',  f'{conda_file}', '-c', 'defaults', '-c', 'conda-forge'])
except FileNotFoundError as e:
    if RAISE_CONDA_ERROR:
        raise e
    else:
        warnings.warn('Conda environment not found. It is highly recommended to install'
                      ' TidalPy in an Anaconda environment. Attempting to install with pip alone.')
        subprocess.run(['pip', 'install', 'numpy'])
        subprocess.run(['pip', 'install', 'ecos'])
        subprocess.run(['pip', 'install', 'scs'])

# Install pip requirements
print('Installing additional components using Pip.')
subprocess.run(['pip', 'install', '-r', f'{pip_file}'])

# Install python
if INSTALL_AS_EDITABLE:
    subprocess.run(['pip', 'install', '-e', '.'])
else:
    subprocess.run(['pip', 'install', '.'])
