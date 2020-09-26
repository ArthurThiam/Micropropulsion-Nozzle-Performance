from cx_Freeze import setup, Executable
import os
import scipy

includefiles_list = []
scipy_path = os.path.dirname(scipy.__file__)
includefiles_list.append(scipy_path)

build_exe_options = {

    "packages": ['scipy.signal', 'scipy.sparse', 'scipy.integrate', 'scipy.sparse.csgraph_validation'],
}

setup(name="Microprop Performance",
      version="1.0",
      description="Performance curves for various micro-nozzle and chamber configurations",
      executables=[Executable("Performance.py")])
