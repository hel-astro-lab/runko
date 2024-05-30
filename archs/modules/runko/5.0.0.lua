help ([[Runko (v5.0.0-dev)

For detailed instructions, go to:
https://runko.readthedocs.io/en/latest/installation.html
]])

whatis("Version: 5.0.0")
whatis("Keywords: PIC code")
whatis("URL: https://runko.readthedocs.io/")
whatis("Description: C++/Python CPU/GPU-enabled Plasma toolkit")

-- gcc13.2.0 compiler pack (TODO no hdf5)
-- load("GCCcore/13.2.0")
-- load("git/2.42.0-GCCcore-13.2.0")
-- load("CMake/3.27.6-GCCcore-13.2.0")
-- load("Python/3.11.5-GCCcore-13.2.0")
-- load("OpenMPI/4.1.6-GCC-13.2.0")

-- gcc12.3.0 compiler pack
load("GCCcore/12.3.0")
load("OpenMPI/4.1.5-GCC-12.3.0")
load("make/4.4.1-GCCcore-12.3.0")
load("CMake/3.26.3-GCCcore-12.3.0")
load("git/2.41.0-GCCcore-12.3.0-nodocs")
load("HDF5/1.14.0-gompi-2023a")
load("Python/3.11.3-GCCcore-12.3.0")
load("virtualenv/20.23.1-GCCcore-12.3.0")

-- extras
-- load("matplotlib/3.7.2-gfbf-2023a")
-- load("FFmpeg/6.0-GCCcore-12.3.0")
-- load("FFTW.MPI/3.3.10-gompi-2023a")
-- load("sympy/1.12-gfbf-2023a")
-- load("scikit-build/0.17.6-GCCcore-12.3.0")

local user = os.getenv("USER")

-- activate python virtualenv
execute {cmd="source /home/" .. user .. "/venvs/runko/bin/activate", modeA={"load"}}

-- smarter version
-- if { [module-info mode load] || [module-info mode switch2] } {
--     puts stdout "source /home/$USER/venvs/runko/bin/activate;"
-- } elseif { [module-info mode remove] && ![module-info mode switch3] } {
--     puts stdout "deactivate;"
-- }

setenv("RUNKODIR",        "/wrk-vakka/users/" .. user .. "/runko")
prepend_path("PYTHONPATH","/wrk-vakka/users/" .. user .. "/runko")
prepend_path("PYTHONPATH","/wrk-vakka/users/" .. user .. "/runko/lib")
prepend_path("PYTHONPATH","/wrk-vakka/users/" .. user .. "/runko/external/corgi/lib")

setenv("CC",  "mpicc")
setenv("CXX", "mpic++")
