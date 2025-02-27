help ([[Runko (v5.0.0-dev)

For detailed instructions, go to:
https://runko.readthedocs.io/en/latest/installation.html
]])

whatis("Version: 5.0.0")
whatis("Keywords: PIC code")
whatis("URL: https://runko.readthedocs.io/")
whatis("Description: C++/Python CPU/GPU-enabled Plasma toolkit")

-- # Cray 
load("PrgEnv-cray")
load("craype-x86-milan") -- # or rome
load("cce")
load("craype")

-- # OFI
load("cray-mpich")
load("craype-network-ofi")
load("libfabric")

-- # shared memory support
load("cray-pmi")
load("cray-dsmml")
load("cray-openshmemx")

-- # other tools
load("cray-hdf5")
load("cray-python")
load("cray-libsci")
load("perftools-base")


local user = os.getenv("USER")

-- activate python virtualenv
-- execute {cmd="source /home/" .. user .. "/venvs/runko/bin/activate", modeA={"load"}}
-- new
-- source /wrk-kappa/users/jnattila/venvs/runko-cray2/bin/activate

execute {cmd="source /wrk-vakka/users/" .. user .. "/venvs/runko-cray/bin/activate", modeA={"load"}}


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

setenv("CC",  "cc")
setenv("CXX", "CC")
