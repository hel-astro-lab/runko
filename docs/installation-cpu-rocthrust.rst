The cpu backend uses rocThrust's cpu backend which most likely will not be available as a package.
Thus we have to install it manually

.. code:: shell

   git clone --no-checkout --depth=1 --filter=tree:0 https://github.com/ROCm/rocm-libraries.git
   cd rocm-libraries
   git sparse-checkout init --cone
   git sparse-checkout set projects/rocthrust
   git checkout develop
   cd projects/rocthrust
   export ROCTHRUST_INSTALL_PREFIX="$PWD/rocthrust-install"
   cmake -Bbuild -DTHRUST_DEVICE_SYSTEM=CPP -DLINK_HIP_DEVICE_LIBS=OFF -DCMAKE_INSTALL_PREFIX="$ROCTHRUST_INSTALL_PREFIX" .
   make -C build install
