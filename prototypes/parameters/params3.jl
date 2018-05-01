cfl = 0.45
#c_omp = 100.0 #skin depth in cell units
#dx = cfl/c_omp

delgam = 1.0e-3

lambda_per_cell = 68
#Lx_per_lambda = 19
Lx_per_lambda = 4
modes = Lx_per_lambda
Nx = Lx_per_lambda * lambda_per_cell


lambda_cu = 0.442*cfl #khat=0.45
#lambda_cu = 0.397 #khat=0.5
h = lambda_cu/lambda_per_cell
Lx = dx*Nx

println(h)
dx = h

c_omp = cfl/dx


vth = sqrt(delgam)*cfl
deb = vth
deb_per_cell = deb/dx


km = 2*pi/lambda_cu
khat = km*deb


println("dx                   : $dx")
#println("dt                   : $dt")
println("cells per skin depth : $c_omp")
println("thermal velocity     : $vth")
println("Debye length         : $deb")
println("debye per cell       : $deb_per_cell")

println("modes                : $modes")
println("lambda               : $lambda_cu")
println("lambda per cell      : $lambda_per_cell")
println("khat                 : $khat")

println("Lx                   : $Lx")
println("Nx                   : $Nx")
#println("T                    : $T")
#println("Nt                   : $Nt")

ompr = 1 + (3/2)*khat^2 + (15/8)*khat^4 + (147/16)*khat^6
ompi = -0.5*sqrt(pi/2)*(1/khat^3 - 6*khat)*exp(-0.5/khat^2 - 3/2 -3*khat^2 -12*khat^4)

println("omp Re: $ompr")
println("omp Im: $ompi")

