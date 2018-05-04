cfl = 0.45
c_omp = 100.0 #skin depth in cell units
delgam = 2.05e-3
Nx = 400
Rm = 7
Nt = 2500


dx = 1/c_omp
dt = cfl*dx

Lx = Nx*dx
lambda = Lx/Rm
lambda_per_cell = lambda/dx

vth = sqrt(delgam)
deb = vth
deb_per_cell = deb/dx

khat = (2pi/lambda)*vth

T = Nt*dt



println("dx                   : $dx")
println("dt                   : $dt")
println("cells per skin depth : $c_omp")
println("thermal velocity     : $vth")
println("Debye length         : $deb")
println("debye per cell       : $deb_per_cell")

println("modes                : $Rm")
println("lambda               : $lambda")
println("lambda per cell      : $lambda_per_cell")
println("khat                 : $khat")

println("Lx                   : $Lx")
println("Nx                   : $Nx")
println("T                    : $T")
println("Nt                   : $Nt")

ompr = 1 + (3/2)*khat^2 + (15/8)*khat^4 + (147/16)*khat^6
ompi = -0.5*sqrt(pi/2)*(1/khat^3 - 6*khat)*exp(-0.5/khat^2 - 3/2 -3*khat^2 -12*khat^4)

println("omp Re: $ompr")
println("omp Im: $ompi")

