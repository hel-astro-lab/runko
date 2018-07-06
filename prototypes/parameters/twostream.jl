cfl = 0.45
c_omp = 100.0 #skin depth in cell units
delgam = 1.0e-3
Nx = 400
Rm = 2
Nt = 10000
gamma = 0.25


dx = 1/c_omp
dt = cfl*dx

#box size and lambda length
Lx = Nx*dx
lambda = Lx/Rm
lambda_per_cell = lambda/dx

#thermal speed & debye length
vth = sqrt(delgam)
deb = vth
deb_per_cell = deb/dx

#beam velocity
vb = gamma/vth

kbar = 2pi/lambda
khat = (2pi/lambda)*vth

T = Nt*dt



println("dx                   : $dx")
println("dt                   : $dt")
println("cells per skin depth : $c_omp")
println("gamma                : $gamma")
println("thermal velocity     : $vth")
println("Debye length         : $deb")
println("debye per cell       : $deb_per_cell")
println("vb                   : $vb")

println("modes                : $Rm")
println("lambda               : $lambda")
println("lambda per cell      : $lambda_per_cell")
println("kbar                 : $kbar")
println("khat                 : $khat")

println("Lx                   : $Lx")
println("Nx                   : $Nx")
println("T                    : $T")
println("Nt                   : $Nt")

ompr = 1 + (3/2)*khat^2 + (15/8)*khat^4 + (147/16)*khat^6
ompi = -0.5*sqrt(pi/2)*(1/khat^3 - 6*khat)*exp(-0.5/khat^2 - 3/2 -3*khat^2 -12*khat^4)

println("omp Re: $ompr")
println("omp Im: $ompi")

