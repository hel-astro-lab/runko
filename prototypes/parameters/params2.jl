cfl = 0.45
c_omp = 100.0 #skin depth in cell units
dx = cfl/c_omp

#omp = cfl/c_omp #plasma reaction

delgam = 1.0e-3


vth = sqrt(delgam)
deb = vth*cfl
deb_per_cell = deb/dx

Lx = 50.0*dx
Nx = Lx/dx
T  = 20.0
Nt = T/(cfl*dx)


modes = 19
lambda = Lx/modes
lambda_per_cell = lambda/dx
km = 2*pi/lambda
khat = km*deb*dx/cfl


println("dx                   : $dx")
#println("dt                   : $dt")
println("cells per skin depth : $c_omp")
println("thermal velocity     : $vth")
println("Debye length         : $deb")
println("debye per cell       : $deb_per_cell")

println("modes                : $modes")
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

