
c   = 1.0
cfl = 0.45


#grid dimensions
#dx = 0.0032
#dx = 0.0048
#dx = 0.05
#dx = 1.0/700.0
dx = 0.002
dx = 0.01

dt = cfl*dx

Nx = 23*12
Nt = 10000

#plasma parameters
theta = 0.001
vb = 0.0
n0 = 1.0

# problem specific parameters
modes = 2


##################################################

Lx = Nx*dx
T  = Nt*dt

#skin depth
skin = cfl/dx
#skin = c/dx
#skin = cfl*c/dx

#Debye length
vth = sqrt(theta)*cfl
#deb = sqrt(theta/n0)
deb = vth/sqrt(n0)

deb_per_cell = deb/dx

#perturbation wavelengths and modes
lambda = Lx/modes
km = 2*pi/lambda
khat = km*deb
khat = km*vth


println("dx                   : $dx")
println("dt                   : $dt")
println("cells per skin depth : $skin")
println("thermal velocity     : $vth")
println("Debye length         : $deb")
println("debye per cell       : $deb_per_cell")
println("lambda               : $lambda")
println("khat                 : $khat")

println("Lx                   : $Lx")
println("T                    : $T")

ompr = 1 + (3/2)*khat^2 + (15/8)*khat^4 + (147/16)*khat^6
ompi = -0.5*sqrt(pi/2)*(1/khat^3 - 6*khat)*exp(-0.5/khat^2 - 3/2 -3*khat^2 -12*khat^4)

println("omp Re: $ompr")
println("omp Im: $ompi")


if true
    println("##################################################")
    #qe for khat=0.3 = -3.65354426E-05
    #qe for khat=0.45= -2.06708339E-08

    #dx = 0.0048
    #dx = 0.0032
    #dx = 0.1
    #dx = 0.00144
    dx = 0.002

    c = 1.0
    c_omp = 1.0/dx
    cfl = 0.45
    #c_omp = cfl/dx

    println("c_omp $c_omp")
    #omp = c/c_omp
    omp = cfl/c_omp
    println("omp $omp")

    gamma0 = 1.0
    ppc0 = 1
    me = 1.0
    mi = 1.0

    #no sqrt because qe/me = 1 
    qe=-(omp^2*gamma0)/((ppc0*.5)*(1+me/mi)) 
    println("qe $qe")
    	
    qi=-qe 
    println("qi $qi")
    	
    me=me*abs(qi)
    mi=mi*abs(qi)
    println("me/mi $me / $mi")
    	
    qme=qe/me	# equal to -1/me in input
    qmi=qi/mi	# equal to  1/mi in input
    println("qme/qmi $qme / $qmi")

end



