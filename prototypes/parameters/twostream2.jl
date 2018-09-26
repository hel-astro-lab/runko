
cfl = 0.25
dx = 0.05

T_per_Gi = 8 #how many folding times
L_per_l = 2

NxMesh = 20

#-------------------------------------------------- 

u2gam(u) = sqrt(1+u^2)
u2v(u) = u/u2gam(u)

Gi(gam) = 1/sqrt(8*gam^2)

#-------------------------------------------------- 
function setup_twostream(u)

    gam = u2gam(u)
    v = u2v(u)

    println("u:   $u")
    println("v:   $v")

    println("gam: $gam")

    G = Gi(gam)
    println("G:   $G")

    dt = cfl*dx
    Nt = (T_per_Gi / G)/dt
    println("Nt:  $Nt")

    km = sqrt( (3/8) / (gam^3 * v^2) )
    lm = 2pi/km
    println("km:  $km")
    println("lm:  $lm")

    L = lm * L_per_l
    Nx = L/dx

    println("L:   $L")
    println("Nxf: $Nx")
    println("Nx:  $(Nx/NxMesh)")

    Nxsug = signif(Nx/NxMesh,1)
    println("suggested size ")
    println("Nx:$Nxsug NxMesh:$NxMesh dx:$(L/(NxMesh*Nxsug))")

    optimal_size(L, Nxsug, NxMesh)

    return
end


# Brute force search for best grid size
function optimal_size(L, Nx, NxMesh)

    err_min = Inf
    dxn_min = 0
    Nx2  = 0
    NxM2 = 0

    for NxM1 in range(Int(NxMesh) - 50, Int(NxMesh) + 50)
        (NxM1 < 5) && continue
        for Nx1 in range(Int(Nx) - 50, Int(Nx) + 50)
            (Nx1 < 5) && continue

            dxn = L/(NxM1*Nx1)
            err = abs(dxn - dx)/dx

            #println("NxM1: $NxM1 Nx1:$Nx1 dxn:$dxn err:$err")

            if err < err_min
                Nx2     = Nx1
                NxM2    = NxM1
                err_min = err
                dxn_min = dxn
                #println("better!")
            end
        end
    end

    println("optimal")
    println("Nx2:$Nx2 NxM2: $NxM2 dxn:$dxn_min err:$err_min")

end


#-------------------------------------------------- 

println("--------------------------------------------------")
setup_twostream(1.0)

println("--------------------------------------------------")
setup_twostream(2.0)

println("--------------------------------------------------")
setup_twostream(3.0)

println("--------------------------------------------------")
setup_twostream(4.0)
