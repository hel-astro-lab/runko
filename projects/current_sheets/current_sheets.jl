#Basic study of Harris Current sheets

n0 = 1.0
delta = 0.05
B0 = 1.0

ymin = 0.0
ymax = 1.0
Ly = ymax-ymin

function cosh_distribution(delta)
    
    y = ymin + (ymax-ymin)*rand()
    n = 0.0
    nmax = 1.1*n0
    if y < Ly/2
        n = n0*cosh((y-Ly/4)/delta)^(-2) + 0.1*n0
    else
        n = n0*cosh((y-3*Ly/4)/delta)^(-2) + 0.1*n0
    end

    x = nmax*rand()
    if x > n
        return cosh_distribution(delta)
    end

    return y
end


function cosh_func(y, delta)

    n = 0.0
    if y < Ly/2
        n = n0*cosh((y-Ly/4)/delta)^(-2) + 0.1*n0
    else
        n = n0*cosh((y-3*Ly/4)/delta)^(-2) + 0.1*n0
    end

    return n
end

Ny = 100
ygrid = linspace(ymin, ymax, Ny)
ngrid = zeros(Ny)
for i = 1:Ny
    ngrid[i] = cosh_func(ygrid[i], delta)
end
p = plot(ygrid, ngrid, "b-")
display(p)

N = 10000
ngrid2 = zeros(N)
for i = 1:N
    ngrid2[i] = cosh_distribution(delta)
end
xhist, yhist = hist(ngrid2, 50)
yhist = yhist ./ maximum(yhist)
p = oplot(xhist, yhist, "k--")


###################################################
# B field solution to Harris sheet

function Harris_B(y, delta, alpha)

    Bx = 0.0
    if y < Ly/2
        Bx = -B0*tanh((y-Ly/4)/delta)
    else
        Bx = B0*tanh((y-3*Ly/4)/delta)
    end

    return Bx
end


Bgrid = zeros(Ny)
for i = 1:Ny
    Bgrid[i] = Harris_B(ygrid[i], delta, 0.0)
end
p = oplot(ygrid, Bgrid, "r-")
display(p)



