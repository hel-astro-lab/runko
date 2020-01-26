from math import sqrt



def rk4(f, t0, y0, t1, n):
    vt = [0] * (n + 1)
    vy = [0] * (n + 1)

    dt = (t1 - t0) / float(n)

    vt[0] = t = t0
    vy[0] = y = y0
    for i in range(1, n + 1):
        k1 = dt * f(t         , y)
        k2 = dt * f(t + 0.5*dt, y + 0.5*k1)
        k3 = dt * f(t + 0.5*dt, y + 0.5*k2)
        k4 = dt * f(t + dt    , y + k3)

        vt[i] = t = t0 + i*dt
        vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6

    return vt, vy
 
def f(t, y):
    return t * sqrt(y)
 
vt, vy = rk4(f, 0, 1, 10, 100)
for t, y in list(zip(vt, vy))[::10]:
    print("%4.1f %10.5f %+12.4e" % (t, y, y - (4 + t * t)**2 / 16))
 




