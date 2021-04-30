import numpy as np
import matplotlib.pyplot as plt


def W1st(d, coeff):
    #coeff[0] = 0.0    #// W1_im1
    coeff[0] = 1.0-d  #// W1_i
    coeff[1] = d      #// W1_ip1
    return coeff



def W2nd(d, coeff):
    coeff[0] = 0.50*(0.5 - d)*(0.5 - d) #W2_im1 
    coeff[1] = 0.75 - d*d               #W2_i   
    coeff[2] = 0.50*(0.5 + d)*(0.5 + d) #W2_ip1 
    return coeff

def W3rd(d, coeff):
    d2 = d*d
    d3 = d*d*d

    coeff[0] = ( 1.0 - d3 )/6.0 - 0.5*( d-d2 )  # (1-d^3)/6 - (d-d^3)/2
    coeff[1] = 2.0/3.0 - d2 + 0.5*d3            # 2/3 - d^2 + d^3/2
    coeff[2] = 1.0/6.0 + 0.5*( d+d2-d3 )        # 1/6 - (d + d^2 - d^3)/2
    coeff[3] = d3/6.0                           # d^3/6
    return coeff

def W4th(d, coeff):
    d2 = d*d
    d3 = d*d*d
    d4 = d*d*d*d

    coeff[0] = 1.0  /384.0 - 1.0 /48.0*d  + 1.0/16.0*d2 - 1.0/12.0*d3 + 1.0/24.0*d4 # W4_im2
    coeff[1] = 19.0/96.0   - 11.0/24.0*d  + 1.0/4.0* d2 + 1.0/6.0* d3 - 1.0/6.0* d4 # W4_im1
    coeff[2] = 115.0/192.0 - (5.0/8.0)*d2 +(1.0/4.0)*d4                             # W4_i  //ok
    coeff[3] = 19.0 /96.0  + 11.0/24.0*d  + 1.0/4.0* d2 - 1.0/6.0* d3 - 1.0/6.0* d4 # W4_ip1
    coeff[4] = 1.0  /384.0 + 1.0 /48.0*d  + 1.0/16.0*d2 + 1.0/12.0*d3 + 1.0/24.0*d4 # W4_ip2
    return coeff






if __name__ == "__main__":
    fig = plt.figure(1, figsize=(3.487, 2.35))
    
    plt.rc('font',  family='sans-serif')
    #plt.rc('text',  usetex=True)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.rc('axes',  labelsize=8)

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.0)
    gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    #axs.append( plt.subplot(gs[1,0]) )
    
    for ax in axs:
        ax.minorticks_on()

    axs[0].set_xlabel(r'$i$')
    axs[0].set_ylabel(r'$W_n$')

    axs[0].set_xlim((-3, 3))
    axs[0].set_ylim((0.0, 1.0))

    #--------------------------------------------------
    #d = 0.261488
    d = 0.2

    c1 = np.zeros(2)
    x1 = np.array([0, 1])
    c1 = W1st(d, c1)
    print("W1", c1, np.sum(c1))


    dx1 = 0.2 - 0.5
    fx1 = 0.5 - dx1 
    fx2 = 0.5 + dx1 
    print('F1', fx1, fx2)


    c2 = np.zeros(3)
    x2 = np.array([-1, 0, 1])
    c2 = W2nd(d, c2)
    print("W2", c2, np.sum(c2))

    c3 = np.zeros(4)
    x3 = np.array([-1, 0, 1, 2])
    #x3 = np.array([-2, -1, 0, 1]) #WRONG
    c3 = W3rd(d, c3)
    print("W3", c3, np.sum(c3))

    c4 = np.zeros(5)
    x4 = np.array([-2, -1, 0, 1, 2])
    c4 = W4th(d, c4)
    print("W4", c4, np.sum(c4))



    if False:
        dx=0.261488
        dy=-0.393028
        dz=-0.380409

        c4x = np.zeros(5)
        c4x = W4th(dx, c4x)

        c4y = np.zeros(5)
        c4y = W4th(dy, c4y)

        c4z = np.zeros(5)
        c4z = W4th(dz, c4z)

        res = 0.0
        for il in [-2,-1,0,1,2]:
            for jl in [-2,-1,0,1,2]:
                for kl in [-2,-1,0,1,2]:
                    res += c4x[il+2]*c4y[jl+2]*c4z[kl+2]
        print('res', res)



    axs[0].plot(x1, c1, marker='.', color='C0', label='1st')
    axs[0].plot(x2, c2, marker='.', color='C1', label='2nd')
    axs[0].plot(x3, c3, marker='.',  color='C2', label='3rd')
    axs[0].plot(x4, c4, marker='.',  color='C3', label='4th')


    plt.savefig('shapes.pdf')

