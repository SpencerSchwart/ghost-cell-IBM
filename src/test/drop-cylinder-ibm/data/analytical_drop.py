
import matplotlib.pyplot as plt
import numpy as np
import sys

### constants ###
lc0 = np.sqrt(2)/2
rc  = 0.5               # cylinder radius
r0  = 0.5               # inital drop radius
thetas = 120*np.pi/180

### solving variables ###

# s0 = initial area of the droplet
def get_s0 (r, lc):
    term1 = np.pi*(r0**2)
    term2 = (rc**2)*np.arccos(((rc**2)-(r**2)+(lc**2))/(2*r*lc))
    term3 = (r0**2)*np.arccos((-(rc**2)+(r**2)+(lc**2))/(2*r0*lc))
    term4 = np.sqrt(((rc+r-lc)*(rc-r+lc)*(rc+r+lc)*(-r+r+lc))/4)

    return term1 - term2 - term3 + term4

def alpha (r, lc):
    return np.arccos(((rc**2)-(r**2)+(lc**2))/(2*r*lc))

# R = final radius of droplet
def R (lc):
    term1 = rc*np.cos(thetas)
    term2 = np.sqrt((rc**2)*(np.cos(thetas)**2)-(rc**2)+(lc**2))

    return term1 + term2

def ffunc (lc):
    r = R(lc)
    #print(r, lc)
    term1 = get_s0 (r, lc0)
    term2 = (alpha(r, lc) + thetas)*(r**2)
    term3 = alpha(r, lc) * (rc**2)
    term4 = r*rc*np.sin(thetas)

    return term1 - term2 + term3 - term4

def find_lc():
    itr = 0 
    error = 1e30
    maxItr = 100 
    maxError = 1e-10
    minlc = 1.90001*rc
    maxlc = 1.9*rc

    # (min*rc, max*rc)
    # theta = 30  (0.68, 0.64)
    # theta = 60  (1.1, 1)
    # theta = 90  (,)
    # theta = 120 (1.7,1.65)
    # theta = 150 (1.91,1.90001

    newlc = 0
    while np.abs(error) > maxError and itr < maxItr:
        newlc = (maxlc + minlc)/2
        error = ffunc(newlc)
        print (itr, " error = ", error)
        if error > 0:
            maxlc = newlc
        else:
            minlc = newlc
        itr += 1

    return R(newlc), newlc

        


#rf, lcf = find_lc()
lcf = float(sys.argv[1])*rc

#x_vals = np.linspace(0.64*rc,0.68*rc,1000)
#y_vals = [ffunc(x) for x in x_vals]

#plt.plot (x_vals, y_vals)
#plt.show()

error = ffunc(lcf)
rf = R(lcf)

print ("Error = ", error)
print ("R =", rf, "Lc =", lcf, "S0 =", get_s0(rf, lcf), "@Theta = ", thetas)
