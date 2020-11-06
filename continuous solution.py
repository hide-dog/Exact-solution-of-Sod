import numpy 
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

g    = 1.4    # specific heat ratio

rho1 = 0.125  # density at area 1
u1   = 0.0    # velocity at area 1
p1   = 0.1    # pressure at area 1

rho4 = 1.0    # density at area 4
p4   = 1.0    # velocity at area 4
u4   = 0.0    # presuure at area 4

t    = 0.2    # time
dt   = 0.001
nt   = int(t/dt)

lenx = 1.0          # length of shock tube
dx   = 0.01         # dx
x0   = 0.5          # position of shock (t = 0.0)
nx   = int(lenx/dx) # number of cell

a1 = (g*p1/rho1) **0.5  # speed of sound (sos) at area 1
a4 = (g*p4/rho4) **0.5  # speed of sound (sos) at area 4

# constant 
a = p4/p1
b = 2*g/(g+1)
c = (g-1)/(g+1) * a1/a4
alpha = 2*g/(g-1)

# initial value (all value is ok, but you should choose 1<m<10)
m_init = 5.0

# output
dir_rho = "rho"
dir_u = "u"
dir_p = "p"

# Function to calculate p3
def cal_p3(m):
    temp = 1 + (2*g*(m**2-1))/(g+1)
    ans  = temp*p1
    return ans

# Function to calculate p2
def cal_p2(p3):
    ans = p3
    return ans

# Function to calculate rho3
def cal_rho3(m, p3):
    temp = 1-c*(m-1/m)
    a3   = a4*temp
    ans  = g*p3 / a3**2
    return ans

# Function to calculate rho2
def cal_rho2(m):
    temp = (g+1)*m**2 / (2+(g-1)*m**2)
    ans  = rho1*temp
    return ans

# Function to calculate u2
def cal_u2(m):
    temp = 2*a1/(g+1)
    ans  = temp*(m-1/m)
    return ans

# Function to calculate u3
def cal_u3(m):
    temp = 1-c*(m-1/m)
    a3   = a4*temp

    temp = a4/a3-1
    ans  = temp * 2*a3 /(g-1)
    return ans, a3

# Function of Mach for newton method
def fx(m):
    ans = b*(m**2 -1) +1 -a*(1-c*(m-1/m))**alpha
    return ans

# f prime of Mach
def fx_prime(m):
    ans = 2*b*m -a*alpha*(1-c*(m-1/m))**(alpha-1) * (-c*(1+1/(m**2)))
    return ans

# Function of main
def main():
    cre_dir(dir_rho)
    cre_dir(dir_u)
    cre_dir(dir_p)

    m = m_init             # Mach number
    
    # newton method
    for i in range(1,20):  
        f  = fx(m)
        fp = fx_prime(m)
        m  = m-f/fp
        #print(m)           # output to check if it is converged

    p3 = cal_p3(m)         # presuure at area 3
    p2 = cal_p2(p3)        # presuure at area 2
    
    rho2 = cal_rho2(m)     # density at area 2
    rho3 = cal_rho3(m, p3) # density at area 3

    u2 = cal_u2(m)         # velocity at area 3
    u3, a3 = cal_u3(m)     # velocity and sos at area 3

    us = m*a1              # velocity of shock wave
    
    # initialization of array
    rho = numpy.zeros(nx)
    u   = numpy.zeros(nx)
    p   = numpy.zeros(nx)
    x   = numpy.zeros(nx)
    for i in range(nx):
        x[i] = dx/2 + dx*i


    for tt in range(nt):
        time = dt*tt
        
        # Storage of values
        for i in range(nx):
            if x0+us*time < x[i]:
                rho[i] = rho1
                u[i]   = u1
                p[i]   = p1
            elif x0+u2*time < x[i] and x[i] <= x0+us*time:
                rho[i] = rho2
                u[i]   = u2
                p[i]   = p2
            elif x0+(u3-a3)*time < x[i] and x[i] <= x0+u2*time:
                rho[i] = rho3
                u[i]   = u3
                p[i]   = p3
            elif x0-a4*time < x[i] and x[i] <= x0+(u3-a3)*time:
                temp   = a4*(2/(g+1) - (g-1)/(g+1)*((x[i]-x0)/(a4*t)))
                rho[i] = rho4 * (temp/a4)**(2/(g-1))
                u[i]   = 2*a4/(g+1)*(1+((x[i]-x0)/(a4*t)))
                p[i]   = p4 * (temp/a4)**(2*g/(g-1))
            else:
                rho[i] = rho4
                u[i]   = u4
                p[i]   = p4
        
        # plot density
        pyplot.plot(x, rho, '-o')                      # plot
        t_in = 'rho_' + '{:.05f}'.format(time) + 's'   # title
        pyplot.title(t_in, fontsize='18')              # title
        pic_name = dir_rho + "/" + str(tt+1) + '.png'  # output name
        pyplot.savefig(pic_name)                       # output
        pyplot.clf()                                   # reset
        
        # plot velocity
        pyplot.plot(x, u, '-o')                      # plot
        t_in = 'u_' + '{:.05f}'.format(time) + 's'   # title
        pyplot.title(t_in, fontsize='18')              # title
        pic_name = dir_u + "/" + str(tt+1) + '.png'  # output name
        pyplot.savefig(pic_name)                       # output
        pyplot.clf()                                   # reset
        
        # plot presuure
        pyplot.plot(x, p, '-o')                      # plot
        t_in = 'p_' + '{:.05f}'.format(time) + 's'   # title
        pyplot.title(t_in, fontsize='18')              # title
        pic_name = dir_p + "/" + str(tt+1) + '.png'  # output name
        pyplot.savefig(pic_name)                       # output
        pyplot.clf()                                   # reset

        print(str(tt+1) + '/' + str(nt))


def cre_dir(dirname):
    try:
        os.mkdir(dirname)
    except:
        pass

# main
main()



