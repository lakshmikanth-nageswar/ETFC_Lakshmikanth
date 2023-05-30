
import numpy as np
from numpy import arange
import matplotlib
import matplotlib.pyplot as plt

'''
The system parameters are:
m1, m2, m3;
k1, k2, k3;
c1, c2, c3;
f1, f2, f3;

'''
#random values of the system parameters

m1=1
m2=1
m3=1

k1=1
k2=1
k3=1

c1=1
c2=1
c3=1

#The input for the system is f1, f2, f3
#Input values:
f1=2
f2=2
f3=2


#Dynamics of the system F_Xi representing the acceleration of the blocks 
def F_X1(f1, x1, x2, x1_dot, x2_dot):
    return (1/m1)*(f1-k1*x1+k2*(x2-x1)-c1*x1_dot+c2*(x2_dot-x1_dot))

def F_X2(f2, x1, x2, x3, x1_dot, x2_dot, x3_dot):
    return (1/m2)*(f2-k3*(x2-x3)-k2*(x2-x1)-c2*(x2_dot-x1_dot)-c3*(x2_dot-x3_dot))

def F_X3(f3, x2, x3, x2_dot, x3_dot):
    return (1/m3)*(f3-k3*(x3-x2)-c3*(x3_dot-x2_dot))

dt = 0.1
tpoints = arange(0,20,dt)
x1points = []
x1_dotpoints = []
x2points = []
x2_dotpoints = []
x3points = []
x3_dotpoints = []
x2_dotpoints = []
x1_dotdotpoints = []
x2_dotdotpoints = []
x3_dotdotpoints = []

#initialising the variables of the system 
x1 = 0
x1_dot = 0
x1_dotdot = 0 

x2 = 0
x2_dot = 0
x2_dotdot = 0 

x3 = 0
x3_dot = 0
x3_dotdot = 0 

  

#Solver for Z coordinate considering all parameters in the dynamics equations
for t in tpoints:
    #adding the state variable values to corresponding lists (appending)
    x1points.append(x1)
    x1_dotpoints.append(x1_dot)
    x2points.append(x2)
    x2_dotpoints.append(x2_dot)
    x3points.append(x3)
    x3_dotpoints.append(x3_dot)
    x1_dotdotpoints.append(x1_dotdot)
    x2_dotdotpoints.append(x2_dotdot)
    x3_dotdotpoints.append(x3_dotdot)

    #Numerical solution using Eulers method

    x1 += dt*x1_dot
    x1_dot += dt*x1_dotdot
    x1_dotdot = F_X1(f1, x1, x2, x1_dot, x2_dot)

    x2 += dt*x2_dot
    x2_dot += dt*x2_dotdot
    x2_dotdot = F_X2(f2, x1, x2, x3, x1_dot, x2_dot, x3_dot)

    x3 += dt*x3_dot
    x3_dot += dt*x3_dotdot
    x3_dotdot = F_X3(f3, x2, x3, x2_dot, x3_dot)


plt.plot(tpoints, x1points, color='r', label='x1')
plt.plot(tpoints, x2points, color='g', label='x2')
plt.plot(tpoints, x3points, color='b', label='x3')
plt.xlabel("Time")
plt.ylabel("Magnitude")
plt.title("Displacement vs time for step input")
plt.legend()
plt.show()

plt.plot(tpoints, x1_dotpoints, color='r', label='x1')
plt.plot(tpoints, x2_dotpoints, color='g', label='x2')
plt.plot(tpoints, x3_dotpoints, color='b', label='x3')
plt.xlabel("Time")
plt.ylabel("Magnitude")
plt.title("Velocity vs time for step input")
plt.legend()
plt.show()

plt.plot(tpoints, x1_dotdotpoints, color='r', label='x1')
plt.plot(tpoints, x2_dotdotpoints, color='g', label='x2')
plt.plot(tpoints, x3_dotdotpoints, color='b', label='x3')
plt.xlabel("Time")
plt.ylabel("Magnitude")
plt.title("Acceleration vs time for step input")
plt.legend()
plt.show()


