import numpy as np
from numpy import arange
import matplotlib
import matplotlib.pyplot as plt

'''
The system parameters are:
m1, m2;
j1, j2, J1, J2, J0;
r1, r2;
t1, t2: (representing thetha1, theta2 - the parameters of pendulum)
t1_dot, t2_dot;
t1_dotdot, t2_dotdot;
L1, L2, l1, l2;

'''
#Values of the parameters are taken from a reference

m1 = 0.3
m2 = 0.075
j1 = 0.0248
j2 = 0.00386
L1 = 0.278
L2 = 0.3
l1 = 0.15
l2 = 0.148
#damping constants 
b1 = 0.0001
b2 = 0.00028
#define gravity 
g = 9.8

J1 = j1 + m1*(l1**2)
J2 = j2 + m2*(l2**2)
J0 = J1 + m2*(L1**2)


#Dynamics of the system F_Xi representing the acceleration of the blocks 
def F_T1(tow_1, tow_2, t1, t2, t1_dot, t2_dot):

    A1 = np.array([-J2*b1, m2*L1*l2*b2*np.cos(t2), -(J2**2)*np.sin(2*t2),
                    -0.5*J2*m2*L1*l2*np.cos(t2)*np.sin(2*t2), J2*m2*L1*l2*np.sin(t2)])
    A2 = np.array([[t1_dot], [t2_dot], [t1_dot*t2_dot], [t1_dot**2], [t2_dot**2]])
    B1 = np.array([J2, -m2*L1*l2*np.cos(t2), 0.5*(m2**2)*(l2**2)*L1*np.sin(2*t2)])
    B2 = np.array([[tow_1], [tow_2], [g]])

    den = J0*J2 + (J2**2)*((np.sin(t2))**2) - (m2**2)*(L1**2)*(l2**2)*((np.cos(t2))**2)

    A = np.dot(A1, A2)
    B = np.dot(B1, B2)

    return np.float64(((A*B)/den))

def F_T2(tow_1, tow_2, t1, t2, t1_dot, t2_dot):

    A1 = np.array([m2*L1*l2*b1*np.cos(t2), -b2*(J0+J2*((np.sin(t2))**2)), 
                   m2*L1*l2*J2*np.cos(t2)*np.sin(2*t2), -0.5*np.sin(2*t2)*(J0*J2+((J2*np.sin(2*t2))**2)),
                   -0.5*((m2*L1*l2)**2)*np.sin(2*t2) ])
    A1 = A1.reshape(5, 1)
    A2 = np.array([t1_dot, t2_dot, t1_dot*t2_dot, t1_dot**2, t2_dot**2]) 
    A2 = A2.reshape(5, 1)
    B1 = np.array([-m2*L1*l2*np.cos(t2), J0+J2*((np.sin(t2)**2)), -m2*l2*np.sin(t2)*(J0+J2*((np.sin(t2))**2))])
    B2 = np.array([[tow_1], [tow_2], [g]])

    den = J0*J2 + (J2**2)*((np.sin(t2))**2) - (m2**2)*(L1**2)*(l2**2)*((np.cos(t2))**2)

    # A = np.dot(A1, A2)
    # B = np.dot(B1, B2)

    A = (A1*A2).sum()
    B = (B1*B2).sum()

    return np.float64(((A*B)/den))

dt = 0.1 
tpoints = arange(0,20,dt)
t1points = []
t2points = []
t1_dotpoints = []
t2_dotpoints = []
t1_dotdotpoints = []
t2_dotdotpoints = []

#initialising the variables of the system 
t1 = 0
t2 = 0 
t1_dot = 0
t2_dot = 0
t1_dotdot = 0
t2_dotdot = 0

#Assuming that constant torque is provided to the motors 
tow_1 = 1
tow_2 = 1

#Solver for Z coordinate considering all parameters in the dynamics equations
for t in tpoints:
    #adding the state variable values to corresponding lists (appending)
    t1points.append(t1)
    t2points.append(t2)
    t1_dotpoints.append(t1_dot)
    t2_dotpoints.append(t2_dot)
    t1_dotdotpoints.append(t1_dotdot)
    t2_dotdotpoints.append(t2_dotdot)
    
    #Numerical solution using Eulers method

    t1 += dt*t1_dot
    t1_dot += dt*t1_dotdot
    t1_dotdot = F_T1(tow_1, tow_2, t1, t2, t1_dot, t2_dot)

    t2 += dt*t2_dot
    t2_dot += dt*t2_dotdot
    t2_dotdot = F_T2(tow_1, tow_2, t1, t2, t1_dot, t2_dot)