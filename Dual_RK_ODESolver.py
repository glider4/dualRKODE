"""
Dual ODE Solving via Classical RK4 Method
Two ODE's and two ODE solutions pre-built in and plotted
"""


import matplotlib.pyplot as plt
import math


def main(a, b, n):  # a lower bound, b upper bound, n number of steps

    p = 0   # flag
    h = (b-a) / n   # step size
    
    # Initialize Matrices
    x_val = []
    y_val = []
    t_val = []
    
    x_exact = []
    y_exact = []
    
    # Initial Conditions at t = 0
    t = 0
    x = 1
    y = 0
    
    # Define Functions for Calculation
    # first ODE
    def x_prime(x, y, t):
        sol = x - y + (2 * t) - t**2 - t**3
        return sol
    
    # second ODE
    def y_prime(x, y, t):
        sol = x + y - (4 * (t**2)) + t**3
        return sol
    
    # if exact solution of first ODE kwown, enter here
    def x_exact_func(t):
        sol = math.exp(t) * math.cos(t) + t**2
        return sol
    
    # if exact solution of second ODE known, enter here
    def y_exact_func(t):
        sol = math.exp(t) * math.sin(t) - t**3
        return sol

    while p <= n:  # will stop when reached "n" number of iterations
         
        #RK4 Calcuations below for 2 ODE System
        # 0 part
        k0 = h * x_prime(x, y, t)
        l0 = h * y_prime(x, y, t)
        
        # 1 part
        k1 = h * x_prime( (x + (0.5*h)) , (y + (0.5*k0)) , (t + (0.5*l0)) )
        l1 = h * y_prime( (x + (0.5*h)) , (y + (0.5*k0)) , (t + (0.5*l0)) )
        
        # 2 part
        k2 = h * x_prime( (x + (0.5*h)) , (y + (0.5*k1)) , (t + (0.5*l1)) )
        l2 = h * y_prime( (x + (0.5*h)) , (y + (0.5*k1)) , (t + (0.5*l1)) )
        
        # 3 part
        k3 = h * x_prime( (x+h) , (y+k2) , (t+l2) )
        l3 = h * y_prime( (x+h) , (y+k2) , (t+l2) )
        
        # append current values before changing below
        x_val.append(x)
        y_val.append(y)
        t_val.append(t)
        
        # append current EXACT values for ODE solutions at current "t" val
        x_exact.append(x_exact_func(t))
        y_exact.append(y_exact_func(t))
        
        # change variables for next iteration!
        x += h
        y = y + ( (1/6) * (k0 + (2*k1) + (2*k2) + k3) )
        t = t + ( (1/6) * (l0 + (2*l1) + (2*l2) + l3) )
        
        p += 1  # flag variable to keep while loop to "n" iterations
        
    # Plot results of RK4
    plt.plot(t_val,x_val, '-', label = "X Estimation")  # result
    plt.plot(t_val, x_exact, label = "X Exact")  # exact
    
    plt.plot(t_val, y_val, '-', label = "Y Estimation")  # result
    plt.plot(t_val, y_exact, label = "Y Exact")  # exact
    
    # Plot settings
    plt.xlim(0.0, 1.0)
    plt.ylim(0, 3)
    plt.title("RK4 Dual ODE Method Result, n=10, t in [0,1]")
    plt.xlabel("t")
    plt.ylabel("Value")
    plt.show()
    plt.legend()
    
        
main(0,1,10)  # remember main(lower bound, upper bound, number of steps)
    
        
        
        
        
        
        
        
        
