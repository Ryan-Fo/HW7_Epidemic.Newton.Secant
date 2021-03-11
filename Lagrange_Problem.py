import numpy as np
from scipy.misc import derivative

#Global variables
G = 6.674e-11
M = 1.989e30
m = 5.972e24
R = 1.496e11
omega = 1.991e-7

def f_1(r):
    """equation of the location of L1"""
    return G*M*(R-r)**2- G*m*r**2- omega**2*r**3*(R-r)**2 

def f_1prime(r):#overdid the math for this and f_2 prime
    """equation for the derivative of location of L1"""
    return r**2*omega**2*(-5*r**2 + 8*r*R - 3*R**2) - 2*G * (m*r+M*(R - r))

def f_2(r):
    """equation of the location of L2"""
    return G*M*(R-r)**2 + G*m*r**2 - omega**2*r**3*(R-r)**2 

def f_2prime(r):
    """equation for the derivative of location of L2"""
    return 2*G*(m*r - M*(R - r)) + r**2 * omega**2 *(-5*r**2 + 8*r*R - 3 * R**2)

def Newtons_Method(func = None, deriv = None, firstguess = None, max_iterations = 100):
    """Uses Newton method to find root of function
    Input: Function,Derivative, a guess for location
    Output: root value
    """
    guess = firstguess
    for i in range(max_iterations):
        guess = guess - (func(guess) / deriv(guess))
    return guess

def Secant_Method(func = None, deriv = None, firstguess = None, secondguess = None, max_iterations = 100):
    """Uses Secant method to find root of function
    Input: Function,Derivative, 1st and 2nd guess
    Output: root value
    """
    x1, x2 = firstguess, secondguess
    
    for i in range(max_iterations):
        term1 = x2 - func(x2)
        numerator = (x2 - x1)
        denominator = func(x2) - func(x1)
        newpoint = x2 - (func(x2) * (x2 - x1))/(func(x2) - func(x1))
        if denominator <= 0:
            x1 = x2
            x2 = newpoint 
        else:
            x2 = x1
            x1 = newpoint
    return newpoint

    
if __name__ == "__main__":
    L1_newtons = Newtons_Method(func = f_1, deriv = f_1prime, firstguess = 1)
    L2_newtons = Newtons_Method(func = f_2,deriv = f_2prime, firstguess = 1)
    
    L1_secant = Secant_Method(func = f_1, deriv = f_1prime, firstguess = 1e9, secondguess = 1e11)
    L2_secant = Secant_Method(func = f_2, deriv = f_2prime, firstguess = 1e9, secondguess = 1e11)
    
    print('values:',L1_newtons,L2_newtons,L1_secant,L2_secant,'\n')
    print("Using Newtons Method the location at L1 in km:",(R - L1_newtons)/1000)
    print("Using Newtons Method the location at L2 in km:",(L2_newtons - R)/1000)
    
    print("Using the Secant Method the location at L1 in km:",((R - L1_secant)/1000))
    print("Using the Secant Method the location at L2 in km:",((L2_secant - R)/1000))