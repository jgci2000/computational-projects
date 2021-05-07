import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import *


a = 2

def fNoSpring(theta, T, m, g):
    d2x = (T / m) * np.cos(theta)
    d2y = -g * (T / m) * np.sin(theta)
    return np.array([[d2x(theta, T, m)], [d2y(theta, T, m, g)]])

def animationCricle(pos, i):
    return plt.Circle((pos[1, i], pos[0, i]), radius=100)

def animationArrow(pos, i):
    return Arrow(0, 20, pos[0, i], pos[1, i], width=2)

def main():
    tf = 30;
    N = 1000;
    time = np.linspace(0, tf, N, retstep=True);
    
    t = time[0]
    dt = time[1]
    
    T = 10
    g = 9.81
    m = 2.5
    
    pos = np.zeros((2, N))
    vel = np.zeros((2, N))
    theta = np.zeros((1, N))
    
    # Euler Cromer Method
    for i in range(0, N):
        theta[1, i] = np.tan(pos[1, i] / pos[2, i])
        
        vel[:, i + 1] = v[:, i] + fNoSpring(theta, T, m, g) * dt
        pos[:, i + 1] = pos[:, i] + v[:, i + 1] * dt
    
    for i in range(0, N):
        plt.Circle(pos[0, i], pos[1, i])
        plt.Arrow(0, 20, pos[0, i], pos[1, i], width=2)
        
        
    
    

if __name__ == '__main__':
    main()
    
