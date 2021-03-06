from __future__ import division
from matplotlib import animation
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

fig = plt.figure()
ax = fig.add_subplot(1,1,1,xlim=(0,1),ylim=(-1,1))
line, = ax.plot([],[],lw=2)

print 'Input the position on x-axis of disturbution 1 ->',
wxp_1 = float(raw_input())
print 'Input the position on x-axis of disturbution 2 ->',
wxp_2 = float(raw_input())
print 'Input the amplitude of disturbution 2 (amplitude of disturbution 1 as 1)->',
Amp = float(raw_input())

def init():
    x = np.linspace(0,1,101)
    y = np.exp(-1000*(x-wxp_1)**2) + Amp*np.exp(-1000*(x-wxp_2)**2) 
    y[0] = 0
    y[-1] = 0
    line.set_data(x,y)
    
    return line,

def iteration(num):

    x = np.linspace(0,1,101)

    y_now = np.exp(-1000*(x-wxp_1)**2) + Amp*np.exp(-1000*(x-wxp_2)**2)
    y_now[0] = 0
    y_now[-1] = 0

    y_before = deepcopy(y_now)
    y_after = np.zeros(101)

    for j in range(num):
        for i in range(101):
            if i== 0 or i== 100:
                y_after[i] = 0
            else:
                y_after[i] = - y_before[i] + y_now[i+1] + y_now[i-1]
        y_before = deepcopy(y_now)
        y_now = deepcopy(y_after)

    return y_now

    
def animate(i):

    x = np.linspace(0,1,101)
    y = iteration(i)
    line.set_data(x,y)
    
    return line,

anim1=animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=5)
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.title('Wave on a string')
plt.grid(True)
plt.show()