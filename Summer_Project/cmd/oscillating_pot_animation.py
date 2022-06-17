from matplotlib import animation
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# x_data = []
# y_data = []

fig, ax = plt.subplots()
ax.set_xlim(-20,20)
ax.set_ylim(-100,100)
ax.set_xlabel("Postion")
ax.set_ylabel("Potential (Epsilon_0 = 5, frequency(omega) = 2)")
line, = ax.plot(0, 0)

def potential(time):
    x = np.arange(-20, 20, 0.1)
    y = -1/abs(x) + 5*x*np.cos(2*time)
    x_data = list(x)
    y_data = list(y)
    ax.set_title("time={}".format(time))
    line.set_xdata(x_data)
    line.set_ydata(y_data)

    return line,

animation = FuncAnimation(fig, potential, np.arange(0, 20, 0.25))
animation.save('potential_animation.gif')