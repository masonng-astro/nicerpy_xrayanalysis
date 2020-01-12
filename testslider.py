import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)
fig.canvas.set_window_title('Reaktionsfortschritt')

t0 = 0
t = np.arange(0, t0, .5)
k0 = 0.17
a = np.exp(- k0 * t)
min_t = 5

l, = ax.plot(t, a, lw=3, color='crimson')
plt.axis([0, 20, 0, 1])

axrs = plt.axes([0.25, 0.1, 0.65, 0.03])
srs = Slider(axrs, 'Reaktionsschritte', 0, 20, valinit=min_t)

def update(x):
    base_t = np.arange(0,20,0.5)
    base_y = np.exp(-k0*base_t)
    ax.plot(base_t,base_y,lw=0.1,alpha=0.1)

    t0 = x
    t = np.arange(min_t, t0, .5)
    ax.lines.pop(0)  # remove previous line plot
    ax.plot(t, np.exp(- k0 * t), lw=3, color='crimson')  # plot new one
    fig.canvas.draw()

srs.on_changed(update)

plt.show()
