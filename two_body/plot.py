
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = np.loadtxt("orbit.txt")

ax = data[:,0]
ay = data[:,1]

bx = data[:,3]
by = data[:,4]

fig, ax_plot = plt.subplots()

ax_plot.set_aspect('equal')
ax_plot.grid(True)

margin = 0.2
xmin = min(ax.min(), bx.min()) - margin
xmax = max(ax.max(), bx.max()) + margin
ymin = min(ay.min(), by.min()) - margin
ymax = max(ay.max(), by.max()) + margin

ax_plot.set_xlim(xmin, xmax)
ax_plot.set_ylim(ymin, ymax)

# Trajectory traces
trail_a, = ax_plot.plot([], [], label="Body A")
trail_b, = ax_plot.plot([], [], label="Body B")

# Current positions
point_a, = ax_plot.plot([], [], 'o')
point_b, = ax_plot.plot([], [], 'o')

ax_plot.legend()

def update(frame):
    trail_a.set_data(ax[:frame+1], ay[:frame+1])
    trail_b.set_data(bx[:frame+1], by[:frame+1])

    point_a.set_data([ax[frame]], [ay[frame]])
    point_b.set_data([bx[frame]], [by[frame]])

    return trail_a, trail_b, point_a, point_b

ani = FuncAnimation(
    fig,
    update,
    frames=len(ax),
    interval=20,
    blit=False
)

ani.save("orbit.gif", fps=30)