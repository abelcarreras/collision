import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from elastic import collision_two_masses_1d, collision_two_masses_nd
from itertools import combinations


class Particle:
    def __init__(self,
                 mass=1,
                 size=0.1,
                 position=(0, 0),
                 velocity=(0, 0)):
        self.mass = mass
        self.size = size
        self._position = np.array(position)
        self._velocity = np.array(velocity)

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = np.array(value)

    @property
    def velocity(self):
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        self._velocity = np.array(value)


class ParticleBox:
    """
    System class

    bounds is the size of the box: [xmin, xmax, ymin, ymax]
    """

    def __init__(self,
                 init_state=(Particle(position=[-1, -1.],
                                      velocity=[1.0, 0.31]),
                             Particle(position=[-1, -1.4],
                                      velocity=[0.11, 0.2],
                                      mass=3),
                             Particle(position=[0.5, -1.],
                                      velocity=[0.4, 0.2]),
                             Particle(position=[1, -0.5],
                                      velocity=[0.1, 0.2]),
                             Particle(position=[-1, 1.],
                                      velocity=[0.1, 0.11]),
                             Particle(position=[1., 1.0],
                                      velocity=[-0.14, -0.1],
                                      mass=2),
                             ),
                 bounds=(-2, 2, -2, 2)):
        self.init_state = init_state
        self.state = list(self.init_state)
        self.time_elapsed = 0
        self.bounds = bounds
        self.size = 0.1

    def step(self, dt):
        """step once by dt seconds"""
        self.time_elapsed += dt

        # update positions
        for p in self.state:
            p.position += dt * p.velocity

        # collision on particles
        for pair in combinations(self.state, 2):
            d = np.linalg.norm(pair[0].position - pair[1].position)
            if d < pair[0].size + pair[1].size:
                print('---')
                # mass
                m1 = pair[0].mass
                m2 = pair[1].mass

                # velocity vector
                v1 = pair[0].velocity
                v2 = pair[1].velocity

                # calculate collision
                pair[0].velocity, pair[1].velocity = collision_two_masses_nd(m1, m2, v1, v2)

        # collision in boundaries
        for p in self.state:
            if p.position[0] < self.bounds[0] + p.size:
                p.velocity = np.multiply(p.velocity, [-1, 1])
            if p.position[0] > self.bounds[1] - p.size:
                p.velocity = np.multiply(p.velocity, [-1, 1])
            if p.position[1] < self.bounds[2] + p.size:
                p.velocity = np.multiply(p.velocity, [1, -1])
            if p.position[1] > self.bounds[3] - p.size:
                p.velocity = np.multiply(p.velocity, [1, -1])

        # add gravity
        # self.state[:, 3] -= self.M * self.G * dt


# ------------------------------------------------------------
# set up initial state
box = ParticleBox()
dt = 1. / 15  # 30fps

# ------------------------------------------------------------
# set up figure and animation
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-3.2, 3.2), ylim=(-2.4, 2.4))

# particles holds the locations of the particles
particles, = ax.plot([], [], 'bo', ms=6)

# rect is the box edge
rect = plt.Rectangle(box.bounds[::2],
                     box.bounds[1] - box.bounds[0],
                     box.bounds[3] - box.bounds[2],
                     ec='none', lw=2, fc='none')
ax.add_patch(rect)


def init():
    """initialize animation"""
    global box, rect
    particles.set_data([], [])
    rect.set_edgecolor('none')
    return particles, rect


def animate(i):
    """perform animation step"""
    global box, rect, dt, ax, fig
    box.step(dt)

    ms = int(fig.dpi * 2 * box.size * fig.get_figwidth()
             / np.diff(ax.get_xbound())[0])

    #ms = [int(fig.dpi * 2 * p.size * fig.get_figwidth()/ np.diff(ax.get_xbound()))
    #      for p in box.state]

    # update pieces of the animation
    rect.set_edgecolor('k')
    particles.set_data([p.position[0] for p in box.state],
                       [p.position[1] for p in box.state])
    particles.set_markersize(ms)
    return particles, rect


ani = animation.FuncAnimation(fig, animate, frames=600,
                              interval=10, blit=True, init_func=init)

plt.show()