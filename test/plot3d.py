import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from numpy import genfromtxt
x = genfromtxt('rk4_output.dat', delimiter=',')

fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
ax = fig.gca(projection='3d')

ax.plot(x[0],x[1],x[2])

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Lorenz")

plt.show()
