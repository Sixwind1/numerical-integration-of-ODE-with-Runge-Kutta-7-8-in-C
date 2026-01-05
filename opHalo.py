import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('ophalo.txt')

x1, x2, x3 = data[:, 1], data[:, 2], data[:, 3]

fig = plt.figure(figsize=(12, 8))

ax1 = fig.add_subplot(221, projection='3d')
ax1.plot(x1, x2, x3)
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.legend()

ax2 = fig.add_subplot(222)
ax2.plot(x1, x2)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.legend()

ax3 = fig.add_subplot(223)
ax3.plot(x1, x3)
ax3.set_xlabel('X')
ax3.set_ylabel('Z')
ax3.legend()

ax4 = fig.add_subplot(224)
ax4.plot(x2, x3)
ax4.set_xlabel('Y')
ax4.set_ylabel('Z')
ax4.legend()

plt.tight_layout()
plt.show()
