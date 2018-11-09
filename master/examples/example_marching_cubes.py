import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from lib.read_write import *
from skimage import measure
from skimage.draw import ellipsoid
from skimage.data import binary_blobs
from lib.render import *
import sys

# Generate a level set about zero of two identical ellipsoids in 3D
# ellip_base = ellipsoid(6, 10, 16, levelset=True)
# ellip_double = np.concatenate((ellip_base[:-1, ...],
#                                ellip_base[2:, ...]), axis=0)

data = binary_blobs(length=50, blob_size_fraction=0.4, n_dim=3, volume_fraction=0.1, seed=None)
data = scipy.io.loadmat("L:\\Users\\gordon\\00000004 - Running Projects\\20180126 Mito quantification for Gordon\\20180306_results\\3D_seg\\CSM_0b29da6715724d089eb08c6ec15ad193_3DS.mat")['data']
data[data>0] = 1
stack_viewer(data)

# print type(data[0,0,0])
print data.shape
# sys.exit()
# Use marching cubes to obtain the surface mesh of these ellipsoids
verts, faces, normals, values = measure.marching_cubes_lewiner(data, level=None, spacing=(1, 1.0, 1), gradient_direction='descent', step_size=1, allow_degenerate=True, use_classic=False)
print verts
# save_data(verts, "verts", ".\\")
# save_data(faces, "faces", ".\\")
# save_data(normals, "normals", ".\\")
# save_data(values, "values", ".\\")
# save_data(data, "blobs", ".\\")
# Display resulting triangular mesh using Matplotlib. This can also be done
# with mayavi (see skimage.measure.marching_cubes_lewiner docstring).
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Fancy indexing: `verts[faces]` to generate a collection of triangles
mesh = Poly3DCollection(verts[faces])
mesh.set_edgecolor('k')
ax.add_collection3d(mesh)

ax.set_xlabel("x-axis: a = 6 per ellipsoid")
ax.set_ylabel("y-axis: b = 10")
ax.set_zlabel("z-axis: c = 16")

ax.set_xlim(0, 13)  # a = 6 (times two for 2nd ellipsoid)
ax.set_ylim(0, 512)  # b = 10
ax.set_zlim(0, 512)  # c = 16

plt.tight_layout()
plt.show()
