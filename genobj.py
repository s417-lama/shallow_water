import numpy as np

# filename="out3d/h_0047.txt"
filename="out3d/h_0060.txt"
# filename="out3d/h_0099.txt"
points = np.loadtxt(filename, delimiter=",")

xs = sorted(set(points[:, 0]))
ys = sorted(set(points[:, 1]))

nx = len(xs)
ny = len(xs)

center_x = np.median(xs)
center_y = np.median(ys)

heights = points[:, 2].reshape((ny, nx))

# vertex
for y in range(ny):
    for x in range(nx):
        print("v {x} {y} {z}".format(x=xs[x] - center_x, y=heights[y][x], z=ys[y] - center_y))

def normalize(v):
    norm = np.linalg.norm(v, ord=1)
    return v / norm

def normal(v1, v2, v3):
    n = np.cross(v2 - v1, v3 - v1)
    z = np.array([0, 1, 0])
    if np.dot(z, n) < 0:
        return -normalize(n)
    else:
        return normalize(n)

def vertex_vec(x, y):
    return np.array([xs[x], heights[y][x], ys[y]])

normals = np.zeros((ny - 1, nx - 1, 2, 3))
for y in range(ny - 1):
    for x in range(nx - 1):
        normals[y][x][0] = normal(vertex_vec(x, y), vertex_vec(x, y + 1), vertex_vec(x + 1, y + 1))
        normals[y][x][1] = normal(vertex_vec(x, y), vertex_vec(x + 1, y), vertex_vec(x + 1, y + 1))

# normal
for y in range(ny):
    for x in range(nx):
        # take average of normals of triangles around a vertex
        ns = []
        if x > 0 and y < ny - 1:
            ns.append(normals[y][x-1][1])
        if x < nx - 1 and y > 0:
            ns.append(normals[y-1][x][0])
        if x > 0 and y > 0:
            ns.append(normals[y-1][x-1][0])
            ns.append(normals[y-1][x-1][1])
        if x < nx - 1 and y < ny - 1:
            ns.append(normals[y][x][0])
            ns.append(normals[y][x][1])
        n = normalize(np.sum(ns, axis=0))
        print("vn {x} {y} {z}".format(x=n[0], y=n[1], z=n[2]))

for y in range(ny - 1):
    for x in range(nx - 1):
        v1 = y*nx+x+1
        v2 = (y+1)*nx+x+1
        v3 = (y+1)*nx+(x+1)+1
        print("f {v1}//{n1} {v2}//{n2} {v3}//{n3}".format(v1=v1, n1=v1, v2=v2, n2=v2, v3=v3, n3=v3))
        v1 = y*nx+x+1
        v2 = y*nx+(x+1)+1
        v3 = (y+1)*nx+(x+1)+1
        print("f {v1}//{n1} {v2}//{n2} {v3}//{n3}".format(v1=v1, n1=v1, v2=v2, n2=v2, v3=v3, n3=v3))
