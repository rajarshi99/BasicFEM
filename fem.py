import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import triangle as tr

import sys
from time import time

def u(x,y):
    # return x*y
    # return np.sin(np.pi*x)*np.sin(np.pi*y)
    return x**2*y**2
    # val = np.where(y == 1, 1, 0)
    # return val

def f(x,y):
    # return 0
    # return 2*(np.pi)**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    return -2*(x**2 + y**2)
    # return 0

# We wish to solve -laplacian(u)(x,y) = f(x,y)

def plot_on_grid(u_inp, title, fname, triang):
    cplot = plt.tricontourf(triang, u_inp, levels = 100)
    plt.colorbar(cplot)
    # plt.triplot(triang, 'ko-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.savefig(fname)
    plt.close()

checkpoints = [ time() ]

num_points = 2**int(sys.argv[1])
x_1d = np.linspace(-1,1,num_points)
y_1d = np.linspace(-1,1,num_points)
x,y = np.meshgrid(x_1d,y_1d)
x = x.flatten()
y = y.flatten()
domain = {'vertices' : np.column_stack((x,y))}
for key, val in tr.triangulate(domain).items():
    domain[key] = val

checkpoints.append(time())

triang = mtri.Triangulation(x, y, domain['triangles'])
u_exct = u(x,y)

plot_on_grid(u_exct, "u exact", "u_exact.png", triang)

checkpoints.append(time())

K_glob_dict = {}
f_glob_dict = {}
for vert_ids in domain['triangles']:
    x_vert = x[vert_ids]
    x02 = x_vert[0] - x_vert[2]
    x12 = x_vert[1] - x_vert[2]

    y_vert = y[vert_ids]
    y02 = y_vert[0] - y_vert[2]
    y12 = y_vert[1] - y_vert[2]

    detJ = x02*y12 - x12*y02
    B = np.array([[y12, -y02, -y12+y02],
                  [-x12, x02, x12-x02]]) / detJ
    K_elem = np.matmul(B.T,B) * detJ / 2  # Integrated over the elem?

    x_mid = np.sum(x_vert) / 3
    y_mid = np.sum(y_vert) / 3
    f_mid = f(x_mid,y_mid)
    f_int = f_mid * detJ / 6

    for v_l_id,v_g_id in enumerate(vert_ids):
        if v_g_id in f_glob_dict:
            f_glob_dict[v_g_id] += f_int
        else:
            f_glob_dict[v_g_id] = f_int
        for u_l_id,u_g_id in enumerate(vert_ids):
            if (v_g_id,u_g_id) in K_glob_dict:
                K_glob_dict[(v_g_id,u_g_id)] += K_elem[v_l_id,u_l_id]
            else:
                K_glob_dict[(v_g_id,u_g_id)] = K_elem[v_l_id,u_l_id]

checkpoints.append(time())

vert_known_list = np.where(domain['vertex_markers'].reshape(-1) == 1)[0].tolist()
n_known = len(vert_known_list)
u_known = u(x[vert_known_list], y[vert_known_list])

vert_unknown_list = np.where(domain['vertex_markers'].reshape(-1) == 0)[0].tolist()
n_unknown = len(vert_unknown_list)

checkpoints.append(time())

K_glob = np.zeros((n_unknown,n_unknown))
K_known = np.zeros((n_unknown,n_known))
for v_id,u_id in K_glob_dict:
    if v_id in vert_unknown_list:
        v_ind = vert_unknown_list.index(v_id)
        if u_id in vert_unknown_list:
            u_ind = vert_unknown_list.index(u_id)
            K_glob[v_ind,u_ind] = K_glob_dict[(v_id,u_id)]
        else:
            u_ind = vert_known_list.index(u_id)
            K_known[v_ind,u_ind] = K_glob_dict[(v_id,u_id)]

f_glob = np.zeros(n_unknown)
for v_id in f_glob_dict:
    if v_id in vert_unknown_list:
        v_ind = vert_unknown_list.index(v_id)
        f_glob[v_ind] = f_glob_dict[v_id]
f_glob -= np.dot(K_known,u_known)

checkpoints.append(time())
    
u_unknown = np.linalg.solve(K_glob,f_glob)        

checkpoints.append(time())

u_sol = np.zeros(n_known + n_unknown)
u_sol[vert_known_list] = u_known
u_sol[vert_unknown_list] = u_unknown

plot_on_grid(u_sol, f"u fem {num_points}", "u_fem.png", triang)

plot_on_grid(u_sol - u_exct, f"u err {num_points}", "u_err.png", triang)

print(num_points, np.linalg.norm(u_sol - u_exct), end = "")

t_beg = checkpoints[0]
for t in checkpoints:
    print(" ", t - t_beg, end = "")
print(" ", t_beg)
    

    
