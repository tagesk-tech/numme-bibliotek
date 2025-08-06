import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # enables 3D plotting

def heatfd(xl, xr, yb, yt, M, N, D=1.0):
    """
    Forward‐difference solver for u_t = D·u_xx on [xl,xr]×[yb,yt]
    M = number of space steps, N = number of time steps.
    Returns w (including boundaries), x, t.
    """
    # 1) step sizes
    h     = (xr - xl) / M      # Δx
    k     = (yt - yb) / N      # Δt
    m     = M - 1              # interior points in x
    n     = N                  # time‐levels (0…N)
    sigma = D * k / h**2       # CFL number
    if sigma > 0.5:
        raise ValueError(f"Unstable: σ={sigma:.3f} > 0.5")

    # 2) initial & boundary functions
    f = lambda x: np.sin(2*np.pi*x)**2    # u(x,0)
    l = lambda t: 0*t                     # left BC
    r = lambda t: 0*t                     # right BC

    # 3) build coefficient matrix A (size m×m)
    main = (1 - 2*sigma) * np.ones(m)
    off  = sigma        * np.ones(m-1)
    A    = np.diag(main) \
         + np.diag(off, 1) \
         + np.diag(off,-1)

    # 4) grid vectors
    x_interior = xl + np.arange(1, m+1)*h
    t_grid     = yb + np.arange(0, n+1)*k

    # 5) boundary data
    lside = l(t_grid)   # length n+1
    rside = r(t_grid)

    # 6) initialize interior solution W (m×(n+1))
    W = np.zeros((m, n+1))
    W[:, 0] = f(x_interior)

    # 7) time‐marching loop
    for j in range(n):
        # explicit update on interior
        W[:, j+1] = A @ W[:, j]

        # add boundary contributions at each step
        bc = np.zeros(m)
        bc[0]   = sigma * lside[j]
        bc[-1]  = sigma * rside[j]
        W[:, j+1] += bc

    # 8) assemble full solution w with boundaries
    w = np.zeros((m+2, n+1))
    w[0, :]    = lside        # left boundary row
    w[1:-1, :] = W            # interior rows
    w[-1, :]   = rside        # right boundary row

    # 9) full grids including boundaries
    x_full = xl + np.arange(0, m+2)*h
    t_full = t_grid

    return w, x_full, t_full

# — example usage & 3D mesh plot —
w, x, t = heatfd(0, 1, 0, 1, M=10, N=250)
X, T = np.meshgrid(x, t)        # T(i,j)=t[i], X(i,j)=x[j]

fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, w.T, edgecolor='k', linewidth=0.2)
ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('u')
plt.title('Forward‐Difference Solution of Heat Equation')
plt.show()
