import numpy as np

class FluidCube:
    def __init__(self, size=256, diffusion=0.0, viscosity=0.0, dt=0.1):
        self.size = size
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity

        N = size
        self.s = np.zeros(N*N, dtype=np.float32)
        self.density = np.zeros(N*N, dtype=np.float32)

        self.Vx = np.zeros(N*N, dtype=np.float32)
        self.Vy = np.zeros(N*N, dtype=np.float32)

        self.Vx0 = np.zeros(N*N, dtype=np.float32)
        self.Vy0 = np.zeros(N*N, dtype=np.float32)

    def index(self, x, y):
        res = x + y * self.size
        if res < 0:
            return 0
        elif res >= self.size * self.size:
            return self.size * self.size - 1
        return res

    def add_density(self, x, y, amount):
        self.density[self.index(x, y)] += amount
    
    def add_velocity(self, x, y, amountX, amountY):
        index = self.index(x, y)
        self.Vx[index] += amountX
        self.Vy[index] += amountY
    
    def step(self):
        self.diffuse(1, self.Vx0, self.Vx, self.visc, self.dt, 4, self.size)
        self.diffuse(2, self.Vy0, self.Vy, self.visc, self.dt, 4, self.size)
        
        self.project(self.Vx0, self.Vy0, self.Vx, self.Vy, 4, self.size)

        self.advect(1, self.Vx, self.Vx0, self.Vx0, self.Vy0, self.dt, self.size)
        self.advect(2, self.Vy, self.Vy0, self.Vx0, self.Vy0, self.dt, self.size)

        self.project(self.Vx, self.Vy, self.Vx0, self.Vy0, 4, self.size)

        self.diffuse(0, self.s, self.density, self.diff, self.dt, 4, self.size)
        self.advect(0, self.density, self.s, self.Vx, self.Vy, self.dt, self.size)

    def set_bnd(self, b, x, N):
        for i in range(1, N-1):
            x[self.index(0, i)] = -x[self.index(1, i)] if b == 1 else x[self.index(1, i)]
            x[self.index(N-1, i)] = -x[self.index(N-2, i)] if b == 1 else x[self.index(N-2, i)]

            x[self.index(i, 0)] = -x[self.index(i, 1)] if b == 2 else x[self.index(i, 1)]
            x[self.index(i, N-1)] = -x[self.index(i, N-2)] if b == 2 else x[self.index(i, N-2)]

        x[self.index(0, 0)] = 0.5 * (x[self.index(1, 0)] + x[self.index(0, 1)])
        x[self.index(0, N-1)] = 0.5 * (x[self.index(1, N-1)] + x[self.index(0, N-2)])
        x[self.index(N-1, 0)] = 0.5 * (x[self.index(N-2, 0)] + x[self.index(N-1, 1)])
        x[self.index(N-1, N-1)] = 0.5 * (x[self.index(N-2, N-1)] + x[self.index(N-1, N-2)])

    def lin_solve(self, b, x, x0, a, c, iter, N):
        c_recip = 1.0 / c
        x = x.reshape((N, N))
        x0 = x0.reshape((N, N))
        for k in range(iter):
            x[1:-1, 1:-1] = (x0[1:-1, 1:-1] + a * (
                x[2:, 1:-1] + x[:-2, 1:-1] + x[1:-1, 2:] + x[1:-1, :-2]
            )) * c_recip
            self.set_bnd(b, x.ravel(), N)
        # Reshape back to 1D for consistency
        x.shape = (N*N,)

    def diffuse(self, b, x, x0, diff, dt, iter, N):
        a = dt * diff * (N - 2) * (N - 2)
        self.lin_solve(b, x, x0, a, 1 + 4 * a, iter, N)

    def project(self, velocX, velocY, p, div, iter, N):
        # Reshape to 2D for easier indexing
        velocX2d = velocX.reshape((N, N))
        velocY2d = velocY.reshape((N, N))
        p2d = p.reshape((N, N))
        div2d = div.reshape((N, N))

        # Compute divergence and initialize p
        div2d[1:-1, 1:-1] = -0.5 * (
            velocX2d[2:, 1:-1] - velocX2d[:-2, 1:-1] +
            velocY2d[1:-1, 2:] - velocY2d[1:-1, :-2]
        ) / N
        p2d[1:-1, 1:-1] = 0

        self.set_bnd(0, div2d.ravel(), N)
        self.set_bnd(0, p2d.ravel(), N)
        self.lin_solve(0, p2d.ravel(), div2d.ravel(), 1, 6, iter, N)

        # Subtract gradient of p from velocity
        velocX2d[1:-1, 1:-1] -= 0.5 * (p2d[2:, 1:-1] - p2d[:-2, 1:-1]) * N
        velocY2d[1:-1, 1:-1] -= 0.5 * (p2d[1:-1, 2:] - p2d[1:-1, :-2]) * N

        self.set_bnd(1, velocX2d.ravel(), N)
        self.set_bnd(2, velocY2d.ravel(), N)


    def clamp(self, x, min_val, max_val):
        return np.clip(x, min_val, max_val)

    def advect(self, b, d, d0, velocX, velocY, dt, N):
        dt0 = dt * (N - 2)
        I, J = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        x = I - dt0 * velocX.reshape((N, N))
        y = J - dt0 * velocY.reshape((N, N))

        x = self.clamp(x, 0.5, N - 1.5)
        y = self.clamp(y, 0.5, N - 1.5)

        i0 = np.floor(x).astype(int)
        i1 = i0 + 1
        j0 = np.floor(y).astype(int)
        j1 = j0 + 1

        s1 = x - i0
        s0 = 1 - s1
        t1 = y - j0
        t0 = 1 - t1

        i0 = self.clamp(i0, 0, N - 1)
        i1 = self.clamp(i1, 0, N - 1)
        j0 = self.clamp(j0, 0, N - 1)
        j1 = self.clamp(j1, 0, N - 1)

        d0_2d = d0.reshape((N, N))
        d_new = (
            s0 * (t0 * d0_2d[i0, j0] + t1 * d0_2d[i0, j1]) +
            s1 * (t0 * d0_2d[i1, j0] + t1 * d0_2d[i1, j1])
        )

        # Only update the inner cells
        d2d = d.reshape((N, N))
        d2d[1:-1, 1:-1] = d_new[1:-1, 1:-1]
        self.set_bnd(b, d, N)
