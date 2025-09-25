class FluidCube:

    def __init__(self, size=256, diffusion=0.0, viscosity=0.0, dt=0.1):
        self.size = size
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity

        N = size
        self.s = [0.0 for _ in range(N*N)]
        self.density = [0.0 for _ in range(N*N)]

        self.Vx = [0.0 for _ in range(N*N)]
        self.Vy = [0.0 for _ in range(N*N)]

        self.Vx0 = [0.0 for _ in range(N*N)]
        self.Vy0 = [0.0 for _ in range(N*N)]

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
        cRecip = 1.0 / c

        for j in range(1, N-1):
            for i in range(1, N-1):
                x[self.index(i, j)] = (x0[self.index(i, j)] 
                                        + a * (
                                            x[self.index(i+1, j)] 
                                            + x[self.index(i-1, j)] 
                                            + x[self.index(i, j+1)] 
                                            + x[self.index(i, j-1)]
                                            )
                                        ) * cRecip
        self.set_bnd(b, x, N)
    
    def diffuse(self, b, x, x0, diff, dt, iter, N):
        a = dt * diff * (N - 2) * (N - 2)
        self.lin_solve(b, x, x0, a, 1 + 6 * a, iter, N)

    def project(self, velocX, velocY, p, div, iter, N):
        for j in range(1, N-1):
            for i in range(1, N-1):
                div[self.index(i, j)] = -0.5 * (
                    velocX[self.index(i+1, j)] 
                    - velocX[self.index(i-1, j)] 
                    + velocY[self.index(i, j+1)] 
                    - velocY[self.index(i, j-1)]
                ) / N
                p[self.index(i, j)] = 0

        self.set_bnd(0, div, N)
        self.set_bnd(0, p, N)
        self.lin_solve(0, p, div, 1, 6, iter, N)

        for j in range(1, N-1):
            for i in range(1, N-1):
                velocX[self.index(i, j)] -= 0.5 * (p[self.index(i+1, j)] - p[self.index(i-1, j)]) * N
                velocY[self.index(i, j)] -= 0.5 * (p[self.index(i, j+1)] - p[self.index(i, j-1)]) * N

        self.set_bnd(1, velocX, N)
        self.set_bnd(2, velocY, N)

    def clamp(self, x, min_val, max_val):
        return max(min(x, max_val), min_val)

    def advect(self, b, d, d0, velocX, velocY, dt, N):
        dt0 = dt * (N - 2)
        for j in range(1, N-1):
            for i in range(1, N-1):
                tmp1 = dt0 * velocX[self.index(i, j)]
                tmp2 = dt0 * velocY[self.index(i, j)]
                x = i - tmp1
                y = j - tmp2

                x = self.clamp(x, 0.5, N + 0.5)
                i0 = int(x)
                i1 = i0 + 1

                y = self.clamp(y, 0.5, N + 0.5)
                j0 = int(y)
                j1 = j0 + 1

                s1 = x - i0
                s0 = 1 - s1
                t1 = y - j0
                t0 = 1 - t1

                d[self.index(i, j)] = (
                    s0 * (t0 * d0[self.index(i0, j0)] + t1 * d0[self.index(i0, j1)]) +
                    s1 * (t0 * d0[self.index(i1, j0)] + t1 * d0[self.index(i1, j1)])
                )
        self.set_bnd(b, d, N)
        








