import constant
import random

# velocity is a 2D vector field [u,v] where :
#   - u is the horizontal component
#   - v is the vertical component

u = [[random.random()-0.5 for _ in range(constant.WIDTH)] for _ in range(constant.HEIGHT)]
v = [[random.random()-0.5 for _ in range(constant.WIDTH)] for _ in range(constant.HEIGHT)]

# If 1 then the cell is empty, if 0 then the cell is filled
grid = [[constant.EMPTY for _ in range(constant.WIDTH)] for _ in range(constant.HEIGHT)]

# Set borders as WALL
for i in range(constant.HEIGHT):
    grid[i][0] = constant.WALL
    grid[i][constant.WIDTH - 1] = constant.WALL
for j in range(constant.WIDTH):
    grid[0][j] = constant.WALL
    grid[constant.HEIGHT - 1][j] = constant.WALL

# Initialize density field
density = [[0.0 for _ in range(constant.WIDTH)] for _ in range(constant.HEIGHT)]
pressure = [[0.0 for _ in range(constant.WIDTH)] for _ in range(constant.HEIGHT)]

def init_density():
    for i in range(constant.HEIGHT):
        for j in range(constant.WIDTH):
            # Example: density proportional to velocity magnitude if not wall
            if grid[i][j] != constant.WALL:
                density[i][j] = (u[i][j]**2 + v[i][j]**2) ** 0.5
            else:
                density[i][j] = 0.0


def advect():
    for i in range(constant.HEIGHT):
        for j in range(constant.WIDTH):
            ... 

def add_force():
    for i in range(constant.HEIGHT):
        for j in range(constant.WIDTH):
            v[i][j] += constant.GRAVITY * constant.TIMESTEP

def boundary_conditions():
    for i in range(constant.HEIGHT):
        for j in range(constant.WIDTH):
            if grid[i][j] == constant.WALL:
                u[i][j] = 0
                v[i][j] = 0


def project():
    for i in range(1, constant.HEIGHT-1):
        for j in range(1, constant.WIDTH-1):
            d = constant.overrelaxation*(u[i+1][j] - u[i][j] + v[i][j+1] - v[i][j])
            s = grid[i+1][j] + grid[i][j] + grid[i][j+1] + grid[i][j]
            u[i][j] += d * grid[i-1][j]/s
            v[i][j] += d * grid[i][j-1]/s
            u[i+1][j] -= d * grid[i+1][j]/s
            v[i][j+1] -= d * grid[i][j+1]/s
    
            # pressure[i][j] += density[i][j] * d / s / constant.TIMESTEP
            # pressure[i][j] = max(0, pressure[i][j])  # Ensure pressure is non-negative

            # # Update velocity according to pressure gradients
            # if grid[i][j] != constant.WALL:
            #     if grid[i+1][j] != constant.WALL:
            #         u[i][j] -= (pressure[i+1][j] - pressure[i][j]) * constant.TIMESTEP
            #     if grid[i][j+1] != constant.WALL:
            #         v[i][j] -= (pressure[i][j+1] - pressure[i][j]) * constant.TIMESTEP
            

init_density()

def update():
    # add_force()
    boundary_conditions()

    advect()
    project()
