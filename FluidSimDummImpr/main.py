import FluidCube
import pygame

# Simulation parameters
N = 64
diffusion = 0.0001
viscosity = 0.0001
dt = 0.01

# Initialize fluid cube
fluid = FluidCube.FluidCube(N, diffusion, viscosity, dt)

# Pygame setup
pygame.init()
screen = pygame.display.set_mode((512, 512))
pygame.display.set_caption("Fluid Simulation")
clock = pygame.time.Clock()
running = True

# Interaction key
pressed = {}

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # Add density and velocity at the center
    # fluid.add_density(N//3, N//3, 1000)
    # fluid.add_velocity(N//3, N//3, 100, 100)

    # fluid.add_density(2*N//3, N//3, 1000)
    # fluid.add_velocity(2*N//3, N//3, -100, 100)

    # Step the simulation
    fluid.step()

    # Render density field with white background and grey smoke
    screen.fill((255, 255, 255))
    for i in range(N):
        for j in range(N):
            d = fluid.density[fluid.index(i, j)]
            if d > 255:
                d = 255
            # Smoke is grey: (d, d, d) on white background, so lower d = lighter grey
            color = (255 - d, 255 - d, 255 - d)
            rect = pygame.Rect(i * (512 // N), j * (512 // N), (512 // N), (512 // N))
            pygame.draw.rect(screen, color, rect)

    # If mouse moved, add density and velocity according to mouse deplacement
    mouse_buttons = pygame.mouse.get_pressed()
    if mouse_buttons[0]:  # Left button
        mouse_pos = pygame.mouse.get_pos()
        grid_x = mouse_pos[0] * N // 512
        grid_y = mouse_pos[1] * N // 512
        fluid.add_density(grid_x, grid_y, 500)
        if 'last_mouse_pos' in pressed:
            dx = mouse_pos[0] - pressed['last_mouse_pos'][0]
            dy = mouse_pos[1] - pressed['last_mouse_pos'][1]
            fluid.add_velocity(grid_x, grid_y, dx * 10, dy * 10)
        pressed['last_mouse_pos'] = mouse_pos

    pygame.display.flip()
    clock.tick(60)