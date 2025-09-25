import pygame
import simulator
import constant

pygame.init()
screen = pygame.display.set_mode((800, 600))
pygame.display.set_caption("Vector Field Simulation")

running = True
clock = pygame.time.Clock()

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill((255, 255, 255))  # Clear screen with white

    # Draw vector field
    for i in range(constant.WIDTH):
        for j in range(constant.HEIGHT):
            x = int((i+1/2) * screen.get_width() / constant.WIDTH)
            y = int((j+1/2) * screen.get_height() / constant.HEIGHT)
            u = simulator.u[i][j]
            v = simulator.v[i][j]
            scale = 50  # Adjust scale for visibility
            end_x = int(x + u * scale)
            end_y = int(y + v * scale)
            pygame.draw.aaline(screen, (0, 0, 0), (x, y), (end_x, end_y), 1)
            pygame.draw.circle(screen, (0, 0, 0), (x, y), 2)

    pygame.time.delay(100)
    simulator.update()

            
    pygame.display.flip()
    clock.tick(60)

pygame.quit()