import os

# Grid constant
WIDTH = 10
HEIGHT = 10

SCREEN_WIDTH = os.environ.get('SCREEN_WIDTH', 800)
SCREEN_HEIGHT = os.environ.get('SCREEN_HEIGHT', 400)

# Physics constants
GRAVITY = 9.81
TIMESTEP = 0.01
h = 1  # Grid cell size
overrelaxation = 1.9

# Tile types
WALL = 0
EMPTY = 1