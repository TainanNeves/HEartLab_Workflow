# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 11:35:13 2023

@author: HEartLab
"""

import vedo
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import imageio

# Função para criar cada quadro da animação
def update_frame(frame):
    # Limpe a cena
    scene.clear()

    # Rotação do objeto 3D
    mesh.rotate_y(angle_increment)

    # Adicione o objeto à cena
    scene.add(mesh)

    # Salve a imagem do quadro atual
    image = scene.screenshot(size=(800, 600))
    images.append(image)

# Crie uma malha 3D simples (por exemplo, um cubo)
mesh = vedo.Cube()

# Configure os parâmetros da animação
angle_increment = 1  # Ângulo de rotação a ser incrementado a cada quadro
total_frames = 360  # Número total de quadros

# Crie a cena
scene = vedo.show(mesh, interactive=False)

# Lista para armazenar os quadros da animação
images = []

# Crie a animação
animation = FuncAnimation(plt.figure(), update_frame, frames=total_frames, repeat=False)

# Renderize a animação
plt.show()

# Salve a animação como um arquivo de vídeo
# imageio.mimsave('rotating_cube.gif', images, fps=30)