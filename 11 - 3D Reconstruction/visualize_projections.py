# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 20:35:55 2023

@author: HEartLab
"""

import scipy.io as sio
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import numpy as np
import os

root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
root.update()

filename = tk.filedialog.askopenfilename(title='Open projection file',
                                          filetypes=(('MATLAB files', '*.mat'),
                                                      ('Numpy arrays', '*.npy'),
                                                      ('All files', '*.*')))
root.destroy()

matfile = sio.loadmat(filename)

    
# Criação de um dicionário para armazenar os resultados finais
vertices = matfile["geometry_HR"]["vertices"][0,0]
# vertices[:,2] *= -1
faces = matfile["geometry_HR"]["faces"][0, 0] -1
#%%
signame = ''
def load_item():
    global signame
    signame =  dropdown.get()
    root.destroy()
    
    
# Create the main window
root = tk.Tk()
root.title("Selecting data")

# Create a list of items
items = set(matfile.keys())
for key in ["__globals__",
            "__header__", 
            "__version__",
            "dataA",
            "dataA_path",
            "dataB",
            "dataB_path",
            "dataC",
            "dataC_path",
            "geometry_HR", 
            "geometry_20000",
            "geometry_10000",
            "geometry_2500", 
            "geometry_1200",
            "projected_mask",
            "projected_vertices"]:
    
    items.discard(key)
    
items = list(items) + ["blank"]


# Create a StringVar to store the selected item
selected_item = tk.StringVar()

# Create a Label
label = tk.Label(root, text="Select the desired signal or map:")

# Create a Dropdown Menu
dropdown = ttk.Combobox(root, textvariable=selected_item, values=items)

# Create a Load button
load_button = tk.Button(root, text="Load", command=load_item)

# Pack the widgets
label.pack(pady=10)
dropdown.pack(pady=5)
load_button.pack(pady=10)

root.mainloop()

if signame == "blank":
    signals = np.zeros((vertices.shape[0], 3)) + 191

else:
    signals = matfile[signame] #[matfile["geometry_2500"]["HR_idx"][0,0].ravel() -1]


# import open3d as o3d

# mesh = o3d.geometry.TriangleMesh()
# mesh.vertices = o3d.utility.Vector3dVector(vertices)
# mesh.triangles = o3d.utility.Vector3iVector(faces)
# mesh.vertex_colors = o3d.utility.Vector3dVector(signals)
import vedo
import matplotlib.pyplot as plt

scene = vedo.Plotter().parallel_projection(True)

colormap = "gray"

if signals.shape[1] in (3, 4):
    mesh =  vedo.Mesh([vertices, faces],).lighting(None)  
    
    if signals.max() <= 1:

        signals *= 255
        signals = np.round(signals).astype(int)
    else:
        signals = np.round(signals).astype(int)
    
    if signals.shape[1] == 3:
        signals = np.hstack([signals, np.zeros([signals.shape[0], 1]) + 255]).astype(int)
        
    mesh.pointcolors = signals
    
    
    
elif signals.shape[1] > 0:
    # parula = list(zip(np.arange(parula_colors.shape[0]), *parula_colors.T))
    from parula import parula_colors
    # vmin = np.percentile(signals[:,0][signals[:,0] > 0], 1)
    vmin = 0
    vmax = np.percentile(signals[:,0][signals[:,0] > 0], 99)
    if colormap == 'parula':
        cmap = parula_colors[::-2]
    
    else:
        cmap = plt.cm.get_cmap(colormap, 100)

        # Obtém as cores do colormap em um array numpy
        cmap = cmap(np.linspace(0, 1, 100)) 
        cmap[0] = [0.75, 0.75, 0.75, 1]
        
        
    mesh =  vedo.Mesh(
        [vertices, faces]
        ).cmap(
            input_cmap=cmap[:,:3], 
            input_array=signals[:,0],
            vmin=vmin, vmax=vmax
            ).lighting(
                "plastic"
                ).add_scalarbar(size=(100, 1000), font_size=30)


for mea in matfile["MEAS_HR_IDX"]:
    scene += vedo.Points(vertices[matfile["MEAS_HR_IDX"][mea][0,0][0] - 1].T.tolist(),r=30, c=tuple(np.random.rand(3)))
    
scene += mesh

scene.show(interactive=True)
