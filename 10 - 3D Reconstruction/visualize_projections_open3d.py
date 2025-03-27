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

root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
root.update()

output_folder = tk.filedialog.askdirectory(title='Select a folder where the video would be saved')
root.destroy()



#%%
import vedo
import matplotlib.pyplot as plt
# Crie um colormap "jet" com 100 cores

## Visualize results
cam = dict(pos=(-290.5, -157.2, -270.0),
           focalPoint=(334.1, 384.7, 223.5),
           viewup=(0.3394, 0.3888, -0.8565),
           distance=963.0,
           clippingRange=(465.9, 1591),
)



scene = vedo.Plotter(interactive=True).parallel_projection(True)

colormap = "gray"

# scene += vedo.Mesh((local_sphere_vertices, local_triangles), c='w').lw(2).lc('black')
# scene += vedo.Points(sphere_HR_vertices[HR_known], c='b', r=13)
# scene += vedo.Points(sphere_HR_vertices[to_interp],r=12, c='g')
# scene += vedo.Points(sphere_HR_vertices[to_project],r=10, c='r')
# mesh = vedo.Mesh([vertices, faces],)
# mesh.cellcolors = np.hstack([signals[faces[:,0]], np.ones((faces.shape[0], 1))])

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
        
        
    mesh =  vedo.Mesh([vertices, faces]).cmap(input_cmap=cmap[:,:3], 
                                              input_array=signals[:,0],
                                              vmin=vmin, vmax=vmax
                                              ).lighting("plastic", specular_color=[1, 1, 1],
                                                         ).add_scalarbar(size=(100, 1000), font_size=30)

    
                                         
# Print information about the cell arrays of the mesh and their shape


# Display the mesh with the assigned colors and the docstring

scene += mesh

scene.show()
cam = dict(
    position = scene.camera.GetPosition(),
    focal_point = scene.camera.GetFocalPoint(),
    viewup = scene.camera.GetViewUp(),
    distance = scene.camera.GetDistance(),
    clipping_range = scene.camera.GetClippingRange(),
)

scene.close()



#%%

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

num_frames = 500

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
        
        
    mesh =  vedo.Mesh([vertices, faces]).cmap(input_cmap=cmap[:,:3], 
                                              input_array=signals[:,0],
                                              vmin=vmin, vmax=vmax
                                              ).lighting("plastic"
                                                          ).add_scalarbar(size=(100, 1000), font_size=30)


for mea in ["MEA1", "MEA2", "MEA3"]:
    scene += vedo.Points(vertices[matfile["MEAS_HR_IDX"][mea][0,0][0]].T.tolist(),r=30, c=tuple(np.random.rand(3)))
    
scene += mesh


scene.show(camera=cam, interactive=False)


# Laço para renderizar cada quadro
for frame in range(num_frames):
    # Defina a posição da câmera para cada quadro (ajuste conforme necessário)
    
    scene.camera.Azimuth(360/num_frames)  # Rotacionar 1 grau em torno do objeto

    # Renderize o quadro atual
    scene.screenshot(f"{output_folder}/frame_{frame:04d}.png")
    print("Frame: ", frame)
# Feche a cena Vedo
                                         
scene.interactive()

from moviepy.editor import ImageSequenceClip

# Specify the folder containing your images
image_folder = output_folder

# Specify the frame rate (frames per second)
frame_rate = 30

output_file = f"{output_folder}/{os.path.basename(output_folder)}.mp4"
if os.path.exists(output_file):
    os.remove(output_file)
# Load the images from the folder and create a video
video_clip = ImageSequenceClip(image_folder, fps=frame_rate)

# Specify the output file name

# Write the video to a file
video_clip.write_videofile(output_file, codec='libx264', fps=frame_rate)



scene.close()

#%%

# from scipy.signal import butter, sosfiltfilt, find_peaks

# sos = butter(8, [0.5, 30], btype='bandpass', fs=500, output='sos')
# f_signals = np.zeros_like(signals)
# f_signals = sosfiltfilt(sos, signals)
# f_signals -= np.percentile(f_signals, 1, axis=1, keepdims=True)
# f_signals /= np.percentile(f_signals, 99, axis=1, keepdims=True)

# if not os.path.exists(output_folder):
#     os.makedirs(output_folder)

# num_frames = 2000

# scene = vedo.Plotter().parallel_projection(True)

# colormap = "jet"

    
# if f_signals.shape[1] > 0:
#     # parula = list(zip(np.arange(parula_colors.shape[0]), *parula_colors.T))
#     from parula import parula_colors
#     vmin, vmax = np.percentile(f_signals, [1, 99])
     
#     if colormap == 'parula':
#         cmap = parula_colors[::-2]
    
#     else:
#         cmap = plt.cm.get_cmap(colormap, 100)

#         # Obtém as cores do colormap em um array numpy
#         cmap = cmap(np.linspace(0, 1, 100)) 
#         # cmap[0] = [0.75, 0.75, 0.75, 1]
        
        
#     mesh = vedo.Mesh([vertices, faces]).cmap(input_cmap=cmap[:,:3], 
#                                               input_array=f_signals[:,0],
#                                               vmin=vmin, vmax=vmax
#                                               ).lighting("plastic"
#                                                           ).add_scalarbar(size=(100, 1000), font_size=30)



# scene += mesh
# scene.show(camera=cam, interactive=False)


# # Laço para renderizar cada quadro
# for frame in range(num_frames):
#     # Defina a posição da câmera para cada quadro (ajuste conforme necessário)
    
#     scene.camera.Azimuth(360/num_frames)  # Rotacionar 1 grau em torno do objeto
#     # mesh.pointdata["scalars"] = f_signals[:, frame]
#     mesh.cmap(input_cmap=cmap[:,:3], input_array=f_signals[:,frame + 1000], vmin=vmin, vmax=vmax)
#     scene.render()
#     # Renderize o quadro atual
#     scene.screenshot(f"{output_folder}/frame_{frame:04d}.png")
#     print("Frame: ", frame)
# # Feche a cena Vedo
                                         
# scene.interactive()

# from moviepy.editor import ImageSequenceClip

# # Specify the folder containing your images
# image_folder = output_folder

# # Specify the frame rate (frames per second)
# frame_rate = 30

# output_file = f"{output_folder}/{os.path.basename(output_folder)}.mp4"
# if os.path.exists(output_file):
#     os.remove(output_file)
# # Load the images from the folder and create a video
# video_clip = ImageSequenceClip(image_folder, fps=frame_rate)

# # Specify the output file name

# # Write the video to a file
# video_clip.write_videofile(output_file, codec='libx264', fps=frame_rate)

# scene.close()

#%%


# from scipy.signal import butter, sosfiltfilt, find_peaks

# def get_connectivity(triangles, autoindex=True):
#     """Get 1st and 2nd order connectivity matrix from triangles."""
#     size = triangles.max() + 1
#     connectivity = np.zeros((size, size), bool)
#     for i, j, k in triangles:
#         connectivity[i,j] = True
#         connectivity[i,k] = True
#         connectivity[j,i] = True
#         connectivity[k,i] = True
#         connectivity[j,k] = True
#         connectivity[k,j] = True
    
    
#     if autoindex:
#         connectivity[np.arange(size), np.arange(size)] = True
#     else:
#         connectivity[np.arange(size), np.arange(size)] = False

#     return connectivity

# def mesh_blur(signals, faces, it=1):
#     connectivity = get_connectivity(faces) 
#     outsignals = np.zeros_like(signals)
    
#     for i in range(signals.shape[0]):
#         outsignals[i] = signals[connectivity[i]].mean(axis=0)
        
#     return outsignals

# signals_blur = mesh_blur(signals, faces, it=10)

# sos = butter(8, [0.5, 30], btype='bandpass', fs=500, output='sos')
# f_signals = sosfiltfilt(sos, signals_blur)[:,1000:-1000]
# f_signals -= np.percentile(f_signals, 1, axis=1, keepdims=True)
# f_signals /= np.percentile(f_signals, 99, axis=1, keepdims=True)


# peaks_idx, _ = find_peaks(1 - f_signals.ravel(), height=0.8, distance=500 * 1/120)
# peaks = np.zeros_like(f_signals.ravel())
# peaks[peaks_idx] = 1

# peaks = peaks.reshape(*f_signals.shape)

# if not os.path.exists(output_folder):
#     os.makedirs(output_folder)

# num_frames = 2000

# scene = vedo.Plotter().parallel_projection(True)

# colormap = "jet"

    
# if f_signals.shape[1] > 0:
#     # parula = list(zip(np.arange(parula_colors.shape[0]), *parula_colors.T))
#     from parula import parula_colors

#     cmap = plt.cm.get_cmap(colormap, 100)

#     # Obtém as cores do colormap em um array numpy
#     cmap = cmap(np.linspace(0, 1, 100)) 
#     # cmap[0] = [0.75, 0.75, 0.75, 1]
        
        
#     mesh = vedo.Mesh([vertices, faces]).cmap(input_cmap=cmap[:,:3], 
#                                               input_array=peaks[:,0],
#                                               vmin=0, vmax=1
#                                               ).lighting(None)



# scene += mesh
# scene.show(camera=cam, interactive=False)


# # Laço para renderizar cada quadro
# for frame in range(num_frames):
#     # Defina a posição da câmera para cada quadro (ajuste conforme necessário)
    
#     scene.camera.Azimuth(360/num_frames)  # Rotacionar 1 grau em torno do objeto
#     # mesh.pointdata["scalars"] = f_signals[:, frame]
#     mesh.cmap(input_cmap=cmap[:,:3], input_array=peaks[:,frame], vmin=0, vmax=1)
#     scene.render()
#     # Renderize o quadro atual
#     scene.screenshot(f"{output_folder}/frame_{frame:04d}.png")
#     print("Frame: ", frame)
# # Feche a cena Vedo
                                         
# scene.interactive()

# from moviepy.editor import ImageSequenceClip

# # Specify the folder containing your images
# image_folder = output_folder

# # Specify the frame rate (frames per second)
# frame_rate = 30

# output_file = f"{output_folder}/{os.path.basename(output_folder)}.mp4"
# if os.path.exists(output_file):
#     os.remove(output_file)
# # Load the images from the folder and create a video
# video_clip = ImageSequenceClip(image_folder, fps=frame_rate)

# # Specify the output file name

# # Write the video to a file
# video_clip.write_videofile(output_file, codec='libx264', fps=frame_rate)



# scene.close()