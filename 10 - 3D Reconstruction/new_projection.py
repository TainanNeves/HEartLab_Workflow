# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 17:53:07 2023

@author: HEartLab
"""
#%%
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
root = tk.Tk()
tk.Label(root, text="New Data Projection is running").pack()

import scipy.io as sio
from tkinter import filedialog, messagebox
import numpy as np
import h5py
from PIL import Image
from tkinter import ttk
import matplotlib.pyplot as plt
from cv2 import fillPoly, findContours
# from scipy.interpolate import griddata
from open_strampix import Open_StreamPix
from projection_tools import project_data
from mesh_processing import get_connectivity, erodemorph3d, dilatemorph3d

#%%

window = tk.Toplevel(root)
window.withdraw()
window.attributes("-topmost", True)
window.update()

matfilename = tk.filedialog.askopenfilename(title='Open projection file',
                                          filetypes=(('MATLAB files', '*.mat'),
                                                      ('Numpy arrays', '*.npy'),
                                                      ('All files', '*.*')))
window.destroy()

if matfilename == '':
    exit()

matfile = sio.loadmat(matfilename)
    
# Criação de um dicionário para armazenar os resultados finais
HR_vertices = matfile["geometry_HR"]["vertices"][0,0]
HR_faces = matfile["geometry_HR"]["faces"][0, 0] -1
conn = get_connectivity(HR_faces, autoindex=False)

#%%

#Esta função lida com o carregamento de dados de arquivos e realiza algumas verificações para garantir que os dados sejam utilizáveis 
def load_file(block_num):
    # Execute your generic function here
    if dropdowns[block_num].get() == "---":
        data = np.asarray(files[block_num], np.float64)
    else:
        data = np.asarray(files[block_num][dropdowns[block_num].get()])
    
    if len(data.shape) == 2:
        data = data[np.newaxis]
    elif len(data.shape) == 1:
        print('Invalid data')
        return
    
    datas[block_num] = data.astype(np.float64)
    turners[block_num]["state"] = "normal"
    flippers[block_num]["state"] = "normal"
    rolldims[block_num]["state"] = "normal"


def update_image(block_num):
    # Update the image to white
    axes[block_num].imshow(datas[block_num][0], cmap='gray')
    canvas[block_num].draw()

def open_signal(block_num):
    filename = tk.filedialog.askopenfilename(title="Select the signals or map for the " + camnames[block_num])
    if filename:
        entry[block_num]["state"] = 'normal'
        entry[block_num].delete(0, tk.END)
        entry[block_num].insert(0, filename)
        entry[block_num]["state"] = 'disabled'
        
    paths[block_num] = entry[block_num].get()
    if paths[block_num].endswith(".mat"):
        try:
            signalfile = sio.loadmat(paths[block_num])
        except NotImplementedError:
            signalfile = dict(h5py.File(paths[block_num]))
            for k in signalfile.keys():
                signalfile[k] = signalfile[k][()]
            
        dropdowns[block_num]["values"] = [k for k in signalfile.keys()]
        dropdowns[block_num]["state"] = "normal"
        dropdowns[block_num].set(dropdowns[block_num]["values"][0])
        
    elif paths[block_num].endswith(".npy"):
        signalfile = np.load(paths[block_num], allow_pickle=True)
        if isinstance(signalfile[()], dict):
            signalfile = signalfile[()]
            dropdowns[block_num]["values"] = [k for k in signalfile.keys()]
            dropdowns[block_num]["state"] = "normal"
            dropdowns[block_num].set(dropdowns[block_num]["values"][0])
        
        else:
            dropdowns[block_num]["values"] = ["---",]
            dropdowns[block_num]["state"] = "disabled"
            
    elif paths[block_num].endswith((".bmp", ".dib", ".jpeg", ".jpg", ".jpe", ".jp2", 
                        ".png", ".webp", ".pbm", ".pgm", ".ppm", ".pxm", 
                        ".pnm", ".pfm", ".sr", ".ras", ".tiff", ".tif", 
                        ".exr", ".hdr", ".pic")):
        signalfile = np.array(Image.open(paths[block_num])).transpose(2, 0, 1)
        dropdowns[block_num]["values"] = ["---",]
        dropdowns[block_num]["state"] = "disabled"
    
    
    elif paths[block_num].endswith(".seq"):
        signalfile = Open_StreamPix(paths[block_num], maxframes=10)[1].transpose(2, 0, 1)
        dropdowns[block_num]["values"] = ["---",]
        dropdowns[block_num]["state"] = "disabled"
        
    
    files[block_num] = signalfile
    
    
def turn_image(block_num):
    datas[block_num] = datas[block_num].transpose(0, 2, 1)[:,:,::-1]
    
def flip_image(block_num):
    datas[block_num] = datas[block_num][:,::-1,:]
   
def rolldim_image(block_num):
    datas[block_num] = datas[block_num].transpose(2, 0, 1)
    
is_phase = False
no_data_value = float("nan")
quit_count = 0
def confirm_selection():
    global is_phase, no_data_value, quit_count
    if not all([d is None for d in datas]):
        is_phase = phase_checkbox_var.get() != 0
        no_data_value = float(no_data_entry.get()) if no_data_entry.get() != "" else float("NaN")
        window.destroy()
        root.quit()
    else:
        quit_count += 1
        print("Please load all views before confirm")
        if quit_count >= 10:
            quit()
   
 
# Create the main window
window = tk.Toplevel(root)
window.title("Load Signal Files")
window.attributes('-topmost', 1)
window.protocol("WM_DELETE_WINDOW", confirm_selection)
# window.geometry("800x300")

# Create three blocks
camnames = ["camera C", "camera A", "camera B"]
files = [None, None, None]
datas = [None, None, None]
paths = [ None, None, None]

blocks = []
entry = []
buttons = []
axes = []
canvas = []
selected_options = []
dropdowns = []
turners = []
flippers = []
rolldims = []

for i in range(3):
    # Create a frame for each block
    block_frame = tk.Frame(window, padx=10, pady=10)
    block_frame.pack(side=tk.LEFT)
    
    
    text_label = tk.Label(block_frame, text=camnames[i])
    text_label.pack()
    
    # Create an entry for displaying the selected file name
    entry.append(tk.Entry(block_frame, width=40, justify='right',state='disabled'))
    entry[i].pack()

    # Create a button to open a file dialog
    buttons.append(tk.Button(block_frame, text="Select Signal or Map File", command=lambda i=i: open_signal(i)))
    buttons[i].pack()
    
    # Create a dropdown menu for options A, B, C
    selected_options.append(tk.StringVar(value="---"))
    dropdowns.append(ttk.Combobox(block_frame, textvariable=selected_options[i], values=["---",], state="disabled"))
    dropdowns[i].pack(pady=5)
    
    # Create a button to load the file and update the image
    load_button = tk.Button(block_frame, text="Load", command=lambda i=i: [load_file(i), update_image(i)])
    load_button.pack(pady=5)

    img_frame = tk.Frame(block_frame, padx=10, pady=10)
    img_frame.pack()
    
    # Create a figure and an axis for displaying the image
    fig, ax = plt.subplots(figsize=(3,5))
    # fig.tight_layout(pad=0)
    # ax.spines['left'].set_position('center')
    # ax.spines['bottom'].set_position('center')
    fig.set_facecolor('#d9d9d9')
    ax.set_xticks([])
    ax.set_yticks([])
    axes.append(ax)
    canvas.append(FigureCanvasTkAgg(fig, master=img_frame))
    canvas[i].get_tk_widget().grid(row=0, column=0, columnspan=3)
    
    turners.append(tk.Button(img_frame, text = "⭮", command = lambda i=i: [turn_image(i), update_image(i)], state="disabled", width=10))
    turners[i].grid(row=1, column=0)
    flippers.append(tk.Button(img_frame, text = " ⏛ ", command = lambda i=i: [flip_image(i), update_image(i)], state="disabled", width=10))
    flippers[i].grid(row=1, column=1)
    rolldims.append(tk.Button(img_frame, text = " 123 -> 312 ", command = lambda i=i: [rolldim_image(i), update_image(i)], state="disabled", width=10))
    rolldims[i].grid(row=1, column=2)
    

options_frame = tk.Frame(window, padx=10, pady=10)
options_frame.pack(side=tk.LEFT)

phase_checkbox_var = tk.StringVar(value="0")
phase_checkbox = tk.Checkbutton(options_frame, text="Phase map", variable=phase_checkbox_var, anchor="e")
phase_checkbox.pack()

no_data_entry_frame = tk.Frame(options_frame)
no_data_entry_frame.pack()
no_data_label = tk.Label(no_data_entry_frame, text = "\"No Data\" value:")
no_data_label.pack(side=tk.LEFT)
no_data_entry = tk.Entry(no_data_entry_frame, width=4)
no_data_entry.pack(side=tk.RIGHT)
no_data_entry.insert(0, "NaN")

confirm_button = tk.Button(options_frame, text="Confirm", command=confirm_selection)
confirm_button.pack(side=tk.BOTTOM, pady=100)

window.attributes('-topmost', 0)

# Start the Tkinter event loop
root.mainloop()


#%%

for i in range(3):
    if not datas[i] is None:
        n, c = np.unique(datas[i], return_counts=True)
        if c.max() >= c.sum() / 4:
            val = n[c.argmax()]
        else:
            a, b = np.histogram(datas[i], bins=1024)
            val = c[(c >= b[a.argmax()]) & (c < b[a.argmax() + 1])].mean()

        empty = np.ones_like(datas[i]) * val

for i in range(3):
    if datas[i] is None:
        datas[i] = empty
#%%

# Criação de um dicionário para armazenar os resultados finais
# conn = get_connectivity(HR_faces, autoindex=True)
projected_vertices = matfile["projected_vertices"].copy()
projections_mask = matfile["projected_mask"].copy() > 0.5

projected_vertices = projected_vertices.astype(np.int32)
projected_vertices[:,:,0][projected_vertices[:,:,0] >= datas[0].shape[2]] = datas[0].shape[2] - 1
projected_vertices[:,:,1][projected_vertices[:,:,1] >= datas[0].shape[1]] = datas[0].shape[1] - 1 
projected_vertices[projected_vertices < 0] = 0

#%%

contours = []
print("Computing projection ROIs.", end="")
for cam, img in enumerate(datas):     
    triangles_mask = np.zeros_like(HR_faces, dtype = bool)
    for i in np.where(projections_mask[cam])[0]:
        triangles_mask |= HR_faces == i
    
    triangles_mask = np.all(triangles_mask, axis=1)
    proj_triangles = HR_faces[triangles_mask]
    for i, n in enumerate(np.unique(proj_triangles)):
        proj_triangles[proj_triangles == n] = i
        
    
    mask = np.zeros(img.shape[1:3], dtype=np.uint8)
    mask[projected_vertices[cam,:,1], projected_vertices[cam,:,0]] = 1
    for tri in HR_faces:
        triangle = np.array(projected_vertices[cam][tri,:])
        if np.any(triangle.ptp(axis=0) > 5):
            print("ops")
            break
        
        if all(triangle[0] == triangle[1]) or all(triangle[0] == triangle[2]) or all(triangle[1] == triangle[2]):
            continue
        
        else:
            triangle = ((triangle - triangle.mean(axis=0, keepdims=True))*2 + triangle.mean(axis=0, keepdims=True)).astype(np.int32)
            triangle = np.vstack([triangle, triangle[(0,),]])
            fillPoly(mask, pts = [triangle], color=1)
    
    c = findContours(mask, 0, 1)[0][0][:,0,:]
    
    contours.append(c)
    print(".", end="")
print("")

#%%
plt.close("all")

fig, ax = plt.subplots(nrows=1, ncols=3)
fig.suptitle("Projected geometry over new maps")

for i in range(3):
    ax[i].set_title(camnames[i])
    ax[i].imshow(datas[i][0], cmap='gray')
    ax[i].plot(contours[i][:,0], contours[i][:,1], 'r', linewidth=5)

plt.show()
plt.pause(4)
plt.close("all")


#%%
# Função para subamostragem 2D dos dados usando coordenadas xi
def subsample2D(data, xi):
    
    n, h, w = data.shape
    data = data.copy()
    xi = np.round(xi.copy()).astype(int)
    dataout = np.zeros((xi.shape[0], n), dtype=data.dtype)
    
    # Criação de uma máscara para coordenadas xi fora dos limites da imagem
    mask = (xi[:,0] < 0) | (xi[:,1] < 0) | (xi[:,0] >= w) | (xi[:,1] >= h)
    
    # Preenchimento de locais fora da imagem com NaN
    if dataout.dtype in [np.float16, np.float32, np.float64, float]:
        dataout[mask, :] = np.nan
    elif dataout.dtype in [np.uint8, np.uint16, np.uint32, np.uint64]:
        dataout[mask, :] = 0
    
    else:
        dataout[mask, :] = -1
    
        
    # Preenchimento de locais dentro da imagem com os dados correspondentes
    dataout[~mask] = data[:, xi[~mask,1], 
                             xi[~mask,0]].T
    
        
    return dataout


#%%

for i in range(3):
    if datas[i].ndim == 2:
        datas[i] = datas[i][np.newaxis,:,:]
    
t = datas[0].shape[0]
weights = matfile["projection_weights"]
signals = project_data(datas, projected_vertices, projections_mask, HR_vertices, weights, no_data_value=no_data_value, phase=is_phase)

final_mask = signals[:,0] != no_data_value
final_mask = erodemorph3d(final_mask, conn, r=10)
final_mask = dilatemorph3d(final_mask, conn, r=20)
final_mask = erodemorph3d(final_mask, conn, r=10)

no_data = signals[:,0] == no_data_value

# if no_data.sum() > 0 and no_data.sum()/no_data.size < 0.4:
for _ in range(20):
    if no_data.sum() == 0:
        break
    for i in np.where(no_data * conn[~no_data].any(axis=0))[0]:
        signals[i] = signals[conn[i] * ~no_data].mean()
    no_data = signals[:,0] == no_data_value

signals[~final_mask] = float("NaN")


#%%

class NameEntryWindow:
    def __init__(self, window):
        self.window = window
        self.window.title("Projection")
        self.root = self.window.winfo_toplevel()

        self.name_var = tk.StringVar()
        self.exit_flag = tk.BooleanVar(value=True)

        self.label = tk.Label(window, text="Enter the new projection name:")
        self.entry = tk.Entry(window, textvariable=self.name_var)
        self.create_button = tk.Button(window, text="Create", command=self.create_name)
        self.exit_button = tk.Button(window, text="Exit", command=self.exit_window)

        self.label.grid(row=0, column=0)
        self.entry.grid(row=1, column=0)
        self.create_button.grid(row=1, column=1)
        self.exit_button.grid(row=2, column=1)
        self.name = "new_proj"
        
    def create_name(self):
        self.name = self.name_var.get()
        self.window.destroy()  # Close the window
        self.root.quit()

    def exit_window(self):
        self.window.destroy()  # Close the window
        self.exit_flag.set(False)
        self.root.quit()



#%%
window = tk.Toplevel(root)
app = NameEntryWindow(window)
root.mainloop()
mapname = app.name
root.destroy()

#%%

matfile[mapname] = signals
sio.savemat(matfilename, matfile)

print(f"New projected data saved in the variable {mapname} inside the same file ({matfilename})")

#%%

import vedo
from parula import parula_colors
# vmin = np.percentile(signals[:,0][signals[:,0] > 0], 1)
scene = vedo.Plotter(interactive=True)
if is_phase:
    vmin = -np.pi
    vmax = np.pi
else:
    vmin = np.percentile(signals[:,0][signals[:,0] > 0], 1)
    vmax = np.percentile(signals[:,0][signals[:,0] > 0], 99)
cmap = parula_colors
cmap[0] = [0.75, 0.75, 0.75]
    
    
mesh =  vedo.Mesh([HR_vertices, HR_faces]).cmap(input_cmap=cmap[:,:3], 
                                            input_array=signals[:,0],
                                            vmin=vmin, vmax=vmax,
                                            ).lighting("plastic"
                                                        ).add_scalarbar(size=(100, 1000), font_size=30)


scene += mesh
scene.show().interactive().close()

