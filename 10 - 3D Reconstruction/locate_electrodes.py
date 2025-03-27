# -*- coding: utf-8 -*-
"""
Created on Wed May 31 11:25:52 2023

@author: HEartLab
"""


import scipy.io as sio

import numpy as np
import matplotlib
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog, messagebox
matplotlib.use('tkagg')


#%%
def find_HR_idx(HR_vertices, LR_vertices):
    dist = HR_vertices[:,np.newaxis,:] - LR_vertices[np.newaxis,:,:]
    dist = np.linalg.norm(dist, axis=-1)
    return np.argmin(dist, axis=0)

#%%


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
HR_vertices = matfile["geometry_HR"]["vertices"][0,0]
HR_faces = matfile["geometry_HR"]["faces"][0, 0] - 1
    
datas = [matfile["dataC"][0], matfile["dataA"][0], matfile["dataB"][0]]
projected_vertices = matfile["projected_vertices"]
projections_mask = matfile["projected_mask"]



#%%
class NameEntryWindow:
    def __init__(self, root):
        self.root = root
        self.root.title("Add New MEA")

        self.name_var = tk.StringVar()
        self.exit_flag = tk.BooleanVar(value=True)
        self.validate_cmd = root.register(self.validate_entry)
        self.label = tk.Label(root, text="Enter the MEA name:")
        self.entry = tk.Entry(root, textvariable=self.name_var)
        self.name_var.trace_add(mode="write", callback=self.validate_entry)
        self.create_button = tk.Button(root, text="Create", command=self.create_name, state='disabled')
        self.exit_button = tk.Button(root, text="Exit", command=self.exit_window)
        self.label.grid(row=0, column=0)
        self.entry.grid(row=1, column=0)
        self.create_button.grid(row=1, column=1)
        self.exit_button.grid(row=2, column=1)
        self.MEA = ""
        
    def create_name(self):
        self.MEA = self.name_var.get().replace(" ", "_")
        
        self.root.destroy()  # Close the window

    def exit_window(self):
        self.root.destroy()  # Close the window
        self.exit_flag.set(False)  # Set exit flag to False
    
    def validate_entry(self, *args):
        # Enable or disable the button based on the entry text
        name = self.name_var.get()

        if name == "":
            self.create_button.config(state=tk.DISABLED)

        elif name.replace(" ", "")[0].isnumeric():
            self.create_button.config(state=tk.DISABLED)
        
        elif not name.replace(" ", "").isalnum():
            self.create_button.config(state=tk.DISABLED)

        else:
            self.create_button.config(state=tk.NORMAL)

meas = {}


keeprunning = True
while keeprunning:

    root = tk.Tk()
    app = NameEntryWindow(root)
    root.mainloop()
    keeprunning = app.exit_flag.get()
    
    if not keeprunning:
        break
    
    meaname = app.MEA
    
    electrode_indices = np.zeros((16,3), int) -1
    for cam, img in enumerate(datas):
        fig, ax = plt.subplots()
        ax.imshow(img, cmap='gray')
        ax.set_title("Locate the electrodes for " + meaname)
        sct = ax.scatter([],[], s = 200, c='r',marker='x',zorder=10, alpha=0.75)
        
        points = []
        def onclick(event):
            # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            #       ('double' if event.dblclick else 'single', event.buvtton,
            #         event.x, event.y, event.xdata, event.ydata))
            global points
            global fig
            global ax
            
            if event.dblclick and event.button==1:
                points.append([event.xdata, event.ydata])
                arrpoints = np.asarray(points)
                sct.set_offsets(arrpoints)
                fig.canvas.draw()
        
            if event.button==3:
                arrpoints = np.asarray(points)
                dist = arrpoints - np.array([[event.xdata, event.ydata]])
                dist = np.linalg.norm(dist, axis=1)
                if dist.min() <= 20:
                    arrpoints = arrpoints[dist >= 20]
                    points.clear()
                    points.extend(arrpoints.tolist())
                    arrpoints = np.asarray(points)
                    if arrpoints.shape[0] == 0:
                        arrpoints = np.empty((0,2))
                    sct.set_offsets(arrpoints)
                    fig.canvas.draw()
                    
        # fig_open = True      
        # def fig_closed(event):
        #     global fig_open
        #     fig_open = False
        
        cid1 = fig.canvas.mpl_connect('button_press_event', onclick)
        # cid2 = fig.canvas.mpl_connect('close_event', fig_closed)
        
        
        
        plt.show(block=True)

        
        
        arrpoints = np.asarray(points).reshape(-1,2)
        arrpoints = arrpoints[arrpoints[:,0] != None]
        
        #===========================================================
        if arrpoints.shape[0] == 0:
            continue
        
        fig, ax = plt.subplots()
        ax.imshow(img, cmap='gray')
        ax.set_title("Select the indices of the electrodes for " + meaname)
        ax.scatter(arrpoints[:,0],arrpoints[:,1], s = 200, c='r',marker='x',zorder=10, alpha=0.75,  picker=True)
        sct = ax.scatter([], [], s = 200, c='g',marker='x',zorder=10, alpha=1)
        
        
        ind = 0
        counter = 0
        indices = np.zeros(arrpoints.shape[0], int) -1
        annotations = []
        for i in range(arrpoints.shape[0]):
            annotations.append(ax.annotate(text = "", 
                                           xy = arrpoints[i],
                                           xytext = np.array([-5,10]),
                                           textcoords = 'offset points',
                                           color='red', size=20))    
        def pickpoint(event):
            global ind
            global ax
            global fig
            global counter
            if isinstance(event.artist, matplotlib.collections.PathCollection):         
                # for att in dir(thisline): print(att)
                ind = event.ind[0]      
                sct.set_offsets(arrpoints[ind,:])
        
                indices[ind] = counter
                counter %= 16
                counter += 1
                
                annotations[ind].set_text(str(counter))
                fig.canvas.draw()
                
        def undofunction(event):
            global counter, ind, indices
            
            
            annotations[ind].set_text("")

            indices[ind] = -1
            
            counter -= 1
            if counter < 0:
                counter = 15
            
            sct.set_offsets([])
            fig.canvas.draw()

    
        
        def resetfunction(event):
            global counter, indices
            
            counter = 0
            indices = np.zeros(arrpoints.shape[0], int) -1
            for annot in annotations:
                annot.set_text("")
            
            fig.canvas.draw()
    
            
            
        # fig_open = True
        # def fig_closed(event):
        #     global fig_open
        #     fig_open = False
        
        # cid3 = fig.canvas.mpl_connect('close_event', fig_closed)
        cid4 = fig.canvas.mpl_connect('pick_event', pickpoint)
        
        
        fig.subplots_adjust(bottom=0.2)
        axundo = fig.add_axes([0.7, 0.05, 0.08, 0.075])
        axreset = fig.add_axes([0.79, 0.05, 0.12, 0.075])
        bundo = Button(axundo, 'Undo')
        bundo.on_clicked(undofunction)
        breset = Button(axreset, 'Reset')
        breset.on_clicked(resetfunction)
    
        plt.show(block=True)
            
        proj_electrodes = np.zeros((16,2)) -1
        proj_electrodes[indices] = arrpoints
        
        electrodes_proj_ind = np.zeros(16, int) -1
        for i, coord in enumerate(proj_electrodes):
            if coord[0] != -1:
                electrodes_proj_ind[i] = np.argmin(np.linalg.norm(projected_vertices[cam][projections_mask[cam]>0.5]
                                                                  - coord[np.newaxis,:], axis=1))
            
        proj_inds = np.where(projections_mask[cam])[0]
        electrodes_ind = np.zeros(16, int) -1
        electrodes_ind[electrodes_proj_ind > -1] = proj_inds[electrodes_proj_ind[electrodes_proj_ind > -1]]
        
        electrode_indices[:,cam] = electrodes_ind
        
        
    electrodes_coords = np.zeros((16,3)) -1
    
    for i, inds in enumerate(electrode_indices):
        electrodes_coords[i] = np.mean(HR_vertices[inds[inds>-1]], axis=0)
        
    final_indices = np.zeros(16, int) -1
    for i, coord in enumerate(electrodes_coords):
        if coord[0] != -1:
            final_indices[i] = np.argmin(np.linalg.norm(HR_vertices - coord[np.newaxis,:], axis=1))
        
    meas[meaname] = final_indices + 1
    
matfile["MEAS_HR_IDX"] = meas

for key in list(matfile):
    if key.startswith("geometry_") and not key.endswith("HR"):
        matfile[f"MEAS_IDX_{key[9:]}"] = dict([(measname,
                                                find_HR_idx(matfile[key]['vertices'][0,0],
                                                            HR_vertices[meas[measname] - 1]) + 1)
                                               for measname in meas])


sio.savemat(filename, matfile)
print(f"MEA information saved in projection file ({filename})")

# #%%
# import vedo
# import matplotlib.pyplot as plt

# scene = vedo.Plotter().parallel_projection(True)
# signals = matfile["map3d_HR"]
# vmin, vmax = np.percentile(signals, [5,95])
# mesh =  vedo.Mesh(
#     [HR_vertices, HR_faces]
#     ).cmap(
#         input_cmap="gray", 
#         input_array=signals[:,0],
#         vmin=vmin, vmax=vmax
#         ).lighting(
#             "plastic"
#             ).add_scalarbar(size=(100, 1000), font_size=30)


# for mea in meas:
#     scene += vedo.Points(HR_vertices[matfile["MEAS_HR_IDX"][mea][0,0][0] - 1].tolist(),r=HR_vertices.ptp(axis=0).mean()/20, c=tuple(np.random.rand(3)))
    
# scene += mesh

# scene.show(interactive=True)