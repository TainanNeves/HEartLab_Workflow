#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 03:06:41 2022

@author: italo
"""
#%%
from matplotlib.widgets import Button, Slider
import matplotlib
from skimage.measure import label
import matplotlib.pyplot as plt
from tkinter import filedialog, messagebox
import scipy.signal as sig
import scipy.io as sio
import tkinter as tk
import numpy as np
import cv2
import os
matplotlib.use('tkagg')
from scipy import optimize
import open3d as o3d
from tqdm import tqdm
from edit_segmentation import manual_segmentation_edit
from scipy import signal
import datetime
#%%
# 



def lite_dilate3D(vol, iterations=1, cte=0):
    """
    Do morphologic operations 'dilate' in a memory optimized way
    """
    vol = (vol > 0).astype(np.uint8)
    vol = np.pad(vol, 1, mode='constant', constant_values=cte)
    for i in range(iterations): 
        vol[:, :, 1:-1] += vol[:, :, :-2] + vol[:, :, 2:]
        vol[:, 1:-1, :] += vol[:, :-2, :] + vol[:, 2:, :]
        vol[1:-1, :, :] += vol[:-2, :, :] + vol[2:, :, :]
        vol = (vol>0).astype(np.uint8)
        
    return vol[1:-1,1:-1,1:-1]

def lite_erode3D(vol, iterations=1):
    """
    Do morphologic operations 'erode' in a memory optimized way
    """
    vol = (vol > 0).astype(np.uint8)
    return 1 - lite_dilate3D(1 - vol, iterations, cte=1)

#%%
# Do mean blur 3D in a optimized way

def lite_blur(vol, k=3, cte=0):
    assert k%2==1, 'k should be an odd integer'
    h, w, l = vol.shape
    vol = (vol > 0).astype(np.uint16)
    p = int(np.floor(k/2))
    vol = np.pad(vol, p, mode='constant', constant_values=cte)
    vol_out = vol.copy()
    for i in range(0,p):
        vol_out[:, :, p:-p] += vol[:, :, i:l+i] + vol[:, :, k-1-i:l+k-1-i]
    vol = vol_out.copy()

    for i in range(0,p):
        vol_out[:, p:-p, :] += vol[:, i:w+i, :] + vol[:, k-1-i:w+k-1-i, :]
    vol = vol_out.copy()

    for i in range(0,p):
        vol_out[p:-p, :, :] += vol[i:h+i, :, :] + vol[k-1-i:h+k-1-i, :, :]
    vol = vol_out.copy()

    return vol[p:-p,p:-p,p:-p]
#%%
# Selects the biggest connected component in 2D
def biggest_component(img):
    img = ((img.copy()>0)*255).astype(np.uint8)
    
    
    _, labels = cv2.connectedComponents(img)
    n, counts = np.unique(labels, return_counts=True)
    counts = counts[n!=0]
    n = n[n!=0]
    bigest = n[counts.argsort()[-1]]
    img = ((labels!=bigest)*255).astype(np.uint8)
    
    _, labels = cv2.connectedComponents(img)
    n, counts = np.unique(labels, return_counts=True)
    counts = counts[n!=0]
    n = n[n!=0]
    bigest = n[counts.argsort()[-1]]
    
    return ((labels!=bigest)*255).astype(np.uint8)
#%%
# Selects the biggest connected component in 3D

def biggest_component3D(img):
    img = ((img.copy()>0)*255).astype(np.uint8)
    
    
    labels = label(img, background=0)
    n, counts = np.unique(labels, return_counts=True)
    counts = counts[n!=0]
    n = n[n!=0]
    bigest = n[counts.argsort()[-1]]
    img = ((labels!=bigest)*255).astype(np.uint8)
    
    labels = label(img, background=0)
    n, counts = np.unique(labels, return_counts=True)
    counts = counts[n!=0]
    n = n[n!=0]
    bigest = n[counts.argsort()[-1]]
    
    return ((labels!=bigest)*255).astype(np.uint8)

kernel = np.array([[0, 1, 1, 1, 0],
                   [1, 1, 1, 1, 1],
                   [1, 1, 1, 1, 1],
                   [1, 1, 1, 1, 1],
                   [0, 1, 1, 1, 0]], dtype=np.uint8)

#%%

class draw_polygon:
    def __init__(self, fig, ax, linecolor):
        self.fig = fig
        self.ax = ax
        self.line = ax.plot([],[], color=linecolor, zorder=10, alpha=0.75)[0]
        self.dots = ax.plot([],[], color=linecolor, marker='o', zorder=10, alpha=0.75)[0]
        self.endline = ax.plot([],[], color=linecolor, zorder=10, alpha=0.75)[0]
        self.points = []
        self.fig_open = True
        
    def __call__(self):
        plt.show()
        self.connect()
        while self.fig_open:
            plt.pause(0.01)
            
        # self.fig.canvas.stop_event_loop()
        return self.points
     
    
    
    def onclick(self, event):
        # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
        #       ('double' if event.dblclick else 'single', event.button,
        #         event.x, event.y, event.xdata, event.ydata))
        
        if not event.dblclick and event.button==1:
            self.points.append([event.xdata, event.ydata])
            arrpoints = np.asarray(self.points)
            self.endline.set_xdata([])
            self.endline.set_ydata([])
            self.line.set_xdata(arrpoints[:,0])
            self.line.set_ydata(arrpoints[:,1])
            self.dots.set_xdata(arrpoints[:,0])
            self.dots.set_ydata(arrpoints[:,1])
            self.fig.canvas.draw()
            
        elif event.dblclick:
            self.points.append(self.points[0])
            arrpoints = np.asarray(self.points)
            self.endline.set_xdata([])
            self.endline.set_ydata([])
            self.line.set_xdata(arrpoints[:,0])
            self.line.set_ydata(arrpoints[:,1])
            self.fig.canvas.draw()
            self.fig.canvas.mpl_disconnect(self.cid1)
            self.fig.canvas.mpl_disconnect(self.cid3)
            self.fig.canvas.mpl_disconnect(self.cid4)
            self.fig.canvas.mpl_disconnect(self.cid5)
            self.fig_open = False
            
        elif event.button==3:
            arrpoints = np.asarray(self.points)
            dist = arrpoints - np.array([[event.xdata, event.ydata]])
            dist = np.linalg.norm(dist, axis=1)
            if dist.min() <= 20:
                arrpoints = arrpoints[dist > dist.min()]
                self.points.clear()
                self.points.extend(arrpoints.tolist())
                arrpoints = np.asarray(self.points)
                if arrpoints.shape[0] == 0:
                    arrpoints = np.empty((0,2))
                self.line.set_xdata(arrpoints[:,0])
                self.line.set_ydata(arrpoints[:,1])
                self.dots.set_xdata(arrpoints[:,0])
                self.dots.set_ydata(arrpoints[:,1])
                self.fig.canvas.draw()
          
    def mouse_move(self, event):
        if len(self.points)>0:
            self.endline.set_xdata(np.array([self.points[-1][0], event.xdata]))
            self.endline.set_ydata(np.array([self.points[-1][1], event.ydata]))
            self.fig.canvas.draw()
    
    def mouse_enter(self, event):
        self.cid1 = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.cid3 = self.fig.canvas.mpl_connect('motion_notify_event', self.mouse_move)
        
    def mouse_leave(self, event):
        self.fig.canvas.mpl_disconnect(self.cid1)
        self.fig.canvas.mpl_disconnect(self.cid3)
        self.endline.set_xdata([])
        self.endline.set_ydata([])
        self.fig.canvas.draw()
    
    def escaped(self, event):
        if event.key == "escape" or event.key == "enter":
            self.fig_open = False    
      
    def fig_closed(self, event):
        self.fig_open = False
    
    def connect(self):
        self.cid1 = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.cid2 = self.fig.canvas.mpl_connect('close_event', self.fig_closed)
        self.cid3 = self.fig.canvas.mpl_connect('motion_notify_event', self.mouse_move)
        self.cid4 = self.fig.canvas.mpl_connect('axes_enter_event',self. mouse_enter)
        self.cid5 = self.fig.canvas.mpl_connect('axes_leave_event', self.mouse_leave)
        self.cid6 = self.fig.canvas.mpl_connect('key_press_event', self.escaped)
    

class manual_segmentation():
    def __init__(self, img, roi, im):
        self.img = img
        self.roi = roi
        self.fig_open = True
        self.fig, self.ax = plt.subplots(figsize=(8,8))
        self.im = im
        self.ax.set_title(f"Image {self.im}")

    def __call__(self):
        self.fig.subplots_adjust(bottom=0.2)
        
        imgwroi = self.img.copy()
        imgwroi[self.roi>0, :] = imgwroi[self.roi>0, :]*0.5 + np.array([[0, 0, 127]])
        self.ax.imshow(imgwroi)
        
        axadd = self.fig.add_axes([0.7, 0.05, 0.08, 0.075])
        axsub = self.fig.add_axes([0.79, 0.05, 0.12, 0.075])
        badd = Button(axadd, 'Add')
        badd.on_clicked(self.add_to_roi)
        bsub = Button(axsub, 'Substract')
        bsub.on_clicked(self.sub_from_roi)
        self.fig.canvas.mpl_connect('close_event', self.fig_closed)
        self.fig.canvas.mpl_connect('key_press_event', self.escaped)
        plt.show()
        while True:
            if self.fig_open:
                self.fig.canvas.start_event_loop(0.01)
            else:
                self.fig.canvas.stop_event_loop()
                break
        return self.roi
        
        
    def fig_closed(self, event):
        self.fig_open = False
        
    def escaped(self, event): 
        if event.key == "escape" or event.key == "enter":
            plt.close(self.fig)
            self.fig_open = False
        
    def add_to_roi(self, event):     
        self.fig.canvas.stop_event_loop()
        contour = draw_polygon(self.fig, self.ax, [0, 0.5, 1])()
        addroi = cv2.fillPoly(np.zeros_like(self.roi), np.asarray([contour], dtype=np.int32), color=1)    
        self.roi += addroi
        self.roi = (self.roi>0)*1
        imgwroi = self.img.copy()
        imgwroi[self.roi>0, :] = imgwroi[self.roi>0, :]*0.5 + np.array([[0, 0, 127]])
        self.ax.clear()
        self.ax.imshow(imgwroi)
        self.ax.set_title(f"Image {self.im}")

        plt.show()
        plt.pause(0.01)
        # self.fig.canvas.start_event_loop()

    
    def sub_from_roi(self, event):
        self.fig.canvas.stop_event_loop()
        contour = draw_polygon(self.fig, self.ax, [1, 0.15, 0.15])()
        subroi = cv2.fillPoly(np.zeros_like(self.roi), np.asarray([contour], dtype=np.int32), color=1)    
        self.roi -= subroi
        self.roi = (self.roi>0)*1
        imgwroi = self.img.copy()
        imgwroi[self.roi>0, :] = imgwroi[self.roi>0, :]*0.5 + np.array([[0, 0, 127]])
        self.ax.clear()
        self.ax.imshow(imgwroi)
        self.ax.set_title(f"Image {self.im}")

        plt.show()
        plt.pause(0.01)
        # self.fig.canvas.start_event_loop()

#%%


#### Loading, Preprocessing and getting dimensions ####
root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
root.update()

data_path = tk.filedialog.askdirectory(title='Select a folder with heart images')
root.destroy()

# file_path = 'C:/Users/Joao Salinet/Downloads/silhs1.mat'


root = tk.Tk()
root.withdraw()

thereismask = False

if os.path.exists(data_path + "\mask.npy"):
    result = messagebox.askquestion("Segmentation file", "Segmentation file already exists. Load it? (Otherwise, redo segmentation and overwrite the file)")
    if result == "yes":
        mask = np.load(data_path + "\mask.npy")
        thereismask = True


root.destroy()
        
    
#%%
if not thereismask:
    def refine_roi(roi):
        roi = cv2.erode((roi*255).astype(np.uint8),kern,iterations=2).astype(np.uint8)
        roi = cv2.dilate((roi*255).astype(np.uint8),kern,iterations=4).astype(np.uint8)
        roi = cv2.erode((roi*255).astype(np.uint8),kern,iterations=2).astype(np.uint8)
        
        roi = cv2.connectedComponents(roi)[1]
        
        dists = np.zeros(np.unique(roi).shape[0]) + np.inf
        
        for i in range(roi.max() + 1):
            disti = np.asarray(np.where(roi==i)) - np.array([[roi.shape[1]/2], [roi.shape[0]/2]])
            if disti.size/2 < roi.shape[0]*roi.shape[1]*0.25:
                continue

            disti = np.linalg.norm(disti, axis=0).mean()
            dists[i] = disti
            
        roi = roi == dists.argmin()

        return roi
    
    images = []
    order = []
    for root, dirs, files in os.walk(data_path):
        for file in files:
            if file.endswith('png') or file.endswith('tif'):
                order.append(int(file[-7:-4]))
                images.append(cv2.imread(os.path.join(root, file)))
            elif file.endswith('tiff'):
                order.append(int(file[-8:-5]))
                images.append(cv2.imread(os.path.join(root, file)))
    
    
    images = np.asarray(images)
    order = np.argsort(order) # Order is the sequence of indices that puts images in order according to the number in the file name        
    
    images = images[order]
    images = images.astype(np.uint8)
    
    images = images.transpose(0, 2, 1, 3)[:,::-1,:,:]
    
    
    mask = []
    first_iter = True
    kern = cv2.circle(np.zeros((11,11), int), (5,5), 4, 1, -1).astype(np.uint8)
    for im, img in enumerate(images):
        
        image = img.mean(axis=2).copy().astype(np.uint8)
    
        # ret, th = cv2.threshold(image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)   
        
        
        
        if first_iter:
            hist, bins = np.histogram(image, 256, range=(-0.5, 255.5))
            hist = signal.convolve(hist, np.ones(11)/11, mode="same")
            peaks = signal.find_peaks(hist, height=np.percentile(hist, 90), distance=50)[0]

            th = (peaks[0] + peaks[-1])/2
            first_iter = False
    
            
        initialroi = image > th
        initialroi = refine_roi(initialroi)
        
        imgwroi = img.copy()
        imgwroi[initialroi>0, :] = imgwroi[initialroi>0, :]*0.5 + np.array([[0, 0, 127]])
        
        
        fig, ax = plt.subplots(figsize=(16,10))
        ax.set_title(f"Image {im}")
        fig.subplots_adjust(bottom=0.2)
        ims = ax.imshow(imgwroi)
        axslide = fig.add_axes([0.25, 0.05, 0.5, 0.075])
        th_slider = Slider(ax=axslide,
                             label='threshold',
                             valmin=0,
                             valmax=255,
                             valinit=th,
                             )
        


        def update(val):
            global ims, imgwroi, fig, img, initialroi
                    
            initialroi = image > th_slider.val
            initialroi = refine_roi(initialroi)

            imgwroi = img.copy()
            imgwroi[initialroi>0, :] = imgwroi[initialroi>0, :]*0.5 + np.array([[0, 0, 127]])
            
            ims.set_array(imgwroi)
            # plt.show()
            # plt.pause(0.01)
        
        th_slider.on_changed(update)
        
        resetax = fig.add_axes([0.81, 0.05, 0.15, 0.075])
        button = Button(resetax, 'Reset', hovercolor='0.975')
        
        def reset(event):
            th_slider.reset()
        button.on_clicked(reset)
        
        def fig_closed(event):
            global fig_open
            fig_open = False
            
        def escaped(event): 
            global fig_open
        
            if event.key == "escape" or event.key == "enter":
                plt.close(fig)
                fig_open = False
                
        fig.canvas.mpl_connect('close_event', fig_closed)
        fig.canvas.mpl_connect('key_press_event', escaped)
        
        
        cidclosed = fig.canvas.mpl_connect('close_event', fig_closed)
        cidkey = fig.canvas.mpl_connect('key_press_event', escaped)
        
        
        
        plt.show()
        
        fig_open = True
            
        mask.append(initialroi.astype(np.float64))

        th = th_slider.val

   
    mask = np.asarray(mask).astype(np.uint8)
    np.save(data_path + "\mask.npy", mask)

    segmenter = manual_segmentation_edit(images, mask)
    mask = segmenter()
      
    mask = np.asarray(mask).astype(np.uint8)
    np.save(data_path + "\mask.npy", mask)

#%%

dsm = 2*np.sum(mask[0] * mask[-1])/(mask[0].sum() + mask[-1].sum())
mask = mask[:-1]
n = mask.shape[0]

# Angle variation along each image, considering the number of images and a constant delta-angle
d_deg = 360/n
if dsm < 0.9:
    root = tk.Tk()
    root.withdraw()
    
    rotationerror = False
    
    if os.path.exists(data_path + "\mask.npy"):
        result = messagebox.askquestion("Possible rotation error", "The first and the last segmentation doesn't seem similar (DSM = " + format(dsm*100, '.2f') + "% < 90%). Do you want to estimate the angle progression automaticaly?")
        if result == "yes":
            rotationerror = True
    
    
    root.destroy()

    if rotationerror:
        print("Finding angle progression (1/3)")
        mask_cont = mask * (1 - lite_erode3D(np.concatenate((mask[(0,),],
                                                              mask,
                                                              mask[(-1,),]), axis=0), iterations=2)[1:-1])
        mn, mh, mw = np.where(mask_cont)
        uniquemh = np.unique(mh)
        uniquemn = np.unique(mn)
        
        mean_mask_width = np.zeros((uniquemn.size, uniquemh.size))
        print("Finding angle progression (2/3)")
        for mhi in tqdm(range(uniquemh.size), total = uniquemh.size):
            secmh = mw[mh==uniquemh[mhi]]
            secmn = mn[mh==uniquemh[mhi]]
            if np.unique(secmn).size < n:
                mean_mask_width[:, mhi] = -1
            else:
                for mni in range(mean_mask_width.shape[0]):
                    sec = secmh[secmn==mni]
                    mean_mask_width[mni, mhi] = (sec.max() + sec.min())/2
        
        mean_mask_width = mean_mask_width[:,mean_mask_width[0,:] != -1]
        
        
        
        def test_func(x, a, b, c, d):
            return a * np.sin(b * x + c) + d
        
        
        params = np.zeros((mean_mask_width.shape[1], 4))
        b0 = 0.75*np.pi/n
        print("Finding angle progression (3/3)")
        for i, signal in tqdm(enumerate(mean_mask_width.T), total=mean_mask_width.shape[1]):
            
            a0 = signal.std() / np.sqrt(2)
            c0 = (signal.argmax() * (b0*0.25) + signal.argmin() * (b0*0.75))/2
            d0 = signal.mean()
            # plt.figure()
            # plt.plot(test_func(np.arange(n), a0, b0, c0, d0))
            # plt.pause(0.1)
            params[i] = optimize.curve_fit(test_func, np.arange(n),
                                             signal, p0=[a0, b0, c0, d0], maxfev = 10000, ftol = 1e-8 )[0]
            # plt.figure()
            # plt.plot(signal)
            # plt.plot(test_func(np.arange(n), params[i,0], params[i,1], params[i,2], params[i,3]))
            # plt.pause(0.1)
        
        prog = params[:,1]
        prog = prog[(0.5*np.pi/n < prog) * (prog < 1.5*np.pi/n)]
        prog = np.median(prog)
        
        d_deg = 360 / (2*np.pi/prog)

#%%
root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
root.update()

calibration_file = tk.filedialog.askopenfilename(title='Select the calibration file',
                                                 filetypes=(('Numpy arrays', '*.npy'),
                                                            ('All files', '*.*')))
root.destroy()

mask = np.load(data_path + "\mask.npy")
mask = np.pad(mask, ((0,0), (200,200),(200,200)))
n, h, w = mask.shape
mask = mask.transpose(1, 2, 0)
calibration = np.load(calibration_file, allow_pickle=True).item()["calibration"][:2,:3]
calibration_factor = calibration[0,0]
calibration /= calibration_factor
for sl in range(int(np.ceil(n/64))):
    mask[:,:,sl*64:(sl+1)*64] = cv2.warpAffine(mask[:,:,sl*64:(sl+1)*64], 
                                                 calibration, mask.shape[1::-1],
                                                 flags=cv2.INTER_LINEAR)
mask = mask.transpose(2, 0, 1)

#%%

# Matricial way to compute mean position along dim=2, then consider this the 
# center of rotation.
# rotat_center = int(np.sum(mask * np.arange(w).reshape(1,1,-1), dtype=float)/mask.sum())


# # Padding for the estimated rotational center to be in the center of each image.
# if rotat_center < w//2:
#     mask = np.pad(mask, ((0,0),(0, 0), (w - 2*rotat_center, 0)), 'constant')
# else:
#     mask = np.pad(mask, ((0,0),(0, 0), (0, 2*rotat_center-w)), 'constant')

# # update dimensions

# # Crop unecessary borders to save memory
# pad = 4
# mask = np.pad(mask, ((0,0),(pad+1, pad+1), (pad+1, pad+1)), 'constant')
# n, h, w = mask.shape    

# whereX, whereY = np.where(mask.sum(axis=0)>0)
# cropX = np.min([np.min(whereX), np.min(np.abs(whereX - h)) -1] )-1
# cropY = np.min([np.min(whereY), np.min(np.abs(whereY - w)) -1] )-1

# mask = mask[:, cropX-pad+1:-cropX+pad-1, cropY-pad+1:-cropY+pad-1]

pad = 4
mask = np.pad(mask, ((0,0),(pad+1, pad+1), (pad+1, pad+1)), 'constant')

n, h, w = mask.shape   

whereY, whereX = np.where(mask.sum(axis=0)>0)
cropY = min(whereY.min(), h - whereY.max())


mask = mask[:,(cropY - pad):h - (cropY - pad), whereX.min()-pad:whereX.max()+pad,]


mask[:, (0, -1), :] = 0
mask[:, :, (0, -1)] = 0


# update dimensions
n, h, w = mask.shape   



#%%
#### Getting volume for each mask in its respective angle ####



# The reconstucted volume to be "carved"
volume = np.ones((w, w, h), np.float64)

# fig = plt.figure();
rot_mat = cv2.getRotationMatrix2D((w/2,w/2), -d_deg, 1.0)

# For each silhouette, the volume is carved and then rotated by d_deg angles, 
# then repeat.
print("Carving geometry", end="")
for i, img in tqdm(enumerate(mask[:]), total=mask.shape[0]):
    if i > 0:
        for sl in range(int(np.ceil(h/64))):
            volume[:,:,sl*64:(sl+1)*64] = cv2.warpAffine(volume[:,:,sl*64:(sl+1)*64], 
                                                         rot_mat, volume.shape[1::-1],
                                                         flags=cv2.INTER_LINEAR)
    
    img = (img > 0 ) *1
    
    # fig.clear()
    # im = plt.imshow(volume[w//2].T.astype(np.float64), cmap='gray'); 

        
    # plt.pause(0.01)    
        
    volume *= img.T[np.newaxis,...]*0.9 + 0.2
    
volume = (volume > 0.2).astype(bool) 


del mask
np.save(os.path.join(data_path, 'volume_pre.npy'), volume)

#%%
print("Creating mesh from volume.", end="")
# Carregando o volume de dados de um arquivo e realizando operações iniciais de pré-processamento
#### "Polishing" the volume, then save ####
volume = np.load(os.path.join(data_path, 'volume_pre.npy')).astype(np.uint8)
volume = np.pad(volume, 1, 'constant')

w += 2
h += 2

limx, limy, limz = np.where(volume)
limx = np.min([limx.min(), w-limx.max()])
limy = np.min([limy.min(), w-limy.max()])
limz = np.min([limz.min(), h-limz.max()])

# Reduzindo as dimensões do volume, removendo partes vazias
print(".", end="")
volume = volume[limx:-limx, limy:-limy, limz:-limz]

# Preenchendo o volume novamente com um valor específico
volume = np.pad(volume, pad//2, 'constant')

# Realizando operações de erosão e dilatação no volume
volume = lite_erode3D(volume, 2); print(".", end="")
volume = lite_dilate3D(volume, 4); print(".", end="")
volume = lite_erode3D(volume, 2); print(".", end="")

# Convertendo o volume em uma representação binária, mantendo apenas a maior componente conectada
volume = (biggest_component3D(volume)>0); print(".", end="")

# Salvando o volume processado em um arquivo (esta linha está comentada)
# np.save(os.path.join(data_path, 'volume2.npy'), volume)


#### Getting triangular mesh from volume ####
# Carregando o volume de dados de um arquivo (esta linha está comentada)
# volume = np.load(os.path.join(data_path, 'volume2.npy'))*1.
w, l, h = volume.shape

# Reduzindo a resolução do volume (downsampling)
volume = volume[::2,::2,::2]

# Aplicando um desfoque ao volume (blurring) e obtendo índices de pontos de interesse
padvol = lite_blur(volume, k=3); print(".", end="")
contind = np.asarray(np.where((volume>0) * (padvol<=24))).T
del padvol

# Aplicando outro desfoque ao volume (blurring)
smvol = lite_blur(volume, k=7); print(".", end="")
del volume

# Definindo filtros de convolução para cálculo de gradientes
sobelx = np.array([[[-1,-2,-1],
                    [-2,-4,-2],
                    [-1,-2,-1]],
                   [[ 0, 0, 0],
                    [ 0, 0, 0],
                    [ 0, 0, 0]],
                   [[ 1, 2, 1],
                    [ 2, 4, 2],
                    [ 1, 2, 1]]])

sobely = sobelx.transpose(1,0,2)
sobelz = sobelx.transpose(1,2,0)

# Calculando gradientes nas direções x, y e z
Ix = sig.convolve(smvol, sobelx, mode='same'); print(".", end="")
Iy = sig.convolve(smvol, sobely, mode='same'); print(".", end="")      
Iz = sig.convolve(smvol, sobelz, mode='same'); print(".", end="")      
del smvol

# Calculando normais a partir dos gradientes
grad = np.c_[Ix[...,np.newaxis],
             Iy[...,np.newaxis],
             Iz[...,np.newaxis]]

normals = grad[contind[:,0],contind[:,1],contind[:,2],:].astype(np.float64).copy()
del grad, Ix, Iy, Iz

# Normalizando as normais e ajustando seus valores
normals /= np.linalg.norm(normals, axis=1, keepdims=True) + 1e-16
normals *= np.mean(np.std(contind, axis=0))*0.01

# Criando um objeto PointCloud com coordenadas e normais
print(".", end="")
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(contind)
pcd.normals = o3d.utility.Vector3dVector(normals)

# Criando uma malha triangular usando o algoritmo Poisson
print(".", end="")
mesh, den = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8)
mesh.compute_vertex_normals()

# Mesclando vértices próximos na malha
meanptp = np.asarray(mesh.vertices).ptp(axis=0).mean()
for i in range(10):
    min_dist = (10**((i-5)/5))*meanptp/1000
    mesh.merge_close_vertices(min_dist)
    mesh.remove_degenerate_triangles()
    mesh.remove_duplicated_triangles()
    mesh.remove_duplicated_vertices()
    mesh.remove_unreferenced_vertices()
    
if np.asarray(mesh.triangles).shape[0] < 80000:
    mesh = mesh.subdivide_loop(int(np.ceil(20000/np.asarray(mesh.triangles).shape[0])))
    
mesh = mesh.simplify_quadric_decimation(40000)

# Aplicando filtragem de suavização na malha
print(".", end="")
mesh = mesh.filter_smooth_simple(number_of_iterations=5)
mesh.compute_vertex_normals()
mesh.remove_degenerate_triangles()
mesh.remove_duplicated_triangles()
mesh.remove_unreferenced_vertices()
mesh.compute_vertex_normals()


vertices = np.asarray(mesh.vertices)
vertices -= np.array([[w//2, l//2, h//2]])
vertices *= calibration_factor
vertices = vertices*2
triangles = np.asarray(mesh.triangles)

print(" Complete!")

print("Saving...")

outfile = {}
outfile['__header__'] = str(b'MATLAB 5.0 MAT-file Platform: Python 3, Created on: ', encoding='utf-8')  + datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
outfile['__version__'] = '1.0'
outfile['__globals__'] = []
outfile['vertices'] = vertices
outfile['faces'] = triangles + 1

sio.savemat(data_path + '/geometry.mat', outfile)
print("Mesh saved in:", data_path + '/geometry.mat')

# root.destroy()

# Create a visualization window
import vedo
scene = vedo.Plotter().parallel_projection(True)

scene += vedo.Mesh([vertices * np.array([[1, 1, -1]]), triangles]).cmap(input_cmap="gray", input_array=vertices[:,2],).lighting("plastic")
# scene += vedo.Mesh([vertices * np.array([[1, 1, -1]]), triangles]).cmap(input_cmap="jet", input_array=vertices[:,2],).lighting("plastic").add_scalarbar()
scene.show()
scene.interactive()
scene.close()

print("Script ended successfully")

# %%