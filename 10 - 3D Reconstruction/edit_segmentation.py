# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 19:46:38 2023

@author: HEartLab
"""

import numpy as np
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import cv2



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
    
#%%
class manual_segmentation_edit():
    def __init__(self, imgs, rois):
        self.i = 0
        self.imgs = imgs
        self.rois = rois
        
        self.img = imgs[self.i]
        self.roi = rois[self.i]
        self.fig_open = True
        self.fig, self.ax = plt.subplots(figsize=(10,6))
        self.fig.suptitle("Manual Segmentation")

    def __call__(self):
        self.fig.subplots_adjust(bottom=0.2)

        # self.updatefig()
        
        axadd = self.fig.add_axes([0.70, 0.15, 0.075, 0.075])
        axsub = self.fig.add_axes([0.78, 0.15, 0.075, 0.075])
        axprev = self.fig.add_axes([0.05, 0.05, 0.15, 0.075])
        axconc = self.fig.add_axes([0.425, 0.05, 0.15, 0.075])
        axnext = self.fig.add_axes([0.80, 0.05, 0.15, 0.075])
        
        self.badd = Button(axadd, '+', color=[0, 0.5, 1], hovercolor=[0, 0.7, 1])
        self.badd.on_clicked(self.add_to_roi)
        self.bsub = Button(axsub, '-', color=[1, 0.15, 0.15], hovercolor=[1, 0.21, 0.21])
        self.bsub.on_clicked(self.sub_from_roi)
        self.bprev = Button(axprev, "<< Previous")
        self.bprev.on_clicked(self.previmg)
        self.bconc = Button(axconc, "Conclude")
        self.bconc.on_clicked(self.conclude)        
        self.bnext = Button(axnext, "Next >>")
        self.bnext.on_clicked(self.nextimg)
        
        self.bprev.color = '0.75'
        self.bprev.hovercolor = '0.75'
        # plt.pause(0.1)
        self.bprev.set_active(False)
        self.bprev.enabled = False
        self.bprev.color = '0.75'
        self.bprev.hovercolor = '0.75'
        self.bprev.disabled_color = '0.75'
        self.bnext.disabled_color = '0.75'

        if self.rois.shape[0] == 1:
            self.bnext.color = '0.75'
            self.bnext.hovercolor = '0.75'
            # plt.pause(0.1)
            self.bnext.set_active(False)
            self.bnext.enabled = False

        self.fig.canvas.mpl_connect('close_event', self.fig_closed)
        self.fig.canvas.mpl_connect('key_press_event', self.escaped)
        self.updatefig()
        
        # while True:
        #     if self.fig_open:
        #         self.fig.canvas.start_event_loop(0.01)
        #     else:
        #         self.fig.canvas.stop_event_loop()
        #         break
            
        # plt.close(self.fig)
        return self.rois
        
    def updatefig(self):
        imgwroi = self.img.copy()
        imgwroi[self.roi>0, :] = imgwroi[self.roi>0, :]*0.5 + np.array([[0, 0, 127]])
        self.ax.clear()
        self.ax.imshow(imgwroi)
        self.ax.set_title(f"Image {self.i}")
        plt.show()
        plt.draw()
        plt.pause(0.01)
    
    # def disableall(self):
        
    
    # def enableall(self):
        
        
    
    def nextimg(self, event):
        self.fig.canvas.stop_event_loop()
        self.rois[self.i] = self.roi
        self.i += 1
        self.roi = self.rois[self.i]
        self.img = self.imgs[self.i]
        self.bprev.set_active(True)
        self.bprev.color = '0.85'
        self.bprev.hovercolor = '0.95'

        if self.i == self.rois.shape[0] -1:
            self.bnext.color = '0.75'
            self.bnext.hovercolor = '0.75'
            plt.pause(0.1)
            self.bnext.set_active(False)
        self.updatefig()

    def previmg(self, event):
        self.fig.canvas.stop_event_loop()
        self.rois[self.i] = self.roi
        self.i -= 1
        self.roi = self.rois[self.i]
        self.img = self.imgs[self.i]
        self.bnext.set_active(True)
        self.bnext.color = '0.85'
        self.bnext.hovercolor = '0.95'

        if self.i == 0:
            self.bprev.color = '0.75'
            self.bprev.hovercolor = '0.75'
            plt.pause(0.1)
            self.bprev.set_active(False)
        self.updatefig()

    def conclude(self, event):
        self.fig.canvas.stop_event_loop()
        self.rois[self.i] = self.roi
        plt.close(self.fig)
        self.fig_open = False                
        
    
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
        self.roi[addroi>0] = 1 
                
        self.updatefig()
        # self.fig.canvas.start_event_loop()

    
    def sub_from_roi(self, event):
        self.fig.canvas.stop_event_loop()
        contour = draw_polygon(self.fig, self.ax, [1, 0.15, 0.15])()
        subroi = cv2.fillPoly(np.zeros_like(self.roi), np.asarray([contour], dtype=np.int32), color=1)    
        self.roi[subroi > 0] = 0

        self.updatefig()
        # self.fig.canvas.start_event_loop()


if __name__ == "__main__":
    import os
    from tkinter import filedialog
    import tkinter as tk
    from matplotlib import use
        
    use('tkagg')


    #%%
    #### Loading, Preprocessing and getting dimensions ####
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    root.update()

    data_path = filedialog.askdirectory(title='Select a folder with heart images')

    root.withdraw()
    root.destroy()


    #%%

    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    root.update()

    mask_path = filedialog.askopenfilename(title='Select the segmentation file', filetypes=(("Numpy files", "*.npy"),))

    root.withdraw()
    root.destroy()

    mask = np.load(mask_path).astype(np.uint8)

    #%%
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
    #%%

    segmenter = manual_segmentation_edit(images, mask)
    mask = segmenter()
    mask = np.asarray(mask).astype(np.uint8)


    #%%
    def button_clicked(option,):
        global root
        root.option = option
        root.withdraw()
        root.destroy()

    # Usage
    while True:
        root = tk.Tk()
        root.attributes("-topmost", True)
        root.update()

        root.title("Save segmentation")

        # Message label
        message_label = tk.Label(root, text="Overwrite original file?")
        message_label.pack(pady=20)

        # Option buttons
        button1 = tk.Button(root, text="Yes", command=lambda: button_clicked(1))
        button1.pack(side=tk.LEFT, padx=10)

        button2 = tk.Button(root, text="Save as", command=lambda: button_clicked(2))
        button2.pack(side=tk.LEFT, padx=10)

        button3 = tk.Button(root, text="Discard alterations", command=lambda: button_clicked(3))
        button3.pack(side=tk.LEFT, padx=10)
        # Set the option attribute to None initially
        root.option = None

        root.mainloop()
        
        if root.option == 1:
            np.save(mask_path, mask)
            print("Segmentation saved in: ", mask_path)
            print("Program terminated")
            break
        
        elif root.option == 2:
            
            root = tk.Tk()
            root.withdraw()
            root.attributes("-topmost", True)
            root.update()
            
            saveas_mask_path = filedialog.asksaveasfilename(initialdir = os.path.dirname(mask_path),
                                                            initialfile = "mask.npy",
                                                            title = "Save the segmentation file",
                                                            filetypes = (("Numpy array", "*.npy"),))
            root.withdraw()
            root.destroy()
            
            
            if saveas_mask_path == '':
                print("\"Save as\" canceled")
                continue
            
            else:
                np.save(saveas_mask_path, mask)
                print("Segmentation saved in: ", saveas_mask_path)
                print("Program terminated")

                break
            
        elif root.option == 3:
            root = tk.Tk()
            root.attributes("-topmost", True)
            root.update()
            
            root.title("Discard alterations")

            # Message label
            message_label = tk.Label(root, text="Are you sure you want to discard the changes?")
            message_label.pack(pady=20)

            # Option buttons
            button1 = tk.Button(root, text="Yes", command=lambda: button_clicked(1))
            button1.pack(side=tk.LEFT, padx=10)

            button2 = tk.Button(root, text="No", command=lambda: button_clicked(2))
            button2.pack(side=tk.LEFT, padx=10)
            
            root.mainloop()
            
            if root.option == 1:
                print("Changes discarded, program terminated")
                break
            
            elif root.option == 2:
                continue
            
            else:
                continue
        else:
            continue
            