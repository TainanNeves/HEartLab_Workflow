
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:58:22 2023

@author: HEartLab
"""
#%%
import tkinter as tk
veryroot = tk.Tk()
tk.Label(veryroot, text="Calibration is running").pack()

from tkinter import filedialog
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
import matplotlib
import struct
# import h5py
from PIL import Image
matplotlib.use('tkagg')

#%%

def Open_StreamPix(fileName, showImg=False):
    # Open n-bit Norpix image sequence in MATLAB
    
    # Open file for reading at the beginning
    fid = open(fileName, 'rb')
    if fid.closed:
        print('The specified filename does not exist')
        return
    
    endianType = '<'
    
    # Read header
    fid.seek(28, 0)
    headerInfo = {}
    headerInfo['Version'] = struct.unpack(endianType + 'i', fid.read(4))[0]
    
    fid.seek(32, 0)
    headerInfo['HeaderSize'] = struct.unpack(endianType + 'i', fid.read(4))[0]
    if headerInfo['Version'] >= 5:
        print('Version 5+ detected, overriding reported header size')
        headerInfo['HeaderSize'] = 8192
    
    fid.seek(592, 0)
    DescriptionFormat = struct.unpack(endianType + 'l', fid.read(4))[0]
    
    fid.seek(36, 0)
    headerInfo['Description'] = fid.read(46).decode('utf-16') if DescriptionFormat == 0 else fid.read(46).decode('utf-8')
    
    fid.seek(548, 0)
    tmp = struct.unpack(endianType + '24I', fid.read(24*4))
    headerInfo['ImageWidth'] = tmp[0]
    headerInfo['ImageHeight'] = tmp[1]
    headerInfo['ImageBitDepth'] = tmp[2]
    headerInfo['ImageBitDepthReal'] = tmp[3]
    headerInfo['ImageSizeBytes'] = tmp[4]
    vals = [0, 100, 101, *range(200, 601, 100), 610, 620, 700, 800, 900]
    fmts = ['Unknown', 'Monochrome', 'Raw Bayer', 'BGR', 'Planar', 'RGB', 'BGRx', 'YUV422', 'YUV422_20', 'YUV422_PPACKED', 'UVY422', 'UVY411', 'UVY444']
    headerInfo['ImageFormat'] = fmts[vals.index(tmp[5])]
    
    fid.seek(572, 0)
    headerInfo['AllocatedFrames'] = struct.unpack(endianType + 'H', fid.read(2))[0]
    
    fid.seek(576, 0)
    headerInfo['Origin'] = struct.unpack(endianType + 'H', fid.read(2))[0]
    
    fid.seek(580, 0)
    headerInfo['TrueImageSize'] = struct.unpack(endianType + 'I', fid.read(4))[0]
    
    fid.seek(584, 0)
    headerInfo['FrameRate'] = struct.unpack(endianType + 'd', fid.read(8))[0]
    
    # Reading images
    # Following the header, each image is stored and aligned on a 8192 bytes boundary (when no metadata is included)
    imageOffset = 8192
    
    if showImg:
        import matplotlib.pyplot as plt
    
    imgOut = None
    bitstr = ''
    numFramesToProcess = headerInfo['AllocatedFrames']
    
    if headerInfo['ImageBitDepth'] == 8:
        BitDepth = np.uint8
    elif headerInfo['ImageBitDepth'] == 16:
        BitDepth = np.uint16
    
    if headerInfo['ImageBitDepthReal'] == 8:
        bitstr = 'B'
    elif headerInfo['ImageBitDepthReal'] in [10, 12, 14, 16]:
        bitstr = 'H'
    
    if bitstr == '':
        raise ValueError('Unsupported bit depth')
    
    nread = 0
    for i in range(numFramesToProcess):
        print(f'Reading Image Frame {i+1}')
        fid.seek(imageOffset + nread * headerInfo['TrueImageSize'], 0)
    
        if headerInfo['ImageFormat'] == 'Monochrome':
            numPixels = headerInfo['ImageWidth'] * headerInfo['ImageHeight']
            data = struct.unpack(endianType + bitstr * numPixels, fid.read(numPixels * struct.calcsize(bitstr)))
            image = np.array(data, dtype=BitDepth).reshape(headerInfo['ImageHeight'], headerInfo['ImageWidth'])
            if showImg:
                plt.imshow(image, cmap='gray')
                plt.title(f'frame {nread+1}')
                plt.show()
            if imgOut is None:
                imgOut = np.zeros((headerInfo['ImageHeight'], headerInfo['ImageWidth'], numFramesToProcess), dtype=BitDepth)
            imgOut[:, :, i] = image
    
        elif headerInfo['ImageFormat'] == 'UVY422':
            bitstr = 'B'
            numPixels = headerInfo['ImageWidth'] * headerInfo['ImageHeight'] * 2
            data = fid.read(numPixels)
            
            if not data:
                break
            
            UYVYimg = np.array(struct.unpack(endianType + bitstr * numPixels, data), dtype=np.uint8).reshape(headerInfo['ImageHeight'], headerInfo['ImageWidth'] * 2)
            
            YUV = np.zeros((headerInfo['ImageHeight'], headerInfo['ImageWidth'], 3), dtype=np.uint8)
            YUV[:, :, 0] = UYVYimg[:, 1::2]
            YUV[:, 0::2, 1] = UYVYimg[:, 0::4]
            YUV[:, 1::2, 1] = UYVYimg[:, 0::4]
            YUV[:, 0::2, 2] = UYVYimg[:, 2::4]
            YUV[:, 1::2, 2] = UYVYimg[:, 2::4]
            imgOut = Image.fromarray(YUV, 'YCbCr').convert('RGB')
    
        elif headerInfo['ImageFormat'] == 'YUV422_PPACKED':
            numPixels = headerInfo['ImageWidth'] * headerInfo['ImageHeight'] * 2
            data = fid.read(numPixels)
            imgOut = np.array(struct.unpack(endianType + 'H' * numPixels, data), dtype=np.uint16).reshape(headerInfo['ImageHeight'], headerInfo['ImageWidth'] * 2)
    
        timeSecsPOSIX = struct.unpack(endianType + 'i', fid.read(4))[0]
        subSeconds = struct.unpack(endianType + 'HH', fid.read(4))
        timeDateNum = timeSecsPOSIX / 86400 + 719529  # 719529 is the datenum value for January 1, 1970
        headerInfo['timestamp'] = headerInfo.get('timestamp', []) + [f'{timeDateNum}:{subSeconds[0]}{subSeconds[1]}']
    
        nread += 1
    
    fid.close()
    
    return headerInfo, imgOut





#%%


root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
root.update()
file = tk.filedialog.askopenfilename(title='Select the calibration prism image file',
                                                 filetypes=(("StreamPix Sequence files", '*.seq'),
                                                            ('MATLAB files', '*.mat'),
                                                            ('Numpy Array', '*.npy'),
                                                            ('All files', '*.*')))
root.destroy()

if file.endswith(".mat"): 
    matfile = sio.loadmat(file)
    try:
        del matfile["__header__"], matfile["__version__"], matfile["__"]
    except NameError:
        []
    img = matfile["camera_A"].T[::-1,:]
elif file.endswith(".seq"):
    _, img = Open_StreamPix(file)
    img = img[:,:,0].T[::-1,:]
    
elif file.endswith(".jpg"):
    img = Image.open(file)
    img = np.array(img)
    img = img[:,:,0].T[::-1,:]
#elif file.endswith(".tiff"):
#    img = Image.open(file)
#    img_array = np.array(img)
#    img=img_array
#    img = img[:,:,0].T[::-1,:]
    
img[img > 5*img.mean()] = 5*img.mean()


fig, ax = plt.subplots()
ax.imshow(img, cmap='gray')
sct = ax.scatter([],[], s = 200, c='r',marker='x',zorder=10, alpha=0.75)

points = []
def onclick(event):
    # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
    #       ('double' if event.dblclick else 'single', event.button,
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
            
fig_open = True      
def fig_closed(event):
    global fig_open
    fig_open = False

cid1 = fig.canvas.mpl_connect('button_press_event', onclick)
cid2 = fig.canvas.mpl_connect('close_event', fig_closed)

plt.show()


while fig_open:
    plt.pause(0.2)


arrpoints = np.asarray(points)
arrpoints = arrpoints[arrpoints[:,0]!=None]
        
#%%

# estimate the central point based on every point distance from the average coordinate
dist = arrpoints - arrpoints.mean(axis=0, keepdims=True)
dist = np.linalg.norm(dist, axis=1)
ind = dist.argmin()
 

fig, ax = plt.subplots()
ax.imshow(img, cmap='gray')
ax.scatter(arrpoints[:,0],arrpoints[:,1], s = 200, c='r',marker='x',zorder=10, alpha=0.75,  picker=True)
sct = ax.scatter(arrpoints[(ind,),0],arrpoints[(ind,),1], s = 200, c='g',marker='x',zorder=10, alpha=1)



def pickpoint(event):
    if isinstance(event.artist, matplotlib.collections.PathCollection):
        global ind
        global ax
        global fig
        # for att in dir(thisline): print(att)
        ind = event.ind[0]      
        sct.set_offsets(arrpoints[ind,:])
        fig.canvas.draw()


fig_open = True
def fig_closed(event):
    global fig_open
    fig_open = False

cid3 = fig.canvas.mpl_connect('close_event', fig_closed)
cid4 = fig.canvas.mpl_connect('pick_event', pickpoint)

while fig_open:
    plt.pause(0.2)

#%%

pointsord = np.ones((15,3)) * np.nan
pointsord[:,2] = 1
realcoord = np.array([[-1, 2, 1],
                      [ 0, 2, 1],
                      [ 1, 2, 1],
                      [-1, 1, 1],
                      [ 0, 1, 1],
                      [ 1, 1, 1],
                      [-1, 0, 1],
                      [ 0, 0, 1],
                      [ 1, 0, 1],
                      [-1,-1, 1],
                      [ 0,-1, 1],
                      [ 1,-1, 1],
                      [-1,-2, 1],
                      [ 0,-2, 1],
                      [ 1,-2, 1]])


# central point
pointsord[7,:2] = arrpoints[ind,:]

dist = arrpoints - arrpoints[ind,:]
dist = np.linalg.norm(dist, axis=1)

centralpoints = arrpoints[dist.argsort()[:9]] # find the 9 central points
edgepoints = arrpoints[dist.argsort()[9:]] # the remaining points

centralpoints = centralpoints[centralpoints[:,1].argsort()].reshape(3,3,2)
sortind = centralpoints[:,:,0].argsort(axis=1)
centralpoints[0,:,:] = centralpoints[0, sortind[0,:,],:]
centralpoints[1,:,:] = centralpoints[1, sortind[1,:,],:]
centralpoints[2,:,:] = centralpoints[2, sortind[2,:,],:]

deltax = np.mean(centralpoints[:,1:,0] - centralpoints[:,:-1,0])
deltay = np.mean(centralpoints[1:,:,1] - centralpoints[:-1,:,1])


pointsord[3:12,:2] = centralpoints.reshape(9,2)

#%%

# T*realcoord = pointsord
# T*realcoord*(realcoord^-1) = pointsord*(realcoord^-1)
# T = pointsord*(realcoord^-1)


T = np.matmul(pointsord[3:12].T, np.linalg.pinv(realcoord[3:12]).T)
estimated_points = np.matmul(T, realcoord.T).T

for point in edgepoints:
    err = estimated_points[:,:2] - point
    err = np.linalg.norm(err, axis=1)
    pointsord[err.argmin(),:2] = point

  
#%%

realcoord = realcoord[~np.isnan(pointsord[:,0])]
pointsord = pointsord[~np.isnan(pointsord[:,0])]
pointsord[:,:2] -= np.array([img.shape[1], img.shape[0]])/2


T = np.matmul(pointsord.T, np.linalg.pinv(realcoord).T)
estimated = np.matmul(T, realcoord.T).T

T[:,0] = np.array([np.linalg.norm(T[:,0]), 0, 0])
T[:,1] = np.array([0, np.linalg.norm(T[:,1]), 0])
T[(0, 1),(0,1)] = T[(0, 1),(0,1)].max()
fixed = np.matmul(T, realcoord.T).T

np.save(file[:-4] + "_calibration.npy", {"calibration":np.linalg.pinv(T),
                                         "base_shape":img.shape})


# fig, ax = plt.subplots()
# ax.imshow(img, cmap='gray')
# ax.scatter(pointsord[:,0],pointsord[:,1], s = 200, c='r',marker='x',zorder=10, alpha=0.75,  picker=True)
# ax.scatter(estimated[:,0],estimated[:,1], s = 200, c='g',marker='x',zorder=10, alpha=0.75,  picker=True)
# ax.scatter(fixed[:,0],fixed[:,1], s = 200, c='b',marker='x',zorder=10, alpha=0.75,  picker=True)


#%%

fig, ax = plt.subplots()
ax.imshow(img, cmap='gray')
# ax.scatter(arrpoints[:,0],arrpoints[:,1], s = 200, c='r',marker='x',zorder=10, alpha=0.75,  picker=True)
sct = ax.scatter(pointsord[(0,),0]+img.shape[1]/2,pointsord[(0,),1]+img.shape[0]/2, s = 200, c='r',marker='x',zorder=10, alpha=1)
fig.show()

for i in range(1, 17):
    plt.pause(0.2)
    sct.set_offsets(pointsord[:i,:2]+np.array([img.shape[1], img.shape[0]])/2)
    fig.canvas.draw()
    
plt.close()

print('Calibration saved!')