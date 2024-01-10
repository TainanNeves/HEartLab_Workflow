#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Tainan Cerqueira Neves | Federal University of ABC - HEartLab
@data: 2023-08-01

Internal library for analyzing electrical and optical data from the perfused rabbit heart.
HEartLab
"""





#Librarys
import numpy as np
from ast import literal_eval
from glob import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
from scipy.io import savemat
from msilib.schema import TextStyle
from tkinter import font
import plotly.express as px
from datetime import datetime





"""
Functions from the library Binary.py, from the openEphys developers
"""
def ApplyChannelMap(Data, ChannelMap):
    print('Retrieving channels according to ChannelMap... ', end='')
    for R, Rec in Data.items():
        if Rec.shape[1] < len(ChannelMap) or max(ChannelMap) > Rec.shape[1]-1:
            print('')
            print('Not enough channels in data to apply channel map. Skipping...')
            continue

        Data[R] = Data[R][:,ChannelMap]

    return(Data)





def BitsToVolts(Data, ChInfo, Unit):
    print('Converting to uV... ', end='')
    Data = {R: Rec.astype('float32') for R, Rec in Data.items()}

    if Unit.lower() == 'uv': U = 1
    elif Unit.lower() == 'mv': U = 10**-3

    for R in Data.keys():
        for C in range(len(ChInfo)):
            Data[R][:,C] = Data[R][:,C] * ChInfo[C]['bit_volts'] * U
            if 'ADC' in ChInfo[C]['channel_name']: Data[R][:,C] *= 10**6

    return(Data)





def Load(Folder, Processor=None, Experiment=None, Recording=None, Unit='uV', ChannelMap=[]):
    Files = sorted(glob(Folder+'/**/*.dat', recursive=True))
    InfoFiles = sorted(glob(Folder+'/*/*/structure.oebin'))


    Data, Rate = {}, {}
    for F,File in enumerate(Files):
        File = File.replace('\\', '/') # Replace windows file delims
        Exp, Rec, _, Proc = File.split('/')[-5:-1]
        Exp = str(int(Exp[10:])-1)
        Rec = str(int(Rec[9:])-1)
        Proc = Proc.split('.')[0].split('-')[-1]
        if '_' in Proc: Proc = Proc.split('_')[0]

        if Proc not in Data.keys(): Data[Proc], Rate[Proc] = {}, {}

        if Experiment:
            if int(Exp) != Experiment-1: continue

        if Recording:
            if int(Rec) != Recording-1: continue

        if Processor:
            if Proc != Processor: continue

        print('Loading recording', int(Rec)+1, '...')
        if Exp not in Data[Proc]: Data[Proc][Exp] = {}
        Data[Proc][Exp][Rec] = np.memmap(File, dtype='int16', mode='c')


        Info = literal_eval(open(InfoFiles[F]).read())
        ProcIndex = [Info['continuous'].index(_) for _ in Info['continuous']
                     if str(_['source_processor_id']) == Proc][0] # Changed to source_processor_id from recorded_processor_id

        ChNo = Info['continuous'][ProcIndex]['num_channels']
        if Data[Proc][Exp][Rec].shape[0]%ChNo:
            print('Rec', Rec, 'is broken')
            del(Data[Proc][Exp][Rec])
            continue

        SamplesPerCh = Data[Proc][Exp][Rec].shape[0]//ChNo
        Data[Proc][Exp][Rec] = Data[Proc][Exp][Rec].reshape((SamplesPerCh, ChNo))
        Rate[Proc][Exp] = Info['continuous'][ProcIndex]['sample_rate']

    for Proc in Data.keys():
        for Exp in Data[Proc].keys():
            if Unit.lower() in ['uv', 'mv']:
                ChInfo = Info['continuous'][ProcIndex]['channels']
                Data[Proc][Exp] = BitsToVolts(Data[Proc][Exp], ChInfo, Unit)

            if ChannelMap: Data[Proc][Exp] = ApplyChannelMap(Data[Proc][Exp], ChannelMap)

    print('Done.')

    return(Data, Rate)






"""
Functions from HEartLab UFABC
"""



"""
Eletric Functions
"""

def loadingData(folder, experiment = 1, recording = 1):
    #loading data
    data, sample_rate = Load(folder, Experiment=experiment, Recording=recording)
    #loading sample rate
    sample_rate = sample_rate[list(sample_rate.keys())[0]][list(data[list(sample_rate.keys())[0]])[0]]
    #loading data to a pd.df
    data = data[list(data.keys())[0]][list(data[list(data.keys())[0]].keys())[0]][list(data[list(data.keys())[0]][list(data[list(data.keys())[0]].keys())[0]].keys())[0]]
    data = pd.DataFrame(data)
    #tamanho do arquivo e tempo
    size = data.shape[0]
    n_electrodes = data.shape[1]
    time = size/sample_rate
    
    print("Sample Rate: {}hz \nElectrodes: {} \nRecording time: {}s".format(sample_rate, n_electrodes, time))
    return data, sample_rate, n_electrodes, time
    


    
    
def dataTo_csv(data, file_name="exported_signal.csv"):
    df = pd.DataFrame(data)
    df.to_csv(file_name, index=False)

    



def dataTo_mat(data, file_name="exported_data.mat", sample_rate=4000, n_electrodes=64, recording_time=1):
    transposed_data = data.T  # Transpose the data
    mat_data = {
        'D': {
            'data': transposed_data.to_numpy(),
            'sampleRate': sample_rate,
            'n_electrodes': n_electrodes,
            'recording_time': recording_time
        }
    }
    savemat(file_name, mat_data, format='5')
    
    
    
    
      
def genericPlot(data, sample_rate, startValue=0, endValue=""):
    n_electrodes = data.shape[1] #n of electrodes
    if endValue == "": #If user do not put a endValue
        endValue = data.shape[0]

    time = np.linspace(startValue/sample_rate, endValue/sample_rate, data[0][startValue:endValue].size) #Time array in seconds
    fig, axs = plt.subplots(n_electrodes, dpi=200, figsize=(18, 4*n_electrodes)) #Subplots for all electrodes
    plt.xlabel("Tempo (s)", fontsize=20)
    plt.ylabel("Sinal (uV)", loc="center", fontsize =20)

    for i in range(n_electrodes):
        axs[i].plot(time, data[i][startValue:endValue]) #Individual PLot
        title = "Electrode: {}".format(i+1)
        axs[i].set_title(title)
        
    plt.tight_layout() #To not overlap titles
    plt.show() #To show plots
    

    
    

def specificPlot(data, sample_rate, startValue=0, endValue="", electrodes=[1]):
    if endValue == "": #if the user do not put a endValue
        endValue = data.shape[0]
        
    for i in electrodes:
        time = np.linspace(startValue/sample_rate, endValue/sample_rate, data[i-1][startValue:endValue].size)
        fig = px.line(data[i-1][startValue:endValue], x=time, y=(i-1), title='Electrode: {}'.format(i), labels={"x": "Time (s)", "y": "Electrode Reading (uV)"})
        fig.show()


        
        

def comparasionPlot(data, sample_rate, startValue=0, endValue="", electrodes=[1]):
    if endValue =="": #if the user do not put a endValue
        endValue = data.shape[0]
        
    for i in electrodes:
        time = np.linspace(startValue/sample_rate, endValue/sample_rate, data[i-1][startValue:endValue].size)
        plt.plot(time, data[i-1][startValue:endValue])
    
    plt.title("Electrode comparation", fontsize=15)
    plt.xlabel("Tempo [s]", fontsize=10)
    plt.ylabel("Signal [uV]", fontsize=10)
    
    leg = []
    for el in electrodes:
        leg.append("Electrode: {}".format(el))
    plt.legend(leg)

    plt.tight_layout() #To not overlap titles
    plt.show()





"""
Optic Functions
"""

def norpix2python(file_name, show_img=False, frames=0):
    # Helper function to read data from the file
    def read_from_file(offset, data_type, size):
        fid.seek(offset)
        return np.fromfile(fid, dtype=data_type, count=size)

    # Open file for reading at the beginning
    try:
        fid = open(file_name, 'rb')
    except FileNotFoundError:
        print('The specified filename does not exist')
        return None, None

    # HEADER INFORMATION
    header_info = {}

    # Read version
    header_info['Version'] = read_from_file(28, 'int32', 1)[0]

    # Read header size
    header_info['HeaderSize'] = read_from_file(32, 'uint32', 1)[0]
    if header_info['Version'] >= 5:
        print('Version 5+ detected, overriding reported header size')
        header_info['HeaderSize'] = 8192

    # Read description (Handle UTF-16 and ASCII encoding)
    DescriptionFormat = read_from_file(592, 'uint32', 1)[0]
    header_info['Description'] = read_from_file(36, 'uint8', 512).tobytes()

    try:
        header_info['Description'] = header_info['Description'].decode('utf-16')
    except UnicodeDecodeError:
        try:
            header_info['Description'] = header_info['Description'].decode('ascii')
        except UnicodeDecodeError:
            header_info['Description'] = header_info['Description'].decode('latin-1')

    # Read image information
    tmp = read_from_file(548, 'uint32', 6)
    header_info['ImageWidth'] = tmp[0]
    header_info['ImageHeight'] = tmp[1]
    header_info['ImageBitDepth'] = tmp[2]
    header_info['ImageBitDepthReal'] = tmp[3]
    header_info['ImageSizeBytes'] = tmp[4]
    vals = [0, 100, 101, 200, 300, 400, 500, 600, 610, 620, 700, 800, 900]
    fmts = ['Unknown', 'Monochrome', 'Raw Bayer', 'BGR', 'Planar', 'RGB',
            'BGRx', 'YUV422', 'YUV422_20', 'YUV422_PPACKED', 'UVY422', 'UVY411', 'UVY444']
    header_info['ImageFormat'] = fmts[vals.index(tmp[5])]

    # Read number of allocated frames
    header_info['AllocatedFrames'] = read_from_file(572, 'uint16', 1)[0]

    # Read origin
    header_info['Origin'] = read_from_file(576, 'uint16', 1)[0]

    # Read true image size
    header_info['TrueImageSize'] = read_from_file(580, 'uint32', 1)[0]

    # Read frame rate
    header_info['FrameRate'] = read_from_file(584, 'double', 1)[0]

    # Reading images
    image_offset = 8192
    image_data = None
    timestamps = []

    if frames == 0:
        num_frames_to_process = header_info['AllocatedFrames']
    else:
        num_frames_to_process = frames

    if header_info['ImageBitDepth'] == 8:
        bit_depth = np.uint8
    elif header_info['ImageBitDepth'] == 16:
        bit_depth = np.uint16
    else:
        print('Unsupported bit depth')
        return header_info, None

    img_out = np.zeros((header_info['ImageHeight'], header_info['ImageWidth'], num_frames_to_process), dtype=bit_depth)

    for i in range(num_frames_to_process):
        print('Reading Image Frame', i + 1)
        fid.seek(image_offset + i * header_info['TrueImageSize'])
        img_data = np.fromfile(fid, dtype=bit_depth, count=header_info['ImageSizeBytes'] // np.dtype(bit_depth).itemsize)
        img_out[:, :, i] = np.reshape(img_data, (header_info['ImageHeight'], header_info['ImageWidth']))

        # Read timestamp
        time_secs_posix = np.fromfile(fid, dtype=np.int32, count=1)[0]
        sub_seconds = np.fromfile(fid, dtype=np.uint16, count=2)
        time_date_num = time_secs_posix / 86400 + datetime(1970, 1, 1).toordinal()
        timestamp = f"{datetime.fromordinal(int(time_date_num)).strftime('%Y-%m-%d')}:{sub_seconds[0]:04d}{sub_seconds[1]:04d}"
        timestamps.append(timestamp)

    if show_img:
        import matplotlib.pyplot as plt

        def display_image(img):
            plt.imshow(img, cmap='gray')
            plt.show()

        for i in range(num_frames_to_process):
            display_image(img_out[:, :, i])

    # Close the file
    fid.close()

    return header_info, img_out





def binning(data, bin_size, precision='uint16'):
    return np.round(np.mean(np.mean(data.reshape(data.shape[0] // bin_size, bin_size, 
                                                data.shape[1] // bin_size, bin_size, 
                                                data.shape[2]), (1, 3)), (1, 2))).astype(precision)





def convert_to_python(binSize):
    bin_size = binSize
    precision = 'uint16'  # for binning

    String = ['Camera 1', 'Camera 2', 'Camera 3']

    for N in range(3):
        os.chdir('Optical Mapping')
        os.chdir(String[N])
        DIR = os.listdir('.seq')
        
        for filename in DIR:
            P = os.getcwd()
            header, DATA = norpix2python(filename, 0, 0)
            os.chdir('..//..')
            os.chdir('Optical Mapping to Python')
            
            plt.figure()
            plt.plot(np.mean(np.mean(DATA[500:504, 500:504, :], axis=0), axis=0))
            plt.savefig(f'{filename}_Trace_Cam{N + 1}.png')
            plt.close()

            plt.figure()
            plt.imshow(DATA[:, :, 100], cmap='gray')
            plt.axis('square')
            plt.savefig(f'{filename}_Image_Cam{N + 1}.png')
            plt.close()

            DATA_old = DATA
            DATA = binning(DATA, bin_size, precision)
            #savemat(f'{filename[:-4]}_Bin={bin_size}_Cam{N + 1}.mat', {'DATA': DATA})
            np.save(f'{filename[:-4]}_Bin={bin_size}_Cam{N + 1}.npy', DATA) #To load use np.load("archive_patch.npy")
            os.chdir(P)
            
        os.chdir('..//..')


        
        
        
def pick_up_a_trace(I, DATA, N):
    button = 0
    
    while button != 32: # loop continues until the spacebar is pressed
        plt.figure(99)
        plt.imshow(I, cmap='gray')
        plt.axis('square')
        y, x, button = plt.ginput(1)[0]
        x, y = int(round(x)), int(round(y))

        if button == 3: # right-clicks, the intensity value of the pixel at (x, y) in the image I is set to the maximum intensity value in the entire image.
            I[x, y] = I.max()
            plt.imshow(I)

        if x <= 0 or y <= 0:
            break

        X = np.mean(np.mean(DATA[x - N:x + N, y - N:y + N, :], axis=0), axis=0)
        Trace1 = X

        plt.figure(103)
        plt.plot(X)
        plt.title(f'row: {x} col: {y}')

        display = (X.max() - X.min()) * 100
        print(display)















