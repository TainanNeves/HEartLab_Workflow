# -*- coding: utf-8 -*-
"""
Created on Tue May  3 11:00:02 2022

@author: Joao Salinet
"""

import tkinter as tk
root = tk.Tk()
tk.Label(root, text="Optical Data Projection is running").pack()

import tkinter.filedialog as fd
from tkinter import ttk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import tkinter.filedialog
import scipy.io as sio
import time
from PIL import Image
import os
import h5py
import open3d as o3d
from open_strampix import Open_StreamPix
from mesh_processing import get_connectivity, find_HR_idx, erodemorph3d, dilatemorph3d
from projection_tools import project_data
import datetime

def on_close(event):
    root.quit()

#%%


#Esta função é utilizada para plotar várias linhas independentes em um gráfico usando a biblioteca matplotlib. Cada linha é representada pelos segmentos definidos na matriz segments.    
def plot_independent_lines(ax, segments, **kwargs):
    for i in range(segments.shape[0]):
        ax.add_line(plt.Line2D(segments[i,:,0], segments[i,:,1], **kwargs))

#%%

#Cria uma caixa de diálogo para selecionar um arquivo de geometria.
window = tk.Toplevel(root)
window.withdraw()
window.attributes("-topmost", True)
window.update()

file_path = tk.filedialog.askopenfilename(title='Open geometry file',
                                          filetypes=(('MATLAB files', '*.mat'),
                                                      ('Numpy arrays', '*.npy'),
                                                      ('All files', '*.*')))
window.destroy()

# file_path = 'C:/Users/Joao Salinet/Downloads/silhs1.mat'

data_path = os.path.dirname(file_path)
#%%
#Cria uma malha triangular a partir dos dados lidos e calcula a matriz de conectividade dos triângulos para a malha
if file_path.endswith('.mat'):
    # Carrega o arquivo .mat usando a biblioteca scipy
    matfile = sio.loadmat(file_path)
    # Extrai os dados dos vértices e faces da malha
    vertices = matfile["vertices"]
    faces = matfile["faces"]
    # Normaliza as faces subtraindo o valor mínimo para garantir que estejam em um intervalo válido
    faces -= faces.min()

# Se o arquivo tiver a extensão '.npy'
elif file_path.endswith('.npy'):
    # Carrega o arquivo .npy usando a biblioteca numpy
    file = np.load(file_path, allow_pickle=True)
    # Extrai os dados dos vértices e faces da malha do objeto carregado
    vertices = file.item()["vertices"]
    faces = file.item()["faces"]
    # Normaliza as faces subtraindo o valor mínimo para garantir que estejam em um intervalo válido
    faces -= faces.min()
    
# Se a extensão não for nem '.mat' nem '.npy'
else:
    # Imprime uma mensagem indicando que o tipo de arquivo não é suportado
    print("File type not supported")
    
def fix_normals(vertices, faces):
    norm_vert = vertices - vertices.mean(axis=0, keepdims=True)
    norm_vert /= np.linalg.norm(norm_vert, axis=1, keepdims=True) + 1e-12
    
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(vertices)
    mesh.triangles = o3d.utility.Vector3iVector(faces)
    mesh.compute_vertex_normals()

    vert_normals = np.asarray(mesh.vertex_normals)

    if np.mean(np.sum(norm_vert * vert_normals, axis=1)) < 0:
        print("Correcting face normals")
        faces = faces[:,(0, 2, 1),]
    
    return faces
        
        
# Cria cópias dos dados originais de vértices e faces para uso posterior
all_vertices = [vertices.copy(),]
all_faces = [faces.copy(),]

all_faces[0] = fix_normals(all_vertices[0][:,:3], all_faces[0])

# Calcula a conectividade das faces em alta resolução usando a função personalizada 'get_connectivity'
HR_conn = get_connectivity(all_faces[0], autoindex=True)

# Cria um objeto de malha triangular usando a biblioteca Open3D
mesh = o3d.geometry.TriangleMesh()

# Atribui os vértices e faces à malha
mesh.vertices = o3d.utility.Vector3dVector(vertices)
mesh.triangles = o3d.utility.Vector3iVector(faces)

# Calcula as normais dos vértices para melhorar a iluminação e renderização
mesh.compute_vertex_normals()

# Simplifica a malha usando a decimação quadrica com o objetivo de reduzir para cerca de 1200 vértices
ntriangles = [20000, 10000, 2500, 1200] 
for ntri in ntriangles:
    mesh = mesh.simplify_quadric_decimation(ntri)

    # Atualiza os arrays 'vertices' e 'faces' com os dados da malha simplificada
    all_vertices.append(np.asarray(mesh.vertices))
    all_faces.append(np.asarray(mesh.triangles))
    
for i, vert in enumerate(all_vertices):
    all_vertices[i] = np.hstack((vert, np.ones((vert.shape[0], 1))))
    all_faces[i] = fix_normals(vert[:,:3], all_faces[i])

    
LR_indices = [find_HR_idx(all_vertices[0], vert) for vert in all_vertices[1:]]

vertices = all_vertices[-1].copy()
vertices[:,:3] -= all_vertices[-1][:,:3].mean(axis=0, keepdims=True)
faces = all_faces[-1].copy()

# Calcula a conectividade das faces da malha simplificada usando a função personalizada 'get_connectivity'
conn = get_connectivity(faces, autoindex=True)

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

def open_calibration(block_num):
    filename = fd.askopenfilename(title="Select the calibration file for the " + camnames[block_num], filetypes=(("Numpy Arrays", "*.npy"),))
    if filename:
        calib_entry[block_num]["state"] = 'normal'
        calib_entry[block_num].delete(0, tk.END)
        calib_entry[block_num].insert(0, filename)
        calib_entry[block_num]["state"] = 'disabled'
        
        calib_files[block_num] = np.load(filename, allow_pickle=True)
    
    
    
    
def turn_image(block_num):
    datas[block_num] = datas[block_num].transpose(0, 2, 1)[:,:,::-1]
    
def flip_image(block_num):
    datas[block_num] = datas[block_num][:,::-1,:]
        
    
def confirm_selection():
    if not any([d is None for d in datas]):
        window.destroy()
        root.quit()
    else:
        print("Please load all views before confirm")
   
 
# Create the main window
window = tk.Toplevel(root)
window.bind("<Destroy>", on_close)
window.title("Load Signal Files")
window.attributes('-topmost', 1)
# window.geometry("800x700")

# Create three blocks
camnames = ["camera C", "camera A", "camera B"]
files = [None, None, None]
calib_files = [None, None, None]
datas = [None, None, None]
paths = [ None, None, None]

blocks = []
entry = []
calib_entry = []
buttons = []
calib_buttons = []
axes = []
canvas = []
selected_options = []
dropdowns = []
turners = []
flippers = []

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
    canvas[i].get_tk_widget().grid(row=0, column=0, columnspan=2)
    
    turners.append(tk.Button(img_frame, text = "⭮", command = lambda i=i: [turn_image(i), update_image(i)], state="disabled", width=10))
    turners[i].grid(row=1, column=0)
    flippers.append(tk.Button(img_frame, text = " ⏛ ", command = lambda i=i: [flip_image(i), update_image(i)], state="disabled", width=10))
    flippers[i].grid(row=1, column=1)
    
    
    text_label = tk.Label(block_frame, text="Calibration File")
    text_label.pack()
    
    # Create an entry for displaying the selected file name
    calib_entry.append(tk.Entry(block_frame, width=40, justify='right',state='disabled'))
    calib_entry[i].pack()

    # Create a button to open a file dialog
    calib_buttons.append(tk.Button(block_frame, text="Select Calibration File", command=lambda i=i: open_calibration(i)))
    calib_buttons[i].pack()
    
confirm_button = tk.Button(window, text="Confirm", command=confirm_selection)
confirm_button.pack(side=tk.RIGHT)

window.attributes('-topmost', 0)

# Start the Tkinter event loop
root.mainloop()

#%%
# Obtém as dimensões (altura, largura) da primeira matriz em 'datas'
t, h, w = datas[0].shape

# Obtém as alturas das matrizes de calibração na lista 'calib_files'
calib_h = [arr.item()["base_shape"][0] for arr in calib_files]

# Obtém as larguras das matrizes de calibração na lista 'calib_files'
calib_w = [arr.item()["base_shape"][1] for arr in calib_files]

# Calcula a matriz pseudo-inversa da calibração para cada matriz de calibração em 'calib_files'
calib_mat = [np.linalg.pinv(arr.item()["calibration"]) for arr in calib_files]
calib_mat = [np.hstack((np.array([[arr[0, 0]], [0], [0], [0]]), 
                        np.vstack((np.zeros((1, 3)), arr)))) for arr in calib_mat]

# Calcula a razão entre os tamanhos da imagem atual dos sinais e a imagem utilizada
# na calibração. Normalmente a variação de tamanhos se deve à aplicação de bins,
bin_factor = [int(np.mean((ch/h, cw/w),)) for ch, cw in zip(calib_h, calib_w)]


# Lista com as imagens mostradas na interface
rgb = False
if datas[0].shape[0] in [3, 4]:
    def prep(img):
        img = img.copy()
        p5, p95 = np.percentile(img, [5, 95])
        img[img < p5] = p5
        img[img > p95] = p95
        img -= p5
        img /= p95 - p5
        return img.transpose(1,2,0)
    
    bg = [prep(data[:3])/2 + 0.25 for data in datas]
    rgb = True
else:
    bg = [data[0] for data in datas]

    
# Os ângulos das câmeras B, A e C, respectivamnte
angles = [-np.pi*2/3, 0, np.pi*2/3]

# Definição de uma matriz de rotação 4x4
final_rotmat = np.array([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])

# Criação de uma janela de interface gráfica usando tkinter
window = tk.Toplevel(root)
# window.lift()
window.attributes('-topmost', 1)

# Definição de uma função chamada "get_normal_mask" que recebe coordenadas de vértices,
# conectividade, triângulos e um ângulo como parâmetros
def get_normal_mask(vertices, connectivity, triangles, ang):
    ang = - ang
    # Criação de uma matriz de rotação usando o ângulo dado
    rotmat = np.array([[np.cos(ang), -np.sin(ang), 0, 0],
                       [np.sin(ang),  np.cos(ang), 0, 0],
                       [          0,            0, 1, 0], 
                       [          0,            0, 0, 1]])
    
    # Aplicação da matriz de rotação aos vértices para obter projeções
    projvert = np.matmul(vertices, rotmat.T)
    
    # Criação de um objeto de malha tridimensional usando os vértices projetados e triângulos
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(projvert[:,:3])
    mesh.triangles = o3d.utility.Vector3iVector(triangles)
    mesh.compute_vertex_normals()
    mesh.compute_triangle_normals()

    face_normals = np.array(mesh.triangle_normals)[:, 0]
    
    vert_normals = np.array(mesh.vertex_normals)
    centered = projvert[:,:3] - projvert[:,:3].mean(axis=0, keepdims=True)
    centered /= np.linalg.norm(centered, axis=1, keepdims=True)
    if np.sum(centered * vert_normals, axis=1).mean() < 0:
        face_normals *= -1
    
    # Criação de uma máscara inicialmente preenchida com zeros
    mask = np.zeros(vertices.shape[0], bool)
    
    # Marcação de vértices com base na orientação das faces da malha
    mask[np.unique(triangles[face_normals > 0])] = True
    
    # Criação de uma máscara de vértices com orientação negativa
    neg_mask = np.zeros(vertices.shape[0], bool)
    neg_mask[np.unique(triangles[face_normals <= 0])] = True
    
    # Criação de uma máscara de borda combinando as máscaras positiva e negativa
    border_mask = mask & neg_mask
    
    # Cálculo da matriz de conectividade de borda
    border_conn = connectivity * border_mask.reshape(1, -1) * border_mask.reshape(-1, 1)
    # Obtenção dos índices da matriz de conectividade de borda
    border_idx = np.vstack(np.where(border_conn)).T    
    
    # Aplicação de operações de erosão e dilatação morfológica à máscara
    mask = erodemorph3d(mask, connectivity, r=2)
    mask = dilatemorph3d(mask, connectivity, r=4)
    mask = erodemorph3d(mask, connectivity, r=2)
    
    # Retorno da máscara resultante, das projeções e dos índices de borda
    return mask, projvert, border_idx

# Lista de rótulos para câmeras
labels = ["Camera C", "Camera A", "Camera B"]
# Listas vazias para armazenar informações
frame = []
s = [1/fact for fact in bin_factor]
dx = [w//2] *3
dy = [h//2] *3

normal_mask = [np.ones(vertices.shape[0], dtype=bool)] *3
rotvert = [None, None, None]
border = [None, None, None]
button_up = []
button_dw = []
button_lf = []
button_rt = []
button_lg = []
button_sm = []

# Função para atualizar um gráfico na interface 
def update_plot(ff):
    global dx           # Coordenadas x do deslocamento
    global dy           # Coordenadas y do deslocamento
    global rotvert      # Vértices rotacionados
    global border       # Índices da borda
    global s            # Fatores de escala
    global normal_mask
    global rgb
    figure = plt.Figure(figsize=(3, 6), dpi=100)             # Criação da figura do gráfico
    ax = figure.add_subplot(111)                            # Adiciona um subplot à figura
    chart_type = FigureCanvasTkAgg(figure, frame[ff])       # Criação de um gráfico na interface
    chart_type.get_tk_widget().grid(row=0, column=0)        # Coloca o gráfico na interface
    if rgb:
        ax.imshow(bg[ff])                          # Exibe uma imagem de fundo no gráfico     
    else:
        ax.imshow(bg[ff], cmap="gray")                          # Exibe uma imagem de fundo no gráfico     
    plot_independent_lines(ax, rotvert[ff][border[ff]][:,:,1:3]*s[ff] + np.array([[[dx[ff], dy[ff]]]]), color='red')         ## Plota linhas independentes no gráfico
    ax.set_xlim(0, bg[ff].shape[1])                         # Define limites do eixo x
    ax.set_ylim(bg[ff].shape[0], 0)                         # Define limites do eixo y
    ax.set_xticks([])                                       # Remove as marcações no eixo x
    ax.set_yticks([])                                       # Remove as marcações no eixo y
    
    # ax2 = figure.add_subplot(212)
    # ax2.plot(vertices[normal_mask[ff],1], vertices[normal_mask[ff],0], 'r.')
    # ax2.plot(vertices[~normal_mask[ff],1], vertices[~normal_mask[ff],0], 'k.')
    
    window.update()                                           # Atualiza a interface

# Função para criar matriz de translação    
def translation_matrix(shift):

    transmat = np.array([[1, 0, 0, shift[0]],
                         [0, 1, 0, shift[1]],
                         [0, 0, 1, shift[2]],
                         [0, 0, 0,        1]])
   
    return  transmat

# Função para criar matriz de escala    
def scale_matrix(factor):
    
    scalemat = np.array([[factor[0], 0, 0, 0],
                         [0, factor[1], 0, 0],
                         [0, 0, factor[2], 0],
                         [0, 0,         0, 1]])
   
    return scalemat
    
# Função para criar matriz de rotação
def rotation_matrix(axis, angle):
    
    # Step 1: Normalize the vector
    axis /= np.linalg.norm(axis)
    
    axx, axy, axz = axis

    # Step 2: Compute the skew-symmetric matrix K
    K = np.array([[0, -axz, axy, 0],
                  [axz, 0, -axx, 0],
                  [-axy, axx, 0, 0], 
                  [   0,   0, 0, 0]])

    # Step 3: Compute the rotation matrix R using the Rodrigues' rotation formula
    return np.identity(4) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
   
# Função para aplicar uma transformação às coordenadas dos vértices
def apply_transform(M):
    global vertices, final_rotmat
    center = vertices.mean(axis=0, keepdims=True)[:,:3]
    vertices[:,:3] -= center
    vertices = np.matmul(vertices, M.T)
    vertices[:,:3] += center
    
    final_rotmat = np.matmul(final_rotmat, M.T)
    
    update_proj()                                   # Chama a função para atualizar a projeção

# Função para atualizar a visualização 3D
def update_rot():
    global vertices, faces
    figureR = plt.Figure(figsize=(5,5), dpi=100)        # Criação da figura 3D
    ax = figureR.add_subplot(111, projection='3d')      # Adiciona um subplot 3D à figura
    chart_type = FigureCanvasTkAgg(figureR, framemesh)  # Criação do gráfico 3D na interface
    chart_type.get_tk_widget().grid(row=0, column=0)    # Coloca o gráfico na interface
    ax.plot_trisurf(vertices[:,0], vertices[:,1], -vertices[:,2], triangles=faces, color="red")
    maxr = np.max(vertices.max(axis=0) - vertices.min(axis=0))/2    # Calcula um raio máximo
    cx, cy, cz, _ = (vertices.max(axis=0) + vertices.min(axis=0))/2    # Calcula o centro
    ax.set_xlim([cx-maxr, cx+maxr])                     # Define limites no eixo x
    ax.set_ylim([cy-maxr, cy+maxr])                     # Define limites no eixo y
    ax.set_zlim([-cz-maxr, -cz+maxr])                   # Define limites no eixo z
    ax.set_axis_off()
        
    window.update()                                       # Atualiza a interface
    
# Função para atualizar a projeção dos modelos 3D nas diferentes câmeras
def update_proj():
    global vertices, angles, normal_mask, rotvert, border, faces, conn
    
    # Loop para cada uma das três câmeras
    for ff in range(3):
        # Obtém a máscara de normais, os vértices rotacionados e a borda
        normal_mask[ff], rotvert[ff], border[ff] = get_normal_mask(vertices, conn, faces, angles[ff])
        # Aplica a matriz de calibração aos vértices rotacionados
        rotvert[ff] = np.matmul(rotvert[ff], calib_mat[ff].T)
        update_plot(ff)                                 # Atualiza o gráfico na interface
        
    # update_rot()


# balloon = tk.tix.Balloon(window)
# Criação do rótulo de quadro para a rotação
framemesh = tk.LabelFrame(window, text="Rotation", padx=2, pady=2)
framemesh.grid(row=0, column=0, rowspan=2)

# Criação de um rótulo de quadro para os botões
button_frame = tk.LabelFrame(window)
button_frame.grid(row=1, column=1, columnspan=3, pady=10)

# Criação da matrizes de transformação e botões às quais estão associadas
M = rotation_matrix(np.array([0., 0, 1]), +np.pi/20)
button_cw = tk.Button(button_frame, text='↷', command=lambda M=M: apply_transform(M), width=10)
# balloon.bind_widget(button_cw, balloonmsg="Rotate the heart clock-wise")
M = rotation_matrix(np.array([0., 0, 1]), -np.pi/20)
button_cc = tk.Button(button_frame, text='⤻', command=lambda M=M: apply_transform(M), width=10)
M = translation_matrix(-0.1 * np.array([0, 0, 1]))
button_up = tk.Button(button_frame, text='⇧', command=lambda M=M: apply_transform(M), width=10)
M = translation_matrix(0.1 * np.array([0, 0, 1]))
button_dw = tk.Button(button_frame, text='⇩', command=lambda M=M: apply_transform(M), width=10)
M = scale_matrix(1 + 0.01 * np.array([1, 1, 1]))
button_lg = tk.Button(button_frame, text="Expand", command=lambda M=M: apply_transform(M), width=10)
M = scale_matrix(1/ (1 + 0.01 * np.array([1, 1, 1])))
button_sm = tk.Button(button_frame, text="Shrink", command=lambda M=M: apply_transform(M), width=10)


button_cw.grid(row=0, column=0)
button_cc.grid(row=1, column=0)
button_up.grid(row=0, column=1)
button_dw.grid(row=1, column=1)
button_lg.grid(row=0, column=2)
button_sm.grid(row=1, column=2)

update_rot()


# Loop para cada uma das três câmeras
for f in range(3):
    # Criação do rótulo de quadro para cada câmera com seu respectivo nome
    frame.append(tk.LabelFrame(window, text=labels[f], padx=2, pady=2))
    frame[f].grid(row=0, column=f+1)

    # Obtém a máscara de normais, os vértices rotacionados e a borda
    normal_mask[f], rotvert[f], border[f] = get_normal_mask(vertices, conn, faces, angles[f])
    # Aplica a matriz de calibração aos vértices rotacionados
    rotvert[f] = np.matmul(rotvert[f], calib_mat[f].T)
    
    # Criação de um rótulo de quadro para os botões específicos de cada câmera
    button_frame = tk.LabelFrame(frame[f])
    button_frame.grid(row=1, column=0)
    
    # Criação da matriz de translação, rotação e escala para cada botão
    M = translation_matrix(-0.1 * np.array([np.cos(angles[f] + np.pi/2), np.sin(angles[f] + np.pi/2), 0]))
    button_lf = tk.Button(button_frame, text='⇦', command=lambda M=M: apply_transform(M), width=5)
    M = translation_matrix(0.1 * np.array([np.cos(angles[f] + np.pi/2), np.sin(angles[f] + np.pi/2), 0]))
    button_rt = tk.Button(button_frame, text='⇨', command=lambda M=M: apply_transform(M), width=5)
    M = rotation_matrix(np.array([np.cos(angles[f]), np.sin(angles[f]), 0]), -np.pi/40)
    button_cc = tk.Button(button_frame, text='⟲', command=lambda M=M: apply_transform(M), width=5)
    M = rotation_matrix(np.array([np.cos(angles[f]), np.sin(angles[f]), 0]), +np.pi/40)
    button_cw = tk.Button(button_frame, text='⟳', command=lambda M=M: apply_transform(M), width=5)
    M = scale_matrix(1 + 0.05 * np.array([np.sin(angles[f]), np.cos(angles[f]), 0]))
    button_lg = tk.Button(button_frame, text='|<>|', command=lambda M=M: apply_transform(M), width=5)
    M = scale_matrix(1/ (1 + 0.05 * np.array([np.sin(angles[f]), np.cos(angles[f]), 0])))
    button_sm = tk.Button(button_frame, text='>||<', command=lambda M=M: apply_transform(M), width=5)
    
    # Definição dos layouts dos botões
    button_lf.grid(row=0, column=0)
    button_rt.grid(row=0, column=5)
    button_cc.grid(row=0, column=1)
    button_cw.grid(row=0, column=4)
    button_lg.grid(row=0, column=2)
    button_sm.grid(row=0, column=3)


update_proj()           # Atualiza a projeção

window.attributes('-topmost', 0)
window.bind("<Destroy>", on_close)
root.mainloop()         # Inicia o loop principal da interface

#%%



window = tk.Toplevel(root)              # Criação da instância da interface gráfica
window.bind("<Destroy>", on_close)

labels = ["Camera C", "Camera A", "Camera B"]                   # Nomes das câmeras
frame = []                              # Lista para armazenar os quadros da interface

button_up = []
button_dw = []
button_lf = []
button_rt = []
button_lg = []
button_sm = []

# Função para atualizar o gráfico na interface
def update_plot(ff):
    global dx           # Coordenadas x do deslocamento
    global dy           # Coordenadas y do deslocamento
    global rotvert      # Vértices rotacionados
    global border       # Índices da borda
    global s            # Fatores de escala
    global rgb
    figure = plt.Figure(figsize=(3, 5), dpi=100)             # Criação da figura do gráfico
    ax = figure.add_subplot(111)                            # Adiciona um subplot à figura
    chart_type = FigureCanvasTkAgg(figure, frame[ff])       # Criação de um gráfico na interface
    chart_type.get_tk_widget().grid(row=1, column=1)        # Coloca o gráfico na interface
    if rgb:
        ax.imshow(bg[ff])                          # Exibe uma imagem de fundo no gráfico     
    else:
        ax.imshow(bg[ff], cmap="gray")                          # Exibe uma imagem de fundo no gráfico     

    plot_independent_lines(ax, rotvert[ff][border[ff]][:,:,1:3]*s[ff] + np.array([[[dx[ff], dy[ff]]]]), color='red')         ## Plota linhas independentes no gráfico
    ax.set_xlim(0, bg[ff].shape[1])                         # Define limites do eixo x
    ax.set_ylim(bg[ff].shape[0], 0)                         # Define limites do eixo y
    ax.set_xticks([])                                       # Remove as marcações no eixo x
    ax.set_yticks([])                                       # Remove as marcações no eixo y
    

    
    window.update()   

def upbutton(ff):
    global dy
    
    dy[ff] -= 1
    
    update_plot(ff)


def dwbutton(ff):
    global dy
    
    dy[ff] += 1
    
    update_plot(ff)
    
def lfbutton(ff):
    global dx
    
    dx[ff] -= 1

    update_plot(ff)


def rtbutton(ff):
    global dx
    
    dx[ff] += 1

    update_plot(ff)

    
def lgbutton(ff):
    global s
    
    s[ff] *= 1.05

    update_plot(ff)

    
def smbutton(ff):
    global s
    
    s[ff] /= 1.05
 
    
    update_plot(ff)
    
for f in range(3):
    frame.append(tk.LabelFrame(window, text=labels[f], padx=2, pady=2))
    frame[f].grid(row=0, column=f)
    
    update_plot(f)
    
    
    button_up.append(tk.Button(frame[f], text='▲', command=lambda f=f: upbutton(f), width=10))
    button_dw.append(tk.Button(frame[f], text='▼', command=lambda f=f: dwbutton(f), width=10))
    button_lf.append(tk.Button(frame[f], text='◄', command=lambda f=f: lfbutton(f), height=10))
    button_rt.append(tk.Button(frame[f], text='►', command=lambda f=f: rtbutton(f), height=10))
    button_lg.append(tk.Button(frame[f], text='◄ ►', command=lambda f=f: lgbutton(f), width=5))
    button_sm.append(tk.Button(frame[f], text='► ◄', command=lambda f=f: smbutton(f), width=5))
    
    button_up[f].grid(row=0, column=1)
    button_dw[f].grid(row=2, column=1)
    button_lf[f].grid(row=1, column=0)
    button_rt[f].grid(row=1, column=2)
    button_lg[f].grid(row=2, column=2)
    button_sm[f].grid(row=2, column=0)

window.mainloop()         

#%%


final_rotmat = np.matmul(np.linalg.pinv(all_vertices[-1]), vertices)

for i in range(len(all_vertices)):
    all_vertices[i] = np.matmul(all_vertices[i], final_rotmat)

#%%
# Inicialização das listas para armazenar projeções e máscaras de projeção
projections = []
projections_mask = []

# Iteração sobre cada uma das três câmeras
for f in range(3):
    # Cálculo da máscara de normais, vértices rotacionados e descarte dos índices da borda
    normal_mask[f], rotvert[f], _ = get_normal_mask(all_vertices[0], HR_conn, all_faces[0], angles[f])
    rotvert[f] = np.matmul(rotvert[f], calib_mat[f].T)

    # Criação das projeções dos vértices rotacionados
    projections.append(rotvert[f][:,1:3] * s[f] + np.array([[dx[f], dy[f]]]) )
    # Conversão da máscara booleana em lista para armazenamento
    projections_mask.append(normal_mask[f].tolist())

for i in range(len(all_vertices)):
    all_vertices[i][:,2] *= -1


#%%
# Inicialização de um array para armazenar as projeções dos sinais projetados
# projected_signals = np.zeros((3, all_vertices[0].shape[0], t))
for i in range(len(all_vertices)):
    all_vertices[i] = all_vertices[i][:,:3]

# Cálculo dos pesos para cada vértice, baseado na distância dos vértices aos seus vizinhos na malha
mesh = o3d.geometry.TriangleMesh()
mesh.vertices = o3d.utility.Vector3dVector(all_vertices[0][:,:3])
mesh.triangles = o3d.utility.Vector3iVector(all_faces[0])
mesh.compute_vertex_normals()
vertex_normals = np.asarray(mesh.vertex_normals)
# Cálculo dos pesos para cada vértice, baseado na distância dos vértices aos seus vizinhos na malha
weights = np.stack([- np.sum(vertex_normals * np.array([[np.cos(theta), np.sin(theta), 0]]), axis=1) for theta in [-2*np.pi/3,0, 2*np.pi/3,]])
weights[weights < np.cos(np.pi * 0.5)] = 0

signals = project_data(datas, projections, projections_mask, all_vertices[0], weights)

# Solicitação de um nome de arquivo para salvar os resultados em formato MAT
filename = ""
while filename == "":
    window = tk.Toplevel(root)
    window.withdraw()
    filename = fd.asksaveasfilename(initialfile = "projected_signals_" + os.path.basename(paths[0])[:-4] + ".mat", title="Save projected signals", filetypes=(("MATLAB files","*.mat"), ("All files", "*.*")))
    window.destroy()
    root.quit()
    
# Criação de um dicionário para armazenar os resultados finais
outfile = {}
# outfile['__header__'] = b'MATLAB 5.0 MAT-file Platform: posix, Created on: Thu Dec 22 19:09:01 2022'
outfile['__header__'] = b'MATLAB 5.0 MAT-file Platform: Python 3, Created on: ' + datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y").encode()
outfile['__version__'] = '1.0'
outfile['__globals__'] = []
outfile["map3d_HR"] = signals
outfile["geometry_HR"] = {"vertices":all_vertices[0], "faces":all_faces[0] +1}
for i, ntri in enumerate(ntriangles):
    outfile[f"geometry_{ntri}"] = {"vertices":all_vertices[i+1], 
                                  "faces":all_faces[i+1] +1, 
                                  "HR_idx":LR_indices[i] +1} 
outfile["dataC"] = datas[0]
outfile["dataA"] = datas[1]
outfile["dataB"] = datas[2]
outfile["dataC_path"] = paths[0]
outfile["dataA_path"] = paths[1]
outfile["dataB_path"] = paths[2]
outfile["projected_vertices"] = projections
outfile["projected_mask"] = projections_mask
outfile["projection_weights"] = weights

# Salvamento dos resultados no arquivo MAT
sio.savemat(filename, outfile)
print("Projection saved in:", filename)



#%%
# from mayavi import mlab
# fig = mlab.figure(bgcolor=(1,1,1))
# mesh = mlab.triangular_mesh(all_vertices[0][:,0], all_vertices[0][:,1], all_vertices[0][:,2], 
#                       all_faces[0], colormap="jet", scalars=signals[:,0], 
#                       )

# # mlab.triangular_mesh(HR_vertices[:,0], HR_vertices[:,1], -HR_vertices[:,2], 
# #                       HR_faces, color = (0,0,0), representation = "mesh",tube_radius=0.005
# #                       )
# # mlab.points3d(vertices[261,0]*1.01, vertices[261,1]*1.01, vertices[261,2]*1.01, 
# #               color=(0,0,0), scale_factor=0.5)
# cbar = mlab.colorbar(title='', orientation='vertical', nb_labels=2)
# cbar.label_text_property.color = (0, 0, 0)
# cbar.title_text_property.color = (0, 0, 0)
# cbar.label_text_property.italic = 0
# cbar.title_text_property.italic = 0
# cbar.scalar_bar.unconstrained_font_size = True
# cbar.label_text_property.font_size = 80

# mlab.view(azimuth=85, elevation=50, distance=12)
# # mlab.savefig("D:/Italo/pesquisas/jimena/figura_cinc2023_2.jpg", size=(160,240))
# mlab.show()
# mlab.close()

#%%


root.destroy()

import vedo
import matplotlib.pyplot as plt

scene = vedo.Plotter().parallel_projection(True)

colormap = "gray"

if signals.shape[1] in (3, 4):
    mesh =  vedo.Mesh([all_vertices[0][:,:3], all_faces[0]],).lighting(None)  
    
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
    # vmin = np.percentile(signals[:,0][signals[:,0] > 0], 1)
    vmin, vmax = np.percentile(signals[:,0][signals[:,0] > 0], [1, 99])

    
    cmap = plt.cm.get_cmap(colormap, 100)

    # Obtém as cores do colormap em um array numpy
    cmap = cmap(np.linspace(0, 1, 100)) 
    cmap[0] = [0.75, 0.75, 0.75, 1]
        
        
    mesh =  vedo.Mesh(
        [all_vertices[0][:,:3], all_faces[0]]
        ).cmap(
            input_cmap=cmap[:,:3], 
            input_array=signals[:,0],
            vmin=vmin, vmax=vmax
            ).lighting(
                "plastic"
                ).add_scalarbar(size=(100, 1000), font_size=30)


if "MEAS_HR_IDX" in matfile.keys():
    for mea in matfile["MEAS_HR_IDX"].keys():
        scene += vedo.Points(vertices[matfile["MEAS_HR_IDX"][mea][0,0][0]].T.tolist(),r=30, c=tuple(np.random.rand(3)))
    
scene += mesh

scene.show(interactive=True)
scene.close()