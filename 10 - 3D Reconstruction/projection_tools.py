from scipy.interpolate import griddata
import numpy as np

# Função para subamostragem 2D dos dados usando coordenadas xi
def subsample2D(data, xi):
    
    _, h, w = data.shape
    data = data.copy()
    xi = np.round(xi.copy()).astype(int)
    
    dataout = np.zeros((xi.shape[0], data.shape[0]), dtype=data.dtype)
    
    # Criação de uma máscara para coordenadas xi fora dos limites da imagem
    mask = (xi[:,0] < 0) | (xi[:,1] < 0) | (xi[:,0] >= w) | (xi[:,1] >= h)
    
    # Preenchimento de locais fora da imagem com NaN
    dataout[mask, :] = np.nan
    
    # Preenchimento de locais dentro da imagem com os dados correspondentes
    dataout[~mask] = data[:, xi[~mask,1], 
                             xi[~mask,0]].T
    
        
    return dataout

def project_data(datas, projected_vertices, projections_mask, HR_vertices, weights, no_data_value=float('nan'), phase=False):
    # Conversão das máscaras de projeção em uma matriz NumPy
    is_projected = np.asarray(projections_mask)

    # Verificação se há pelo menos uma projeção para cada vértice
    any_projection = is_projected.sum(axis=0)>0

    t = datas[0].shape[0]

    projected_signals = np.zeros((3, HR_vertices.shape[0], t))
    # Preenchimento dos sinais projetados para cada câmera
    print("Projecting signals...", end="")
    for i in range(3):
        projected_signals[i, is_projected[i]] = subsample2D(datas[i], 
                                                        projected_vertices[i][is_projected[i]])

        print(".", end="")

    # Cálculo dos sinais finais usando os sinais projetados e os pesos   
    weights = weights.copy() * (is_projected>0.5)  * (projected_signals != no_data_value)[...,0]
    if not phase:
        signals = np.sum(projected_signals * 
                        weights[...,np.newaxis],
                        axis=0)/(weights.sum(axis=0).reshape(-1,1) + 1e-12)
    else:
        while np.any(projected_signals > np.pi) or np.any(projected_signals <= -np.pi):
            projected_signals[projected_signals  >  np.pi] -= 2*np.pi
            projected_signals[projected_signals <= -np.pi] += 2*np.pi
        
        projected_signals[1, (projected_signals[1] - projected_signals[0])  >  np.pi] -= 2*np.pi
        projected_signals[1, (projected_signals[1] - projected_signals[0]) <= -np.pi] += 2*np.pi
        projected_signals[2, (projected_signals[2] - projected_signals[0])  >  np.pi] -= 2*np.pi
        projected_signals[2, (projected_signals[2] - projected_signals[0]) <= -np.pi] += 2*np.pi

        signals = np.sum(projected_signals * 
                weights[...,np.newaxis],
                axis=0)/(weights.sum(axis=0).reshape(-1,1) + 1e-12)
        
        while np.any(signals > np.pi) or np.any(signals <= -np.pi):
            signals[signals  >  np.pi] -= 2*np.pi
            signals[signals <= -np.pi] += 2*np.pi


    # Atualização das entradas de projeção que não possuem sinal
    any_projection[np.isnan(signals[:,0])] = False

    # Preenchimento de sinais ausentes usando interpolação
    signals[~any_projection] = griddata(HR_vertices[any_projection], signals[any_projection], HR_vertices[~any_projection], method='nearest')
    print(" Done!")
    
    return signals