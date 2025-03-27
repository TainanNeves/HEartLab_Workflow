import struct
import numpy as np
from PIL import Image

def Open_StreamPix(fileName, maxframes=float('inf'), showImg=False):
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
    numFramesToProcess = min(headerInfo['AllocatedFrames'], maxframes)
    
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
