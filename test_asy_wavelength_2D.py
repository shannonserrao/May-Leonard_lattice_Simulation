import os
import sys
import StringIO
import numpy as np
from scipy.fftpack import fft2 as fft2
from scipy.fftpack import fftshift as fftshift
#import glob2
from scipy.misc import imread
#import pillow
from PIL import Image
import matplotlib
matplotlib.use('agg')

from matplotlib.pyplot import imshow

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.colors import Normalize
#os.chdir('/home/shann87/Control/main_prog/mayleonard_only/test_images_wavelength/') 
def file_read(fname, Lx, Ly, R):
        S_content_array = []
        A_content_array = []
        with open(fname) as f:
                #Content_list is the list that contains the read lines.     
                line1='True'
                for line in f:
                        if line1=='True':
                        	t1 = ' '.join(line.split())
                        	[_, px, py, _]=t1.split()
                        	px=int(px)
                        	py=int(py)
                        	line1='False'
                        else:
                            content=np.reshape(np.fromstring(line, dtype=int,count=3*(1+R)*Lx, sep=' '), ((1+R)*Lx,3))
                            S_content, A_content = content[:-R*Lx, :], content[Lx:, :]
                            #print(S_content.shape, A_content.shape)
                            S_content_array.append(S_content)
                            A_content_array.append(A_content)
                S_content_array=np.asarray(S_content_array)
                A_content_array=np.asarray(A_content_array)
                
                temp=np.all(S_content_array==1,axis=2)
                S_content_array[temp]=0
                temp=np.all(A_content_array==1,axis=2)
                A_content_array[temp]=0
                #print(content_array.shape)
                return S_content_array, A_content_array


#################################
def plot_2dfft_surf(Fx, Fy, abs_img_fft,channel,path):
###############################################################################################

###############                           Display setup            #############################
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    # Make data.
    #X = np.arange(-127, 129, 1)
    #Y = np.arange(-127, 129, 1)
    X, Y = np.meshgrid(Fx, Fy)
    X, Y = np.transpose(X), np.transpose(Y)
    
    # print(len(Fx))
    # print(len(Fy))
    # print(X.shape)
    # print(Y.shape)
    # print(abs_img_fft.shape)
    # # # print(np.maximum(imgfftred))
    # #Plot the surface.
    # surf = ax.plot_surface(X, Y, abs_img_fft, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    # # surf = ax.scatter(X, Y, abs_img_fft, color='green')
    
    # # Customize the z axis.
    # ax.set_xlim(-.5,0.5)
    # ax.set_ylim(-.5,0.5)
    # ax.set_zlim(3400, 4000.01)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    # ax.view_init(45, 45);.
    fft_max, fft_min =np.amax(abs_img_fft), np.amin(abs_img_fft)
    print((fft_max,fft_min))
    colormap = plt.cm.bwr #or any other colormap
    normalize = Normalize(vmin=fft_min, vmax=fft_max)
    #print(abs_img_fft.shape, X.shape, Y.shape)
    area = np.ma.masked_where(abs_img_fft < (abs_img_fft*.75), abs_img_fft)
    area=area/10
    #plt.scatter(X,Y, s=abs_img_fft, c=abs_img_fft, marker='o', alpha=0.5, cmap=colormap,norm=normalize)
    plt.scatter(X,Y, s=area, c=abs_img_fft, marker='o', alpha=0.5, norm=normalize, cmap=colormap)
    plt.xlim(-.25,.25)
    plt.ylim(-.25,.25)
    plt.grid(True)
    plt.xlabel(r'$k_x$',fontsize=25)
    plt.ylabel(r'$k_y$',fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    # plt.title('Wave vector plot for %s' % channel)
    plt.axes().set_aspect('equal')
    # Add a color bar which maps values to colors.
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig(path+channel+'.png',bbox_inches='tight')
    #plt.show()
#####################################################################################################3
def compute_fft(img):
    #################################################################################################3
    #Eliminate the DC Component first 
    img = img - np.mean(img)*np.ones(np.shape(img))
    #imggreen = imggreen - np.mean(imggreen)*np.ones(np.shape(imggreen))
    #imgblue = imgblue - np.mean(imgblue)*np.ones(np.shape(imgblue))
    # Do the Fourier analysis and 
    imgfft = np.fft.fftshift(np.fft.fft2(img))
    #Take the absolute value 
    abs_imgfft=np.abs(imgfft)
    #print(np.abs(imgfftred))
    return abs_imgfft
############################################################################################################
def fft_three_channels(limg):
    #################################################################3
    ##############3 Seperate  in case of two types of images
    # limg = ppm_array[:,:,:] #cropping the pixels
    # rimg = []
    # ##################################################################################################
    # Separate into red/green/blue channels
    imgred = limg[:,:,0]
    imggreen = limg[:,:,1]
    imgblue= limg[:,:,2]
    #################################################################################################3
    abs_imgfftred=compute_fft(imgred)
    abs_imgfftgreen=compute_fft(imggreen)
    abs_imgfftblue=compute_fft(imgblue)
    return abs_imgfftred, abs_imgfftgreen, abs_imgfftblue
###########################################################################################################
###########################################################################################################
def freq_setup(Fs_x=1, Fs_y=1, M=256, N=256):
    dx = 1.0/Fs_x     # centimeters per pixel
    dy = 1.0/Fs_y

    # X = dx*np.linspace(0,M-1)               # centimeters
    # Y = dy*np.linspace(0,N-1)               #   This should be 256 by 256 for a square matrix of size 256 x 256

    dFx = float(Fs_x)/M              # cycles per centimeter
    dFy = float(Fs_y)/N

    twopi=2.0*3.14128
    #Fxt= np.linspace(float(-Fs_x)/2, float(Fs_x)/2 - dFx, num=M)
    #Fyt = np.linspace(float(-Fs_y)/2, float(Fs_y)/2 - dFy, num=N)

    Fx=[-1.0*(Fs_x)/2]
    for i in range(M-1):
        Fx.append(Fx[i]+float(1.0)/M)

    Fx = [i * twopi for i in Fx]

    Fy=[-1.0*(Fs_y)/2]
    for i in range(N-1):
        Fy.append(Fy[i]+float(1.0)/N)

    Fy = [i * twopi for i in Fy]
    # print('Fx,Fy', len(Fx),len(Fy))
    return Fx, Fy 
###########################################################################################################

###########################################################################################################
def find_peaks(abs_fft, Fx, Fy, peaks=4):
    ### Find the specified number of peak in the DFT fourier spectrum and then report their wavelengths.
    # Partition the array into the peaks highest peaks and rest
    ind = np.argpartition(abs_fft, -peaks, axis=None)[-peaks:]
    # Then unravel them for the 2d array
    twod_ind=np.unravel_index(ind[::-1], abs_fft.shape)
    # then resort them in descending order
    ind=ind[np.argsort(abs_fft[twod_ind])[::-1]]
    twod_ind=np.unravel_index(ind[::-1], abs_fft.shape)
    peak_values=abs_fft[twod_ind]
    # D2ind=np.asarray(twod_ind)
    Fx=np.asarray(Fx)
    Fy=np.asarray(Fy)
    # print(type(twod_ind[0][0]))
    # print(type(Fx))
    # print(Fx[:10])
    peakfreq_x=Fx[twod_ind[0]]
    peakfreq_y=Fy[twod_ind[1]]
    peakfreq=np.sqrt(np.square(peakfreq_x)+np.square(peakfreq_x))
    wlen_x=[]
    wlen_y=[]
    for elem in peakfreq_x:
        if elem==0.0:
            wlen_x.append('inf')
        else :
            wlen_x.append(1.0/elem)
    for elem in peakfreq_y:
        if elem==0.0:
            wlen_y.append('inf')
        else :
            wlen_y.append(1.0/elem)    
    wlen=zip(wlen_x,wlen_y)
    return peakfreq, peak_values, wlen

def main():
    ###########################################################################################################
    ################              MAIN PROGRAM   ###########################################################
    ###########################################################################################################

    #ppm_array=file_read('S_256_0.2_0.2_0.8_RPS0_000006500.ppm')
    #ppm_array=file_read('S_256_0.2_0.2_0.8_RPS0_000010900.ppm')
    #ppm_array=file_read('S_256_0.5_0.5_0.5_RPS0_000038300.ppm')


    ppm_array=file_read('S_256_0.2_0.2_0.8_RPS0_000010900.ppm')



    abs_imgfftred, abs_imgfftgreen, abs_imgfftblue =fft_three_channels(ppm_array)

    ###################################################################################################3

    #########    Set up the scale of the wavelengths to be displayed ##############################

    print('The wavelength of the images are :')

    peaks=4
    Fs_x = 1 # sampling freq in x and y 
    Fs_y = 1 # Pixels per distance

    #print(np.shape(imgred))
    (M, N) = np.shape(ppm_array[:,:,0]) 
    #ind = np.argpartition(a, -4)[-4:])      # No of  pixels / datapoints
    Fx, Fy = freq_setup(Fs_x,Fs_y, M, N)

    peak_values, wave_lengths=find_peaks(abs_imgfftred, Fx, Fy, peaks)

    #sys.exit()

    plot_2dfft_surf(Fx, Fy, abs_imgfftred)

if __name__ == "__main__":
    main()



#######################################################################################################
#######################            ROUGH ***********************************************************
#######################################################################################################

#img = imread('S_256_0.2_0.2_0.8_RPS0_000006500.ppm')
#print(img.shape)

#img = Image.open("S_256_0.2_0.2_0.8_RPS0_000006500.ppm")

#print(img.shape)

#imgfftread=np.fft.fftshift(np.fft.fft2(img))

#### Find the maximum of the fourier transform 
# red_max=np.max(np.abs(imgfftred))
# (red_argmaxx,red_argmaxy)=np.argmax(np.abs(imgfftred))
# ind=np.unravel_index(np.argmax(abs_imgfftred,axis=None),abs_imgfftred.shape)
# print(Fx[ind[0]],Fy[ind[1]])
# print(abs_imgfftred[ind])

###################################################################################################
