import sys
import os
import errno
import numpy as np
from numpy import genfromtxt
import csv
import glob
import fnmatch
from StringIO import StringIO
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.fftpack
from scipy.fftpack import fft
from scipy.interpolate import UnivariateSpline
from scipy.signal import correlate 

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

from shutil import copy2
printdir = "/home/shann87/Control/main_prog/mayleonard_only/std_asy_rates_PMC_2020_noydir_new/"

# ######################################################################################################
# def compute_fwhm_old(x,y, interval):
# 	""" Computes the full width half max and res freq 
# 	"""
# 	y_max=np.max(y)
# 	inv_fmax=1/x[np.argmax(y)]
# 	spline = UnivariateSpline(x, y-y_max/2, s=0)
# 	#r1, r2 = spline.roots()
# 	r =spline.roots()
# 	#print(r)
# 	if len(r) > 1 :
# 		if r[1] <0 or r[0] <0:
# 			print("returning half width half max")
# 		rel_time=1/np.abs(r[1]-r[0])
# 	else:
# 		inv_fmax, rel_time = compute_hwhm(x,y, interval)	
# 	return inv_fmax, rel_time

# def compute_hwhm_old(x,y, interval):
# 	""" Computes the half width half max and res freq 
# 	"""
# 	y_max=np.max(y)
# 	inv_fmax=x[np.argmax(y)]
# 	spline = UnivariateSpline(x, y-y_max/2, s=0)
# 	r1 = spline.roots()
# 	# if r1 <0 or r2 <0:
# 	# 	print("returning half width half max")
# 	rel_time=1/r1

# 	return inv_fmax, rel_time
#######################################################################################################
######################################################################################################
# def compute_fwhm(x,y, interval,std_y):

# 	""" Computes the full width half max and res freq 
# 	"""
# 	f=1E-14
# 	# y=np.log(y+f)
# 	y_max=np.max(y)
# 	fmax=x[np.argmax(y)]
# 	spline = UnivariateSpline(x, y-y_max/2, s=0)
# 	#r1, r2 = spline.roots()
# 	r =spline.roots()
# 	rel_time=0.
# 	# Standard deviation for fmax
# 	sigma_fmax=np.abs(fmax*std_y[np.argmax(y)]/y_max)

# 	#print(r)
# 	if len(r)==1:
# 		print("No roots at full width half max, check data or graph.. ")
# 		fmax, rel_time = compute_hwhm(x,y, interval)
	
# 	# case when data has multiple roots or two roots
# 	elif len(r) > 1:
# 		# r_index=[]

# 		sigma_rel_time=0
# 		candidates=[fmax-r_point for r_point in r]
# 		c_array=np.asarray(candidates)
# 		print(r)
# 		# if len(candidates)==0:
# 		# 	print("Check the y min for this by graphing %f" % (fmax))
# 		f_lower=np.min(c_array[c_array>0])
# 		f_lower=fmax-f_lower

# 		candidates=[r_point-fmax for r_point in r]
# 		# if len(candidates)==0:
# 		# 	print("Check the y min for this by graphing %f" % (fmax))
# 		c_array=np.asarray(candidates)
# 		f_higher=np.min(c_array[c_array>0])
# 		f_higher=fmax+f_higher		

# 		rel_time=1/(f_higher-f_lower)
# 		sigma_rel_time = np.abs(rel_time * std_y[np.argmax(y)]/y_max)
	
# 	# 	inv_fmax, rel_time = compute_hwhm(x,y, interval)	
# 	return fmax, rel_time, sigma_fmax, sigma_rel_time

# def compute_hwhm(x,y, interval, std_y):

# 	""" Computes the half width half max and res freq 
# 	"""
# 	sigma_rel_time=0
# 	f=1E-14
# 	inv_e=0.36787944117
# 	y_min=np.min(y)
# 	y_adj=y-y_min
# 	# y_adj=np.log(y_adj+f)
# 	y_max=np.max(y_adj)
# 	fmax=x[np.argmax(y_adj)]
# 	rel_time=0.
# 	spline = UnivariateSpline(x, y_adj-y_max*inv_e, s=0)
# 	r1 = spline.roots()	
# 	if len(r1) >1:
# 		print('returning all roots ' + str(r1) )
# 	if len(r1) >0:
# 		rel_time=r1[0]
# 		xs = [i for i in range(len(x)) if y_adj[i] > (y_max*inv_e)]
# 		xsmin=np.argmin(x[xs])
# 		sigma_rel_time =np.abs(rel_time * std_y[np.argmin(x[xs])]/(y_max*inv_e))
# 	# else:
# 	# 	y_min=np.min(y)
# 	# 	spline = UnivariateSpline(x, y-y_min, s=0)
# 	# 	r1=spline.roots()
# 	# 	rel_time=r1[0]			
	
# 	# if r1 <0 or r2 <0:
# 	# 	print("returning half width half max")
# 	return fmax, rel_time, sigma_rel_time
#######################################################################################################

# ######################################################################################################
def compute_fwhm(x,y, interval):

	""" Computes the full width half max and res freq 
	"""
	f=1E-14
	# y=np.log(y+f)
	y_max=np.max(y)
	fmax=x[np.argmax(y)]
	spline = UnivariateSpline(x, y-y_max/2, s=0)
	#r1, r2 = spline.roots()
	r =spline.roots()
	rel_time=0.
	#print(r)
	if len(r)==1:
		print("No roots at full width half max, check data or graph.. ")
		fmax, rel_time = compute_hwhm(x,y, interval)
	
	# case when data has multiple roots or two roots
	elif len(r) > 1:
		# r_index=[]

		candidates=[fmax-r_point for r_point in r]
		c_array=np.asarray(candidates)
		print(r)
		# if len(candidates)==0:
		# 	print("Check the y min for this by graphing %f" % (fmax))
		f_lower=np.min(c_array[c_array>0])
		f_lower=fmax-f_lower

		candidates=[r_point-fmax for r_point in r]
		# if len(candidates)==0:
		# 	print("Check the y min for this by graphing %f" % (fmax))
		c_array=np.asarray(candidates)
		f_higher=np.min(c_array[c_array>0])
		f_higher=fmax+f_higher		

		rel_time=1/(f_higher-f_lower)
	
	# 	inv_fmax, rel_time = compute_hwhm(x,y, interval)	
	return fmax, rel_time

def compute_hwhm(x,y, interval):

	""" Computes the half width half max and res freq 
	"""
	f=1E-14
	inv_e=0.36787944117
	y_min=np.min(y)
	y_adj=y-y_min
	# y_adj=np.log(y_adj+f)
	y_max=np.max(y_adj)
	fmax=x[np.argmax(y_adj)]
	rel_time=0.
	spline = UnivariateSpline(x, y_adj-y_max*inv_e, s=0)
	r1 = spline.roots()	
	if len(r1) >1:
		print('returning all roots ' + str(r1) )
	if len(r1) >0:
		rel_time=r1[0]
	# else:
	# 	y_min=np.min(y)
	# 	spline = UnivariateSpline(x, y-y_min, s=0)
	# 	r1=spline.roots()
	# 	rel_time=r1[0]			
	
	# if r1 <0 or r2 <0:
	# 	print("returning half width half max")
	return fmax, rel_time
# #######################################################################################################



#df = pd.DataFrame(columns=['Name','age'])
#runs=len(list)
boolean__init_shape=True

#boxsize
#boxsize=256
boxsize=sys.argv[1]

#parameter rates
pred=sys.argv[2]
rep=sys.argv[3]
diff=sys.argv[4]

#Define the corr array size =(species * length *total_time)
printtime=sys.argv[5]
endtime=sys.argv[6]
#25000=sys.argv[6]
#rootdir="/home/shann87/mayleonardmc/mayleonard_only/new_M28_corrl_ML_spirals256/"
rootdir=sys.argv[7]

#ab=pd.DataFrame(columns=['kappa','Vol','MTE'])
#df = pd.DataFrame(columns=['Name','age'])

spawntime=sys.argv[8]
sitesize=sys.argv[9]
number_loc=sys.argv[10]
interval=sys.argv[12]
#print(interval)
printlattice_every=sys.argv[11]
#print(interval)
folder_name=sys.argv[13]

Nx=int(sys.argv[14])
R=int(sys.argv[15])
asy_fact=float(sys.argv[16])
set_filename=sys.argv[17]
printdir=sys.argv[18]

os.chdir(rootdir+folder_name+'/')
list=os.listdir(rootdir+folder_name+'/')
runs=len(list)
# print(runs)
#Check this  [:,:-1]
tot_time=int((int(endtime)-int(printtime))/float(interval) )+1
#print(tot_time)
#sys.exit()
tot_length=(Nx+R*Nx)
species=3
corr_species=(species+2)*(species+1)/2
print(int(float(tot_length)/2), tot_time)

split=int(Nx)

corr_arr_shape=(corr_species,int(float(tot_length)),tot_time)
density_array_shape=(species,tot_time)

# Define corr_func and density arrays
corr_func=np.zeros(corr_arr_shape)
d_array=np.zeros(density_array_shape)

pattern = 'S_CORR_DEN*.txt'
corr_header=['C_aa','C_ab','C_ac','C_ae','C_bb','C_bc','C_be','C_cc','C_ce','C_ee']
#				0		1	  2		3		4	   5	  6 	 7		8		9
den_header=['D_a','D_b','D_c']



n_of_r=np.zeros((3,tot_length,tot_time))
pattern2='S_CumMinAvg*.txt'

# corr_term_shape=(4,int(float(tot_length)),tot_time)
# corr_term=np.zeros(corr_term_shape)

time_start=4000
time_idx_start=int(time_start/50)

# time_corr_term=np.zeros((6,tot_time))

# for run in list:

# 	os.chdir(rootdir+folder_name+'/'+run+'/')
# 	files=os.listdir('.')
# 	for name in files:
# 		if fnmatch.fnmatch(name,pattern):
# 			timestep=int(name[len(name)-4-9:len(name)-4])
# 			#print(timestep)
# 			if timestep < 12501:
# 				time_index=int(int(timestep)/float(interval))
# 				#print(time_index)
# 				temp=np.genfromtxt(name, dtype=float, skip_footer=int(float(tot_length)), delimiter=',')

# 				d_array[:,time_index]=d_array[:,time_index] + temp[:6]
# 				temp=np.genfromtxt(name, dtype=float, skip_header=1, delimiter =' ')
# 				temp=temp.T
# 				# CHANGES @@@@@@@@@@@@@ TO ADDITIONAL TERMS IN CORRELATION
# 				# corr_term[0,:,time_index]+=temp[0,:]/2 + (temp[1,:]+temp[2,:] +temp[3,:])
# 				# corr_term[1,:,time_index]+=temp[4,:]/2 + (temp[1,:]+temp[5,:] +temp[6,:])
# 				# corr_term[2,:,time_index]+=temp[7,:]/2 + (temp[2,:]+temp[5,:] +temp[8,:])
# 				# corr_term[3,:,time_index]+=temp[9,:]/2 + (temp[3,:]+temp[6,:] +temp[8,:])
				
# 				##
# 				if timestep >= time_start:
# 					# print(time_idx_start)
# 					# print(timestep/50)
# 					# print(time_index)
# 					time_corr_term[:,time_index]=d_array[:,time_idx_start-1]*d_array[:,time_index]
# 					# print(d_array[0,time_idx_start-1]*d_array[0,time_index])
# 					# print(d_array[0,time_index])
# 				else :
# 					time_corr_term[:,time_index]=d_array[:,time_index]
# 				# time_corr_term[1]=d_array[1,0]*d_array[0,time_index]

# 				# CHANGES @@@@@@@@@@@@@
# 				corr_func[:,:,time_index]=temp[:10,:] + corr_func[:,:,time_index]
# 		if fnmatch.fnmatch(name,pattern2):
# 			timestep=int(name[len(name)-4-9:len(name)-4])
# 			#print(timestep)
# 			if timestep < 12501:
# 				time_index=int(int(timestep)/float(interval)) 		
# 				# temp_new=np.genfromtxt(name, dtype=float,skip_footer=384, delimiter=' ')
# 				temp_new=np.genfromtxt(name, dtype=float, delimiter=' ')
				
# 				n_of_r[:,:,time_index]=temp_new.T[:3,:] + n_of_r[:,:,time_index]
# 	# a1=n_of_r[0,0,tot_time-1]
# 	# a2=n_of_r[1,0,tot_time-1]
# 	# a3=n_of_r[2,0,tot_time-1]

# 	# b1=n_of_r[0,128,tot_time-1]
# 	# b2=n_of_r[1,128,tot_time-1]
# 	# b3=n_of_r[2,128,tot_time-1]

# 	# n_of_r[:,:split,:]*=n_of_r[:,0,:][:,np.newaxis,:]/runs
# 	# n_of_r[:,split:,:]*=n_of_r[:,split,:][:,np.newaxis,:]/runs			
# 	# print(time_corr_term[0,time_idx_start:])
# 	# print(runs	
# 	break			
# 	# sys.exit()
# runs=1

	

# # split=int(Nx)
# sym_ar=n_of_r[:,:split,:]
# asym_ar=n_of_r[:,split:,:]

# sym_ar/=runs
# asym_ar/=runs





	
# # split=int(float(Nx))
# corr_func_sym=corr_func[:,:split,:]
# corr_func_asym=corr_func[:,split:,:]

# corr_func_sym/=runs
# corr_func_asym/=runs

# # CHANGES @@@@@@@@@@@@@
# time_corr_term/=runs

# d_array/=runs
# # print(time_corr_term[:,time_idx_start:])
# # time_corr_term-=time_corr_term


# # print(d_array[:,time_idx_start-1])
# # print(d_array[:,:])

# # time_corr_term-=np.multiply((d_array[:,time_idx_start-1][:,np.newaxis]),(d_array[:,:]))
# time_corr_term[0,:]-=np.multiply((corr_func_sym[0,0,time_idx_start-1]),(corr_func_sym[0,0,:]))
# time_corr_term[1,:]-=np.multiply((corr_func_sym[4,0,time_idx_start-1]),(corr_func_sym[4,0,:]))
# time_corr_term[2,:]-=np.multiply((corr_func_sym[7,0,time_idx_start-1]),(corr_func_sym[7,0,:]))
# time_corr_term[3,:]-=np.multiply((corr_func_asym[0,0,time_idx_start-1]),(corr_func_asym[0,0,:]))
# time_corr_term[4,:]-=np.multiply((corr_func_asym[4,0,time_idx_start-1]),(corr_func_asym[4,0,:]))
# time_corr_term[5,:]-=np.multiply((corr_func_asym[7,0,time_idx_start-1]),(corr_func_asym[7,0,:]))


# # print(d_array[:,])
# # sys.exit()

# # corr_term/=runs
# # corr_term*=(corr_term[:,0,:][:,np.newaxis,:]/runs)

# # corr_term_sym=corr_term[:,:split,:]
# # corr_term_asym=corr_term[:,split:,:]


# # avg_corr_s1=np.zeros_like(corr_func_sym[0,:,tot_time-1])
# # avg_corr_s2=np.zeros_like(corr_func_sym[4,:,tot_time-1])
# # avg_corr_s3=np.zeros_like(corr_func_sym[7,:,tot_time-1])
# # avg_corr_as1=np.zeros_like(corr_func_asym[0,:,tot_time-1])
# # avg_corr_as2=np.zeros_like(corr_func_asym[4,:,tot_time-1])
# # avg_corr_as3=np.zeros_like(corr_func_asym[7,:,tot_time-1])

# avg_corr_s=np.zeros((3,Nx))
# avg_corr_as=np.zeros((3,R*Nx))


# avg_ar_s=np.zeros_like(sym_ar[:,:,tot_time-1])
# avg_ar_as=np.zeros_like(asym_ar[:,:,tot_time-1])


# for i in range(9):
# 	avg_corr_s[0,:]+=corr_func_sym[0,:,tot_time-1-i*4]
# 	avg_corr_s[1,:]+=corr_func_sym[4,:,tot_time-1-i*4]
# 	avg_corr_s[2,:]+=corr_func_sym[7,:,tot_time-1-i*4]
# 	avg_corr_as[0,:]+=corr_func_asym[0,:,tot_time-1-i*4]
# 	avg_corr_as[1,:]+=corr_func_asym[4,:,tot_time-1-i*4]
# 	avg_corr_as[2,:]+=corr_func_asym[7,:,tot_time-1-i*4]
# 	avg_ar_s+=sym_ar[:,:,tot_time-1-i*4]
# 	avg_ar_as+=asym_ar[:,:,tot_time-1-i*4]

# # avg_corr_s1/=10
# # avg_corr_s2/=10
# # avg_corr_s3/=10

# # avg_corr_as1/=10
# # avg_corr_as2/=10
# # avg_corr_as3/=10

# avg_corr_s/=10
# avg_corr_as/=10

# avg_ar_s/=10
# avg_ar_as/=10




# avg_ar_ao_s=np.zeros_like(sym_ar[:,:,tot_time-1])
# avg_ar_ao_as=np.zeros_like(asym_ar[:,:,tot_time-1])

# avg_ar_ao_s=avg_ar_s[:,0][:,np.newaxis]*avg_ar_s[:,:]
# avg_ar_ao_as=avg_ar_as[:,0][:,np.newaxis]*avg_ar_as[:,:]











# # Denstity plots
# if not os.path.exists(printdir+set_filename+'/'+folder_name+'/'):
#     os.makedirs(printdir+set_filename+'/'+folder_name+'/')
# #os.mkdir(printdir)
# #os.chdir(printdir+set_filename+'/'+folder_name+'/')

# copy2(rootdir+folder_name+'/run1/Sdatafile.txt', printdir+set_filename+'/'+folder_name+'/')


os.chdir(printdir+set_filename+'/'+folder_name+'/')
#shutil.copy2(rootdir+'/run1/Sdatafile.txt', printdir+'/'+folder_name+'/Sdatafile.txt')

d_array=np.genfromtxt("d_array.csv", delimiter=",")

time_corr_term=np.genfromtxt("time_corr_term.csv", delimiter=",")
avg_corr_s=np.genfromtxt("avg_corr_s.csv",  delimiter=",")
# avg_corr_as=np.genfromtxt("avg_corr_as.csv",  delimiter=",")
avg_ar_ao_s=np.genfromtxt("avg_ar_ao_s.csv", delimiter=",")
# avg_ar_ao_as=np.genfromtxt("avg_ar_ao_as.csv", delimiter=",")

std_time_corr_term=np.genfromtxt("std_time_corr_term.csv", delimiter=",")
std_avg_corr_s=np.genfromtxt("std_avg_corr_s.csv",  delimiter=",")
# std_avg_corr_as=np.genfromtxt("std_avg_corr_as.csv",  delimiter=",")
std_avg_ar_ao_s=np.genfromtxt("std_avg_ar_ao_s.csv", delimiter=",")
# std_avg_ar_ao_as=np.genfromtxt("std_avg_ar_ao_as.csv", delimiter=",")


with open('Sdatafile.txt','a+') as f:
	f.write('red : %s \n'  % (str(np.mean(d_array[0,time_idx_start:])) ))
	f.write('green : %s \n'  % (str(np.mean(d_array[1,time_idx_start:])) ))
	f.write('blue : %s \n'  % (str(np.mean(d_array[2,time_idx_start:])) ))


# # time linspace
# time=np.arange(tot_time)
time=np.arange(tot_time)
N=int(d_array.shape[1])
T=int(interval)


# fig1, (ax1, ax2) = plt.subplots(2, 1)
# fig1.subplots_adjust(hspace=0.5)

# ax1.plot(T*time,d_array[0,:],T*time,d_array[1,:],T*time,d_array[3,:],T*time,d_array[4,:])
# ax1.set_xlabel('time in MC steps')
# ax1.set_ylabel('density')
# ax1.grid(True)
# ax1.set_title('Density_plot, S/A=1:%d, D= %s ' % (R, diff) )
# ax1.legend(('sym_a', 'sym_b', 'asym_a', 'asym_b'), loc='upper right')
# ax1.tick_params(axis='both', which='major', labelsize=10)


# #xf=np.linspace(0.0, 1.0/(2.0*T), ceil(N/2+1))
# xf=np.linspace(0.0, 1.0/(2.0*T), np.ceil(N/2+1))
# # Denstity plots
# da1_f = scipy.fftpack.fft(d_array[0,:])
# db1_f = scipy.fftpack.fft(d_array[1,:])
# #dc1_f = scipy.fftpack.fft(d_array[2,:])
# da2_f = scipy.fftpack.fft(d_array[3,:])
# db2_f = scipy.fftpack.fft(d_array[4,:])


# ax2.plot(xf[1:], 2.0/N * np.abs(da1_f[1:N//2+1]), xf[1:], 2.0/N * np.abs(db1_f[1:N//2+1]), xf[1:], 2.0/N * np.abs(da2_f[1:N//2+1]), xf[1:], 2.0/N * np.abs(db2_f[1:N//2+1]))
# ax2.set_xlabel('frequency')
# ax2.set_ylabel('density')
# ax2.grid(True)
# # ax2.set_title('Density_plot, S/A=1:1, D= %s ' % diff )
# ax2.legend(('sym_a', 'sym_b', 'asym_a', 'asym_b'), loc='upper right')
# ax2.tick_params(axis='both', which='major', labelsize=10)

# fig1.savefig(printdir+set_filename+'/'+folder_name+'/'+'den_plot'+"ss"+sitesize+"nloc"+number_loc+"spw"+spawntime+'.png')

# invfreq1, rel_time1= compute_fwhm(xf[1:],2.0/N * np.abs(da1_f[1:N//2+1]), T)
# invfreq2, rel_time2= compute_fwhm(xf[1:],2.0/N * np.abs(db1_f[1:N//2+1]), T)
# invfreq3, rel_time3= compute_fwhm(xf[1:],2.0/N * np.abs(da2_f[1:N//2+1]), T)
# invfreq4, rel_time4= compute_fwhm(xf[1:],2.0/N * np.abs(db2_f[1:N//2+1]), T)

# with open('Sdatafile.txt','a+') as f:
#  	#f.write('roots are : %s %s '  % (str(r1), str(r2)) )
#  	f.write('inv_resonance freq relaxation time for den plot is : %s %s \n'  % (str(invfreq1), str(rel_time1)) )
#  	f.write('inv_resonance freq relaxation time for den plot is : %s %s \n'  % (str(invfreq2), str(rel_time2)) )
#  	f.write('inv_resonance freq relaxation time for den plot is : %s %s \n'  % (str(invfreq3), str(rel_time3)) )
#  	f.write('inv_resonance freq relaxation time for den plot is : %s %s \n'  % (str(invfreq4), str(rel_time4)) )




##############################################  Time corr changes #################################################################

def std_multiply(A,B,std_A, std_B):
	f=1E-14
	mul_AB = np.abs(np.multiply(A,B))*np.sqrt(np.square(np.divide(std_A,A+f)) + np.square(np.divide(std_B,B+f)) + 2*np.divide(np.cov(A,B), np.multiply(A,B)))
	return mul_AB

def std_add(A,B,std_A, std_B):
	add_AB= np.sqrt(np.multiply(np.square(A),np.square(std_A)) + np.multiply(np.square(B),np.square(std_B)))
	return add_AB
















# fig3.subplots_adjust(hspace=0.5)


# ac_sym_0=correlate(corr_func_sym[0,0,:],corr_func_sym[0,0,:])
# ac_sym_4=correlate(corr_func_sym[4,0,:],corr_func_sym[4,0,:])
# ac_sym_7=correlate(corr_func_sym[7,0,:],corr_func_sym[7,0,:])
# ac_asym_0=correlate(corr_func_asym[0,0,:],corr_func_asym[0,0,:])
# ac_asym_4=correlate(corr_func_asym[4,0,:],corr_func_asym[4,0,:])
# ac_asym_7=correlate(corr_func_asym[7,0,:],corr_func_asym[7,0,:])


# print(type(ac_sym_0[0:tot_time]))
# print(np.array(ac_sym_0[0:tot_time][0]).shape)
# print(T.shape)


# time=np.arange(len(ac_sym_0[tot_time-1:]))

# print(len(ac_sym_0[tot_time-1:]))
# print(len(time))

# sys.exit()
fig0, ax0 = plt.subplots()
# # ax1.plot(T*time,d_array[0,:],T*time,d_array[1,:],T*time,d_array[3,:],T*time,d_array[4,:])
ax0.plot(T*time,d_array[0,:],'r-',T*time,d_array[1,:],'g-',T*time,d_array[2,:],'b-')
# ax1.plot(T*time,ac_sym_0[tot_time-1:],T*time,ac_sym_4[tot_time-1:],T*time,ac_sym_7[tot_time-1:],T*time,ac_asym_0[tot_time-1:],T*time,ac_asym_4[tot_time-1:],T*time,ac_asym_7[tot_time-1:])
ax0.set_xlabel(r't (MCS)', fontsize=25)
ax0.set_ylabel(r' $\langle n_{X}(t) \rangle$', fontsize=25)
ax0.tick_params(axis='x', labelsize=20)
ax0.tick_params(axis='y', labelsize=20)

ax0.grid(True)
ax0.legend((r'$a$', r'$b$', r'$c$'), loc='upper right', fontsize=15)

fig0.savefig(printdir+set_filename+'/'+folder_name+'/'+'plot_density'+'.png',bbox_inches='tight')






fig1, ax1 = plt.subplots()
ax1.errorbar(T*time[time_idx_start:],time_corr_term[0,time_idx_start:].T,yerr=std_time_corr_term[0,time_idx_start:].T,c='r',ls='-')
ax1.errorbar(T*time[time_idx_start:],time_corr_term[1,time_idx_start:].T,yerr=std_time_corr_term[1,time_idx_start:].T,c='g',ls='-')
ax1.errorbar(T*time[time_idx_start:],time_corr_term[2,time_idx_start:].T,yerr=std_time_corr_term[2,time_idx_start:].T,c='b',ls='-')
# ,T*time[time_idx_start:],time_corr_term[1,time_idx_start:].T,yerr=std_time_corr_term[1,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[2,time_idx_start:].T,yerr=std_time_corr_term[2,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[3,time_idx_start:].T,yerr=std_time_corr_term[3,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[4,time_idx_start:].T,yerr=std_time_corr_term[4,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[5,time_idx_start:].T,yerr=std_time_corr_term[5,time_idx_start:].T)
# ax1.errorbar(T*time[time_idx_start:],time_corr_term[3,time_idx_start:].T,yerr=std_time_corr_term[3,time_idx_start:].T)
# ax1.errorbar(T*time[time_idx_start:],time_corr_term[4,time_idx_start:].T,yerr=std_time_corr_term[4,time_idx_start:].T)
# ax1.errorbar(T*time[time_idx_start:],time_corr_term[5,time_idx_start:].T,yerr=std_time_corr_term[5,time_idx_start:].T)


# ax1.plot(T*time,ac_sym_0[tot_time-1:],T*time,ac_sym_4[tot_time-1:],T*time,ac_sym_7[tot_time-1:],T*time,ac_asym_0[tot_time-1:],T*time,ac_asym_4[tot_time-1:],T*time,ac_asym_7[tot_time-1:])
ax1.set_xlabel(r't (MCS)',fontsize=25)
ax1.set_ylabel(r' $C_{XX}(t)$',fontsize=25)
ax1.grid(True)
# ax1.set_title('Density_plot, S/A=1:%d, D= %s ' % (R, diff) )
ax1.legend((r'$a_{s}$', r'$b_{s}$', r'$c_{s}$'), loc='upper right',fontsize=25)
ax1.tick_params(axis='both', which='major', labelsize=20)


# # fig3.savefig(printdir+set_filename+'/'+folder_name+'/'+'acorr_plot'+"ss"+sitesize+"nloc"+number_loc+"spw"+spawntime+'.png')
	
N=int(time_corr_term[:,time_idx_start:].shape[1])
print(N)

tf=np.linspace(0.0, 1.0/(2.0*T), np.ceil(N/2+1))
ds_1_f = scipy.fftpack.fft(time_corr_term[0,time_idx_start:])
ds_2_f = scipy.fftpack.fft(time_corr_term[1,time_idx_start:])
ds_3_f = scipy.fftpack.fft(time_corr_term[2,time_idx_start:])
# da_1_f = scipy.fftpack.fft(time_corr_term[3,time_idx_start:])
# da_2_f = scipy.fftpack.fft(time_corr_term[4,time_idx_start:])
# da_3_f = scipy.fftpack.fft(time_corr_term[5,time_idx_start:])

# # check matrix dimensions
# print(np.shape(time))
# print(np.shape(time_corr_term[0,:].T))
fig2, ax2 = plt.subplots()


# ax2.plot(T*time[time_idx_start:],time_corr_term[0,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[1,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[2,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[3,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[4,time_idx_start:].T,T*time[time_idx_start:],time_corr_term[5,time_idx_start:].T)

ax2.plot(tf[1:], 2.0/N * np.abs(ds_1_f[1:N//2+1]),'r-', tf[1:], 2.0/N * np.abs(ds_2_f[1:N//2+1]),'g-', tf[1:], 2.0/N * np.abs(ds_3_f[1:N//2+1]),'b-')
ax2.set_xlabel(r' f',fontsize=25)
ax2.set_ylabel(r' $C_{XX}(f)$',fontsize=25)
ax2.grid(True)
# ax2.set_title('Density_plot, S/A=1:1, D= %s ' % diff )
ax2.legend((r'$a_{s}$', r'$b_{s}$', r'$c_{s}$'), loc='upper right',fontsize=25)

ax2.tick_params(axis='both', which='major', labelsize=20)

fig1.savefig(printdir+set_filename+'/'+folder_name+'/'+'acorr_plot'+'.png',bbox_inches='tight')
fig2.savefig(printdir+set_filename+'/'+folder_name+'/'+'acorr_freq_plot'+'.png',bbox_inches='tight')

# invfreq1sy, rel_time1sy= compute_hwhm(T*time,ac_sym_0[tot_time-1:], T)
# invfreq2sy, rel_time2sy= compute_hwhm(T*time,ac_sym_4[tot_time-1:], T)
# invfreq3sy, rel_time3sy= compute_hwhm(T*time,ac_sym_7[tot_time-1:], T)
# invfreq1as, rel_time1as= compute_hwhm(T*time,ac_asym_0[tot_time-1:], T)
# invfreq2as, rel_time2as= compute_hwhm(T*time,ac_asym_4[tot_time-1:], T)
# invfreq3as, rel_time3as= compute_hwhm(T*time,ac_asym_7[tot_time-1:], T)

# std_ds_1_f = scipy.fftpack.fft(std_time_corr_term[0,time_idx_start:])
# std_ds_2_f = scipy.fftpack.fft(std_time_corr_term[1,time_idx_start:])
# std_ds_3_f = scipy.fftpack.fft(std_time_corr_term[2,time_idx_start:])

std_ds_1_f = scipy.fftpack.fft(time_corr_term[0,time_idx_start:]+std_time_corr_term[0,time_idx_start:])
std_ds_2_f = scipy.fftpack.fft(time_corr_term[1,time_idx_start:]+std_time_corr_term[1,time_idx_start:])
std_ds_3_f = scipy.fftpack.fft(time_corr_term[2,time_idx_start:]+std_time_corr_term[2,time_idx_start:])

invfreq1sy, rel_time1sy= compute_fwhm(tf[1:],2.0/N * np.abs(ds_1_f[1:N//2+1]), T)
invfreq2sy, rel_time2sy= compute_fwhm(tf[1:],2.0/N * np.abs(ds_2_f[1:N//2+1]), T)
invfreq3sy, rel_time3sy= compute_fwhm(tf[1:],2.0/N * np.abs(ds_3_f[1:N//2+1]), T)


std_invfreq1sy, std_rel_time1sy= compute_fwhm(tf[1:],2.0/N * np.abs(std_ds_1_f[1:N//2+1]), T)
std_invfreq2sy, std_rel_time2sy= compute_fwhm(tf[1:],2.0/N * np.abs(std_ds_2_f[1:N//2+1]), T)
std_invfreq3sy, std_rel_time3sy= compute_fwhm(tf[1:],2.0/N * np.abs(std_ds_3_f[1:N//2+1]), T)

# invfreq1sy, rel_time1sy, std_invfreq1sy, std_rel_time1sy= compute_fwhm(tf[1:],2.0/N * np.abs(ds_1_f[1:N//2+1]), T, std_ds_1_f)
# invfreq2sy, rel_time2sy, std_invfreq2sy, std_rel_time2sy= compute_fwhm(tf[1:],2.0/N * np.abs(ds_2_f[1:N//2+1]), T, std_ds_2_f)
# invfreq3sy, rel_time3sy, std_invfreq3sy, std_rel_time3sy= compute_fwhm(tf[1:],2.0/N * np.abs(ds_3_f[1:N//2+1]), T, std_ds_3_f)
# invfreq1as, rel_time1as= compute_fwhm(tf[1:],2.0/N * np.abs(da_1_f[1:N//2+1]), T)
# invfreq2as, rel_time2as= compute_fwhm(tf[1:],2.0/N * np.abs(da_2_f[1:N//2+1]), T)
# invfreq3as, rel_time3as= compute_fwhm(tf[1:],2.0/N * np.abs(da_3_f[1:N//2+1]), T)


with open('Sdatafile.txt','a+') as f:
	f.write('a symm IRF  , rel time is : %s %s \n'  % (str(invfreq1sy), str(rel_time1sy)) )
	f.write('b symm IRF  , rel time is : %s %s \n'  % (str(invfreq2sy), str(rel_time2sy)) )
	f.write('c symm IRF  , rel time is : %s %s \n'  % (str(invfreq3sy), str(rel_time3sy)) )
	# f.write('a asymm IRF  , rel time is : %s %s \n'  % (str(invfreq1as), str(rel_time1as)) )
	# f.write('b asymm IRF  , rel time is : %s %s \n'  % (str(invfreq2as), str(rel_time2as)) )
	# f.write('c asymm IRF  , rel time is : %s %s \n'  % (str(invfreq3as), str(rel_time3as)) )


std_invfreq1sy, std_rel_time1sy= np.abs(std_invfreq1sy-invfreq1sy), np.abs(std_rel_time1sy-rel_time1sy) 
std_invfreq2sy, std_rel_time2sy= np.abs(std_invfreq2sy-invfreq2sy), np.abs(std_rel_time2sy-rel_time2sy)
std_invfreq3sy, std_rel_time3sy= np.abs(std_invfreq3sy-invfreq3sy), np.abs(std_rel_time3sy-rel_time3sy) 



################################################################################
# std_da_1_f = scipy.fftpack.fft(time_corr_term[3,time_idx_start:]+std_time_corr_term[3,time_idx_start:])
# std_da_2_f = scipy.fftpack.fft(time_corr_term[4,time_idx_start:]+std_time_corr_term[4,time_idx_start:])
# std_da_3_f = scipy.fftpack.fft(time_corr_term[5,time_idx_start:]+std_time_corr_term[5,time_idx_start:])




# std_invfreq1sy, std_rel_time1sy= compute_fwhm(tf[1:],2.0/N * np.abs(std_ds_1_f[1:N//2+1]), T)
# std_invfreq2sy, std_rel_time2sy= compute_fwhm(tf[1:],2.0/N * np.abs(std_ds_2_f[1:N//2+1]), T)
# std_invfreq3sy, std_rel_time3sy= compute_fwhm(tf[1:],2.0/N * np.abs(std_ds_3_f[1:N//2+1]), T)



# std_invfreq1as, std_rel_time1as= compute_fwhm(tf[1:],2.0/N * np.abs(std_da_1_f[1:N//2+1]), T)
# std_invfreq2as, std_rel_time2as= compute_fwhm(tf[1:],2.0/N * np.abs(std_da_2_f[1:N//2+1]), T)
# std_invfreq3as, std_rel_time3as= compute_fwhm(tf[1:],2.0/N * np.abs(std_da_3_f[1:N//2+1]), T)



# std_invfreq1sy, std_rel_time1sy= np.abs(std_invfreq1sy-invfreq1sy), np.abs(std_rel_time1sy-rel_time1sy) 
# std_invfreq2sy, std_rel_time2sy= np.abs(std_invfreq2sy-invfreq2sy), np.abs(std_rel_time2sy-rel_time2sy)
# std_invfreq3sy, std_rel_time3sy= np.abs(std_invfreq3sy-invfreq3sy), np.abs(std_rel_time3sy-rel_time3sy) 
# std_invfreq1as, std_rel_time1as= np.abs(std_invfreq1as-invfreq1as), np.abs(std_rel_time1as-rel_time1as)
# std_invfreq2as, std_rel_time2as= np.abs(std_invfreq2as-invfreq2as), np.abs(std_rel_time2as-rel_time2as) 
# std_invfreq3as, std_rel_time3as= np.abs(std_invfreq3as-invfreq3as), np.abs(std_rel_time3as-rel_time3as) 

#################################################################################









fig3, ax3 = plt.subplots()
fig3.subplots_adjust(hspace=0.5)



L_sym=int(float(Nx))
# L_asy=int(float(R*Nx))
dL=1

dist1=np.arange(L_sym)
# dist2=np.arange(L_asy)
#

# avg_ar_ao_s=avg_ar_s[:,:]
# avg_ar_ao_as=avg_ar_as[:,:]


# print(avg_ar_ao_s.shape)
# print(avg_corr_s1.shape)

# dist1[0]=1
# dist2[0]=1
# dist1=np.log(dist1)
# dist2=np.log(dist2)

# sys.exit()
# avg_corr_as1-avg_ar_ao_as[0,:]
# avg_corr_as2-avg_ar_ao_as[1,:]
# avg_corr_as3-avg_ar_ao_as[2,:]

# avg_corr_array_s=np.zeros((3,L_sym))
# avg_corr_array_as=np.zeros((3,L_asym))

avg_corr_array_s=avg_corr_s-avg_ar_ao_s 
# avg_corr_array_as=avg_corr_as-avg_ar_ao_as

std_avg_corr_array_s=std_add(avg_corr_s,avg_ar_ao_s,std_avg_corr_s,std_avg_ar_ao_s) 
# std_avg_corr_array_as=std_add(avg_corr_as,avg_ar_ao_as,std_avg_corr_as,std_avg_ar_ao_as) 

# ax3.plot(dist1,avg_corr_s1, dist1,avg_corr_s2, dist1,avg_corr_s3, dist2,avg_corr_as1, dist2,avg_corr_as2, dist2,avg_corr_as3)

#ax3.plot(dist1[1:],avg_corr_s[0,1:]-avg_ar_ao_s[0,1:], dist1[1:],avg_corr_s[1,1:]-avg_ar_ao_s[1,1:], dist1[1:],avg_corr_s[2,1:]-avg_ar_ao_s[2,1:], dist2[1:],avg_corr_as[0,1:]-avg_ar_ao_as[0,1:], dist2[1:],avg_corr_as[1,1:]-avg_ar_ao_as[1,1:], dist2[1:],avg_corr_as[2,1:]-avg_ar_ao_as[2,1:])

ax3.errorbar(dist1[1:],avg_corr_s[0,1:]-avg_ar_ao_s[0,1:],yerr=std_avg_corr_array_s[0,1:],c='red', ls='-') #, dist1[1:],avg_corr_s[1,1:]-avg_ar_ao_s[1,1:], dist1[1:],avg_corr_s[2,1:]-avg_ar_ao_s[2,1:], dist2[1:],avg_corr_as[0,1:]-avg_ar_ao_as[0,1:], dist2[1:],avg_corr_as[1,1:]-avg_ar_ao_as[1,1:], dist2[1:],avg_corr_as[2,1:]-avg_ar_ao_as[2,1:])
ax3.errorbar(dist1[1:],avg_corr_s[1,1:]-avg_ar_ao_s[1,1:],yerr=std_avg_corr_array_s[1,1:],c='green', ls='-')
ax3.errorbar(dist1[1:],avg_corr_s[2,1:]-avg_ar_ao_s[2,1:],yerr=std_avg_corr_array_s[2,1:],c='blue', ls='-')

# ax3.errorbar(dist2[1:],avg_corr_as[0,1:]-avg_ar_ao_as[0,1:],yerr=std_avg_corr_array_as[0,1:]) #, dist1[1:],avg_corr_s[1,1:]-avg_ar_ao_s[1,1:], dist1[1:],avg_corr_s[2,1:]-avg_ar_ao_s[2,1:], dist2[1:],avg_corr_as[0,1:]-avg_ar_ao_as[0,1:], dist2[1:],avg_corr_as[1,1:]-avg_ar_ao_as[1,1:], dist2[1:],avg_corr_as[2,1:]-avg_ar_ao_as[2,1:])
# ax3.errorbar(dist2[1:],avg_corr_as[1,1:]-avg_ar_ao_as[1,1:],yerr=std_avg_corr_array_as[1,1:])
# ax3.errorbar(dist2[1:],avg_corr_as[2,1:]-avg_ar_ao_as[2,1:],yerr=std_avg_corr_array_as[2,1:])


#ax3.plot(dist1[1:],np.log(avg_corr_s[0,1:]-avg_ar_ao_s[0,1:]), dist1[1:],np.log(avg_corr_s[1,1:]-avg_ar_ao_s[1,1:]), dist1[1:],np.log(avg_corr_s[2,1:]-avg_ar_ao_s[2,1:]), dist2[1:],np.log(avg_corr_as[0,1:]-avg_ar_ao_as[0,1:]), dist2[1:],np.log(avg_corr_as[1,1:]-avg_ar_ao_as[1,1:]), dist2[1:],np.log(avg_corr_as[2,1:]-avg_ar_ao_as[2,1:]))
# ax2.plot(dist1,avg_ar_ao_s[0,:], dist1,avg_ar_ao_s[1,:], dist1,avg_ar_ao_s[2,:], dist2,avg_ar_ao_as[0,:], dist2,avg_ar_ao_as[1,:], dist2,avg_ar_ao_as[2,:])


# Log version
# ax1.plot(dist1,np.log(avg_corr_s1), dist1,np.log(avg_corr_s2), dist1,np.log(avg_corr_s3), dist2,np.log(avg_corr_as1), dist2,np.log(avg_corr_as2), dist2,np.log(avg_corr_as3))
#

ax3.set_xlabel(r' $r$',fontsize=25)
ax3.set_ylabel(r' $C_{XX}(r)$',fontsize=25)
ax3.grid(True)
# ax3.legend((r'a_s', 'b_s', 'c_s', 'a_as',  'b_as', 'c_as'), loc='upper right')
ax3.legend((r'$a_{s}$', r'$b_{s}$', r'$c_{s}$'), loc='upper right',fontsize=25)
ax3.tick_params(axis='both', which='major', labelsize=20)
# ax3.set_title('Correlation plot, S/A=1:%d D = %s' % (R, diff))
ax3.set_xlim(0, 128)

fig3.savefig(printdir+set_filename+'/'+folder_name+'/'+'EQT_corr_plot'+'.png',bbox_inches='tight')

invfreq1a, rel_time1a= compute_hwhm(dist1[1:],avg_corr_s[0,1:]-avg_ar_ao_s[0,1:], 1)
# invfreq2a, rel_time2a= compute_hwhm(dist2[1:],avg_corr_as[0,1:]-avg_ar_ao_as[0,1:], 1)

invfreq1b, rel_time1b= compute_hwhm(dist1[1:],avg_corr_s[1,1:]-avg_ar_ao_s[1,1:], 1)
# invfreq2b, rel_time2b= compute_hwhm(dist2[1:],avg_corr_as[1,1:]-avg_ar_ao_as[1,1:], 1)

invfreq1c, rel_time1c= compute_hwhm(dist1[1:],avg_corr_s[2,1:]-avg_ar_ao_s[2,1:], 1)
# invfreq2c, rel_time2c= compute_hwhm(dist2[1:],avg_corr_as[2,1:]-avg_ar_ao_as[2,1:], 1)


# invfreq1a, rel_time1a, std_invfreq1a, std_rel_time1a= compute_hwhm(dist1[1:],avg_corr_s[0,1:]-avg_ar_ao_s[0,1:], 1, std_avg_corr_array_s[0,1:])
# # invfreq2a, rel_time2a= compute_hwhm(dist2[1:],avg_corr_as[0,1:]-avg_ar_ao_as[0,1:], 1)

# invfreq1b, rel_time1b, std_invfreq1b, std_rel_time1b= compute_hwhm(dist1[1:],avg_corr_s[1,1:]-avg_ar_ao_s[1,1:], 1, std_avg_corr_array_s[1,1:])
# # invfreq2b, rel_time2b= compute_hwhm(dist2[1:],avg_corr_as[1,1:]-avg_ar_ao_as[1,1:], 1)

# invfreq1c, rel_time1c, std_invfreq1c, std_rel_time1c= compute_hwhm(dist1[1:],avg_corr_s[2,1:]-avg_ar_ao_s[2,1:], 1, std_avg_corr_array_s[2,1:])
# # invfreq2c, rel_time2c= compute_hwhm(dist2[1:],avg_corr_as[2,1:]-avg_ar_ao_as[2,1:], 1)

###########################################################################################################################################

std_invfreq1a, std_rel_time1a= compute_hwhm(dist1[1:],avg_corr_s[0,1:]-avg_ar_ao_s[0,1:]+std_avg_corr_array_s[0,1:], 1)
# # std_invfreq2a, std_rel_time2a= compute_hwhm(dist2[1:],avg_corr_as[0,1:]-avg_ar_ao_as[0,1:]+std_avg_corr_array_as[0,1:], 1)

std_invfreq1b, std_rel_time1b= compute_hwhm(dist1[1:],avg_corr_s[1,1:]-avg_ar_ao_s[1,1:]+std_avg_corr_array_s[1,1:], 1)
# # std_invfreq2b, std_rel_time2b= compute_hwhm(dist2[1:],avg_corr_as[1,1:]-avg_ar_ao_as[1,1:]+std_avg_corr_array_as[1,1:], 1)

std_invfreq1c, std_rel_time1c= compute_hwhm(dist1[1:],avg_corr_s[2,1:]-avg_ar_ao_s[2,1:]+std_avg_corr_array_s[2,1:], 1)
# # std_invfreq2c, std_rel_time2c= compute_hwhm(dist2[1:],avg_corr_as[2,1:]-avg_ar_ao_as[2,1:]+std_avg_corr_array_as[2,1:], 1)


std_invfreq1a, std_rel_time1a= np.abs(std_invfreq1a-invfreq1a), np.abs(std_rel_time1a-rel_time1a) 
# # std_invfreq2a, std_rel_time2a= np.abs(std_invfreq2a-invfreq2a), np.abs(std_rel_time2a-rel_time2a)

std_invfreq1b, std_rel_time1b= np.abs(std_invfreq1b-invfreq1b), np.abs(std_rel_time1b-rel_time1b)
# # std_invfreq2b, std_rel_time2b= np.abs(std_invfreq2b-invfreq2b), np.abs(std_rel_time2b-rel_time2b)

std_invfreq1c, std_rel_time1c= np.abs(std_invfreq1c-invfreq1c), np.abs(std_rel_time1c-rel_time1c)
# # std_invfreq2c, std_rel_time2c= np.abs(std_invfreq2c-invfreq2c), np.abs(std_rel_time2c-rel_time2c)










##############################################################################################################################################



# xf_s=np.linspace(0.0, 1.0/(2.0*dL), L_sym/2)
# xf_a=np.linspace(0.0, 1.0/(2.0*dL), L_asy/2)
# # Fourier plots of final time ET
# caaS_f = scipy.fftpack.fft(avg_corr_s1)
# caaA_f = scipy.fftpack.fft(avg_corr_as1)

# cbbS_f = scipy.fftpack.fft(avg_corr_s2)
# cbbA_f = scipy.fftpack.fft(avg_corr_as2)

# cccS_f = scipy.fftpack.fft(avg_corr_s3)
# cccA_f = scipy.fftpack.fft(avg_corr_as3)


# cabS_f = scipy.fftpack.fft(corr_func_sym[4,:,tot_time-1])
# cabA_f = scipy.fftpack.fft(corr_func_asym[4,:,tot_time-1])

#ca0_f = scipy.fftpack.fft(corr_func[3,:,tot_time-1])
#


# corr_func_sym[0,:,:]-corr_term_sym[0,:,:]
# corr_func_asym[0,:,:]-corr_term_asym[0,:,:]
# ax2.plot(dist1,corr_func_sym[0,:,tot_time-1]-corr_term_sym[0,:,tot_time-1], dist1,corr_func_sym[1,:,tot_time-1]-corr_term_sym[1,:,tot_time-1], dist1,corr_func_sym[2,:,tot_time-1]-corr_term_sym[2,:,tot_time-1], dist2,corr_func_asym[0,:,tot_time-1]-corr_term_asym[0,:,tot_time-1], dist2,corr_func_asym[1,:,tot_time-1]-corr_term_asym[1,:,tot_time-1], dist2,corr_func_asym[2,:,tot_time-1]-corr_term_asym[2,:,tot_time-1])

# ax2.plot(dist1,corr_term_sym[0,:,tot_time-1], dist1,corr_term_sym[1,:,tot_time-1], dist1,corr_term_sym[2,:,tot_time-1], dist2,corr_term_asym[0,:,tot_time-1], dist2,corr_term_asym[1,:,tot_time-1], dist2,corr_term_asym[2,:,tot_time-1])


# ax2.plot(dist1,corr_func_sym[0,:,tot_time-1]-corr_term_sym[0,:,tot_time-1], dist1,corr_func_sym[1,:,tot_time-1]-corr_term_sym[1,:,tot_time-1], dist1,corr_func_sym[2,:,tot_time-1]-corr_term_sym[2,:,tot_time-1], dist2,corr_func_asym[0,:,tot_time-1]-corr_term_asym[0,:,tot_time-1], dist2,corr_func_asym[1,:,tot_time-1]-corr_term_asym[1,:,tot_time-1], dist2,corr_func_asym[2,:,tot_time-1]-corr_term_asym[2,:,tot_time-1])




# ax2.plot(xf_s[1:], 2.0/L_sym * np.abs(caaS_f[1:L_sym//2]), xf_s[1:], 2.0/L_sym * np.abs(cbbS_f[1:L_sym//2]), xf_s[1:], 2.0/L_sym * np.abs(cccS_f[1:L_sym//2]), xf_a[1:], 2.0/L_asy * np.abs(caaA_f[1:L_asy//2]), xf_a[1:], 2.0/L_asy * np.abs(cbbA_f[1:L_asy//2]), xf_a[1:], 2.0/L_asy * np.abs(cccA_f[1:L_asy//2]) )


# # plt.plot(xf_a[1:], 2.0/L_asy * np.abs(caa2_f[1:L_asy//2]))
# #plt.plot(xf[1:], 2.0/L * np.abs(ca0_f[1:L//2]))
# #
# #plt.legend(('c_aa', 'c_ab', 'c_a0'), loc='upper right')
# ax2.set_xlabel('wave modes')
# ax2.set_ylabel('EqT corr func')
# ax2.grid(True)
# ax2.legend((r'a_s', r'b_s', r'c_s', r'a_as',  r'a_as', r'a_as'), loc='upper right')
# ax2.tick_params(axis='both', which='major', labelsize=10)
# # ax2.set_title('Correlation plot, S/A=1:1 D = %s' % diff )

# # plt.legend(('c_aa1', 'c_aa2'), loc='upper right')
# # plt.xlabel('wave modes')
# # plt.ylabel('ET corr func')
# #plt.xlim(0.0, 0.35)
# #plt.ylim(0.0, 0.02)


# fig2.savefig(printdir+set_filename+'/'+folder_name+'/'+'EQ_corr_plot'+"ss"+sitesize+"nloc"+number_loc+"spw"+spawntime+'.png')

# invfreq1a, rel_time1a= compute_hwhm(xf_s[1:],2.0/L_sym * np.abs(caaS_f[1:L_sym//2]), dL)
# invfreq2a, rel_time2a= compute_hwhm(xf_a[1:],2.0/L_asy * np.abs(caaA_f[1:L_asy//2]), dL)

# invfreq1b, rel_time1b= compute_hwhm(xf_s[1:],2.0/L_sym * np.abs(cbbS_f[1:L_sym//2]), dL)
# invfreq2b, rel_time2b= compute_hwhm(xf_a[1:],2.0/L_asy * np.abs(cbbA_f[1:L_asy//2]), dL)

# invfreq1c, rel_time1c= compute_hwhm(xf_s[1:],2.0/L_sym * np.abs(cccS_f[1:L_sym//2]), dL)
# invfreq2c, rel_time2c= compute_hwhm(xf_a[1:],2.0/L_asy * np.abs(cccA_f[1:L_asy//2]), dL)





# invfreq3, rel_time3= compute_hwhm(xf_s[1:],2.0/L_sym * np.abs(cabS_f[1:L_sym//2]), dL)
# invfreq4, rel_time4= compute_hwhm(xf_a[1:],2.0/L_asy * np.abs(cabA_f[1:L_asy//2]), dL)


with open('Sdatafile.txt','a+') as f:
 	#f.write('roots are : %s %s '  % (str(r1), str(r2)) )
	f.write('inv_peak wave number, corr length in  symm region is : %s %s \n'  % (str(invfreq1a), str(rel_time1a)) )
	# f.write('inv_peak wave number, corr length in asymm region is : %s %s \n'  % (str(invfreq2a), str(rel_time2a)) )
	f.write('inv_peak wave number, corr length in symm region is : %s %s \n'  % (str(invfreq1b), str(rel_time1b)) )
	# f.write('inv_peak wave number, corr length in asymm region is : %s %s \n'  % (str(invfreq2b), str(rel_time2b)) ) 	
	f.write('inv_peak wave number, corr length in symm region is : %s %s \n'  % (str(invfreq1c), str(rel_time1c	)) )
	# f.write('inv_peak wave number, corr length in asymm region is : %s %s \n'  % (str(invfreq2c), str(rel_time2c)) )


# np.savetxt("d_array"+"ss"+sitesize+"nloc"+number_loc+"spw"+spawntime+".csv", d_array, delimiter=",", fmt="%1.7f")
# with file("corr_func"+"ss"+sitesize+"nloc"+number_loc+"spw"+spawntime+".csv",'w') as outfile:
# 	for slice_2d in corr_func:
# 		np.savetxt(outfile, slice_2d, delimiter=",", fmt= "%s")

###################  changes for correlation length #############################
# take select time region for when we want to analyse the data 
# 

###################  changes for correlation length #############################

RTCL_array=[float(asy_fact),rel_time1sy,rel_time2sy,rel_time3sy,rel_time1a,rel_time1b,rel_time1c]
OSC_array=[float(asy_fact),invfreq1sy,invfreq2sy,invfreq3sy]

with open(printdir+set_filename+'/'+'lengthx.csv', 'a') as csvfile:
    wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)#,quotechar='|', quoting=csv.QUOTE_MINIMAL
    wr.writerow(np.asarray(RTCL_array))
    
with open(printdir+set_filename+'/'+'timex.csv', 'a') as csvfile:
    wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)#,quotechar='|', quoting=csv.QUOTE_MINIMAL
    wr.writerow(np.asarray(OSC_array))

std_RTCL_array=[float(asy_fact),std_rel_time1sy,std_rel_time2sy,std_rel_time3sy,std_rel_time1a,std_rel_time1b,std_rel_time1c]
std_OSC_array=[float(asy_fact),std_invfreq1sy,std_invfreq2sy,std_invfreq3sy]

with open(printdir+set_filename+'/'+'std_lengthx.csv', 'a') as csvfile:
    wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)#,quotechar='|', quoting=csv.QUOTE_MINIMAL
    wr.writerow(np.asarray(std_RTCL_array))
    
with open(printdir+set_filename+'/'+'std_timex.csv', 'a') as csvfile:
    wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)#,quotechar='|', quoting=csv.QUOTE_MINIMAL
    wr.writerow(np.asarray(std_OSC_array))











# fig4, (ax1, ax2) = plt.subplots(2, 1)
# fig4.subplots_adjust(hspace=0.5)



# L_sym=int(float(Nx))
# L_asy=int(float(R*Nx))
# dL=1

# dist1=np.arange(L_sym)
# dist2=np.arange(L_asy)

# # ax1.plot(dist1,sym_ar[0,:,tot_time-1], dist1,sym_ar[1,:,tot_time-1], dist1,sym_ar[2,:,tot_time-1], dist2,asym_ar[1,:,tot_time-1], dist2,asym_ar[1,:,tot_time-1], dist2,asym_ar[2,:,tot_time-1])
# ax1.plot(dist1,sym_ar[0,:,tot_time-1], dist2,asym_ar[1,:,tot_time-1])
# # Log version
# # ax1.plot(dist1,np.log(avg_corr_s1), dist1,np.log(avg_corr_s2), dist1,np.log(avg_corr_s3), dist2,np.log(avg_corr_as1), dist2,np.log(avg_corr_as2), dist2,np.log(avg_corr_as3))
# #

# ax1.set_xlabel(r'r')
# ax1.set_ylabel(r'C(r)')
# ax1.grid(True)
# # ax1.legend((r'a_s', 'b_s', 'c_s', 'a_as',  'a_as', 'a_as'), loc='upper right')
# ax1.tick_params(axis='both', which='major', labelsize=10)
# ax1.set_title('Correlation plot, S/A=1:%d D = %s' % (R, diff))


# ax2.plot(dist1,corr_term_sym[0,:,tot_time-1], dist1,corr_term_sym[1,:,tot_time-1], dist1,corr_term_sym[2,:,tot_time-1], dist2,corr_term_asym[0,:,tot_time-1], dist2,corr_term_asym[1,:,tot_time-1], dist2,corr_term_asym[2,:,tot_time-1])



# ax2.set_xlabel('wave modes')
# ax2.set_ylabel('EqT corr func')
# ax2.grid(True)
# ax2.legend((r'a_s', r'b_s', r'c_s', r'a_as',  r'a_as', r'a_as'), loc='upper right')
# ax2.tick_params(axis='both', which='major', labelsize=10)
# # ax2.set_title('Correlation plot, S/A=1:1 D = %s' % diff )

# fig4.savefig(printdir+set_filename+'/'+folder_name+'/'+'test.png')



# ##############################################  Time corr changes #################################################################

# fig3, (ax1, ax2) = plt.subplots(2, 1)
# fig3.subplots_adjust(hspace=0.5)

# ac_sym_0=correlate(corr_func_sym[0,0,:],corr_func_sym[0,0,:])
# ac_sym_4=correlate(corr_func_sym[4,0,:],corr_func_sym[4,0,:])
# ac_sym_7=correlate(corr_func_sym[7,0,:],corr_func_sym[7,0,:])

# ax1.plot(T*time,ac_sym_0[tot_time-1:],T*time,ac_sym_4[tot_time-1:],T*time,ac_sym_7[tot_time-1:])
# ax1.set_xlabel('time in MC steps')
# ax1.set_ylabel('autocorrelation')
# ax1.grid(True)
# ax1.set_title('Density_plot, S/A=1:%d, D= %s ' % (R, diff) )
# ax1.legend(('sym_a', 'sym_b', 'sym_c'), loc='upper right')
# ax1.tick_params(axis='both', which='major', labelsize=10)

# tf=np.linspace(0.0, 1.0/(2.0*T), np.ceil(N/2+1))
# ds_1_f = scipy.fftpack.fft(ac_sym_0[tot_time-1:])
# ds_2_f = scipy.fftpack.fft(ac_sym_4[tot_time-1:])
# ds_3_f = scipy.fftpack.fft(ac_sym_7[tot_time-1:])

# ax2.plot(tf[1:], 2.0/N * np.abs(ds_1_f[1:N//2+1]), tf[1:], 2.0/N * np.abs(ds_2_f[1:N//2+1]), tf[1:], 2.0/N * np.abs(ds_3_f[1:N//2+1]))
# ax2.set_xlabel('frequency')
# ax2.set_ylabel('autocorrelation')
# ax2.grid(True)
# # ax2.set_title('Density_plot, S/A=1:1, D= %s ' % diff )
# ax2.legend(('sym_a', 'sym_b', 'sym_c'), loc='upper right')

# ax2.tick_params(axis='both', which='major', labelsize=10)

# fig3.savefig(printdir+set_filename+'/'+folder_name+'/'+'acorr_plot'+"ss"+sitesize+"nloc"+number_loc+"spw"+spawntime+'.png')

# invfreq1s, rel_time1s= compute_fwhm(tf[1:],2.0/N * np.abs(ds_1_f[1:N//2+1]), T)
# invfreq2s, rel_time2s= compute_fwhm(tf[1:],2.0/N * np.abs(ds_2_f[1:N//2+1]), T)
# invfreq3s, rel_time3s= compute_fwhm(tf[1:],2.0/N * np.abs(ds_3_f[1:N//2+1]), T)


# with open('Sdatafile.txt','a+') as f:
# 	f.write('inv_resonance freq relaxation time is : %s %s \n'  % (str(invfreq1s), str(rel_time1s)) )
#  	f.write('inv_resonance freq relaxation time is : %s %s \n'  % (str(invfreq2s), str(rel_time2s)) )
#  	f.write('inv_resonance freq relaxation time is : %s %s \n'  % (str(invfreq3s), str(rel_time3s)) )

####################################################################################################################################