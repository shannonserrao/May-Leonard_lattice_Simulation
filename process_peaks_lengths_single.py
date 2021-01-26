import sys
import os
import errno
import numpy as np
from numpy import genfromtxt
import csv
import glob
import fnmatch
# from StringIO import StringIO
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.fftpack
from scipy.fftpack import fft
from scipy.interpolate import UnivariateSpline
from scipy.signal import correlate 
import itertools
from scipy.interpolate import spline

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

from shutil import copy2
#printdir = "/home/shann87/Control/main_prog/mayleonard_only/asy_rates_PMC_May30_decoupled_Jun21/"
printdir = "/home/shann87/Control/main_prog/mayleonard_only/asy_rates_PMC_May30_decoupled_SEP19/"
# ######################################################################################################
# def compute_fwhm(x,y, interval):
# 	""" Computes the full width half max and res freq 
# 	"""
# 	y_max=np.max(y)
# 	fmax=x[np.argmax(y)]
# 	spline = UnivariateSpline(x, y-y_max/2, s=0)
# 	#r1, r2 = spline.roots()
# 	r =spline.roots()
# 	rel_time=0.
# 	#print(r)
# 	if len(r)==1:
# 		print("No roots at full width half max, check data or graph.. ")
# 		fmax, rel_time = compute_hwhm(x,y, interval)
	
# 	# case when data has multiple roots or two roots
# 	elif len(r) > 1:
# 		# r_index=[]

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
# 	# else:
# 	# 	inv_fmax, rel_time = compute_hwhm(x,y, interval)	
# 	return fmax, rel_time

# def compute_hwhm(x,y, interval):
# 	""" Computes the half width half max and res freq 
# 	"""
	
# 	y_min=np.min(y)

# 	y_adj=y-y_min
# 	y_max=np.max(y_adj)
# 	fmax=x[np.argmax(y_adj)]
	

# 	spline = UnivariateSpline(x, y_adj-y_max/2, s=0)
# 	r1 = spline.roots()
	
# 	if len(r1) >1:
# 		print('returning all roots ' + str(r1) )
# 	if len(r1) >0:
# 		rel_time=r1[0]
# 	# else:
# 	# 	y_min=np.min(y)
# 	# 	spline = UnivariateSpline(x, y-y_min, s=0)
# 	# 	r1=spline.roots()
# 	# 	rel_time=r1[0]			
	
# 	# if r1 <0 or r2 <0:
# 	# 	print("returning half width half max")
	

# 	return fmax, rel_time
# #######################################################################################################



#df = pd.DataFrame(columns=['Name','age'])
#runs=len(list)
# boolean__init_shape=True

# #boxsize
# #boxsize=256
# boxsize=sys.argv[1]

# #parameter rates
# pred=sys.argv[2]
# rep=sys.argv[3]
# diff=sys.argv[4]

# #Define the corr array size =(species * length *total_time)
# printtime=sys.argv[5]
# endtime=sys.argv[6]
# #25000=sys.argv[6]
# #rootdir="/home/shann87/mayleonardmc/mayleonard_only/new_M28_corrl_ML_spirals256/"
# rootdir=sys.argv[7]

# #ab=pd.DataFrame(columns=['kappa','Vol','MTE'])
# #df = pd.DataFrame(columns=['Name','age'])

# spawntime=sys.argv[8]
# sitesize=sys.argv[9]
# number_loc=sys.argv[10]
# interval=sys.argv[12]
# #print(interval)
# printlattice_every=sys.argv[11]
# #print(interval)
# folder_name=sys.argv[13]

# Nx=int(sys.argv[14])
# R=int(sys.argv[15])
# asy_fact=float(sys.argv[16])
# set_filename=sys.argv[17]
# printdir=sys.argv[18]

# only take in all the files : print directory

def sort_data(x,y):
	# lists = sorted(itertools.izip(*[x, y]))
	# new_x, new_y = list(itertools.izip(*lists))
	xnew, ynew = zip(*sorted(zip(x, y)))
	# xnew = np.linspace(np.min(x), np.max(x), 14)  
	# ynew = spline(x, y, xnew)
	return xnew,ynew

def trois_couleurs_plt(k_values,values_a,values_b,values_c, print_string_x, print_string_y, R, diff_value,title,filename,yerr1,yerr2,yerr3):
	# k_values, corr_values= sort_data(corr_array[:,0], corr_array[:,10])
	# k_values_new = np.linspace(k_values.min(), k_values.max(), 300)
	# corr_values_new = spline(k_values, corr_values, k_values_new)
	if filename=='corr_asym':
		legend_loc='upper right'
	else:
		legend_loc='lower right'
	fig1, ax1 = plt.subplots()
	# print(len(yerr1))
	if (np.all(np.ceil(values_a)==0)==False):
		print(values_a)
		# k_values=np.asarray(k_values)
		# values_a=np.asarray(values_a)
		# values_b=np.asarray(values_b)
		# values_c=np.asarray(values_c)
		values_i=np.ceil(np.asarray(values_a)).astype(int)
		# print(values_i>0)
		# print(values_i)
		# print(k_values)
		# print(values_a)
		# print(yerr1)
		ax1.errorbar(np.array(k_values)[values_i.astype(int)>0],np.array(values_a)[values_i.astype(int)>0],yerr=np.array(yerr1)[values_i.astype(int)>0],fmt='ro',elinewidth=10.0)
		ax1.errorbar(np.array(k_values)[values_i.astype(int)>0],np.array(values_b)[values_i.astype(int)>0],yerr=np.array(yerr2)[values_i.astype(int)>0],fmt='go',elinewidth=10.0)
		ax1.errorbar(np.array(k_values)[values_i.astype(int)>0],np.array(values_c)[values_i.astype(int)>0],yerr=np.array(yerr3)[values_i.astype(int)>0],fmt='bo',elinewidth=10.0)
		# ax1.errorbar(k_values[values_i>0],values_a[values_i>0],yerr=yerr1[values_i>0],fmt='ro',elinewidth=10.0)
		# ax1.errorbar(k_values[values_i>0],values_b[values_i>0],yerr=yerr2[values_i>0],fmt='go',elinewidth=10.0)
		# ax1.errorbar(k_values[values_i>0],values_c[values_i>0],yerr=yerr3[values_i>0],fmt='bo',elinewidth=10.0)
		# ax1.errorbar(k_values,values_a,yerr=yerr1,fmt='ro',elinewidth=10.0)
		# ax1.errorbar(k_values,values_b,yerr=yerr2,fmt='go',elinewidth=10.0)
		# ax1.errorbar(k_values,values_c,yerr=yerr3,fmt='bo',elinewidth=10.0)
		# ax1.scatter(k_values,values_a,c='r')
		# ax1.scatter(k_values,values_b,c='g')
		# ax1.scatter(k_values,values_c,c='b')
		# ax1.plot(k_values,np.exp(-3*k_values+2)+5)
		ax1.set_xlabel(print_string_x, fontsize=25)
		ax1.set_ylabel(print_string_y, fontsize=25)
		ax1.grid(True)
		ax1.set_title(title+' D = %s ' % (diff_value) )
		ax1.legend((r'$a$', r'$b$', r'$c$'), loc=legend_loc, fontsize=20)
		ax1.tick_params(axis='both', which='major', labelsize=18)
		ax1.xaxis.set_major_locator(plt.MaxNLocator(6))
		fig1.savefig(printdir+filename+'_vs_k'+str(diff_value).replace('.','p')+'.png',bbox_inches='tight')

# def plot_criteria(k_values,xi_a,xi_b,xi_c, om_a, om_b,om_c, dif, mu=0.2, sig=0.2):
# 	''' Plot the graphs that illuminate the criteria of the spiral creation
# 	'''
# 	legend_loc='upper right'
# 	fig2, ax2 = plt.subplots()
# 	xi_a=np.asarray(xi_a)
# 	xi_b=np.asarray(xi_b)
# 	xi_c=np.asarray(xi_c)
# 	om_a=np.asarray(om_a)
# 	om_b=np.asarray(om_b)
# 	om_c=np.asarray(om_c)

# 	xi_mean=(xi_a+xi_b+xi_c)/3
# 	xi_pred=np.sqrt(dif* np.cbrt(k_values)/ mu)
# 	ax2.plot(k_values,xi_mean, 'go', k_values, 10*xi_pred, 'b^')
# 	ax2.legend((r'$\xi_m$', r'$10*\xi_p$'), loc=legend_loc, fontsize=10)
# 	ax2.tick_params(axis='both', which='major', labelsize=10)
# 	fig2.savefig(printdir+'xi_comparison'+str(diff_value).replace('.','p')+'.png')

# 	fig3, ax3 = plt.subplots()
# 	om_mean=(om_a+om_b+om_c)/3
# 	om_pred=np.power(xi_pred, -2)*dif
# 	om_xicomp=np.power(xi_mean, -2)*dif
# 	ax3.plot(k_values,100*om_mean, 'go', k_values, om_pred, 'b^', k_values, 100*om_xicomp, 'r*')
# 	ax3.legend((r'$100*\omega_m$', r'$\omega_p$', r'100*$\omega_\xi$' ), loc=legend_loc, fontsize=10)
# 	ax3.tick_params(axis='both', which='major', labelsize=10)
# 	fig3.savefig(printdir+'om_comparison'+str(diff_value).replace('.','p')+'.png')
def plot_criteria(k_values,xi_a,xi_b,xi_c, om_a, om_b,om_c, dif, mu=0.2, sig=0.2):
	''' Plot the graphs that illuminate the criteria of the spiral creation
	'''
	legend_loc='lower right'
	fig2, ax2 = plt.subplots()
	k_array=np.asarray(k_values)
	xi_a=np.asarray(xi_a)
	xi_b=np.asarray(xi_b)
	xi_c=np.asarray(xi_c)
	om_a=np.asarray(om_a)
	om_b=np.asarray(om_b)
	om_c=np.asarray(om_c)

	xi_ind=np.zeros(k_array.shape[0])
	# for i,k in enumerate(k_array):
	# 	if k <= 1:
	# 		xi_ind[i]=1*np.min((xi_a[i],xi_b[i],xi_c[i]))
	# 	else:
	# 		Arr=np.array([xi_a[i],xi_b[i],xi_c[i]])
	# 		ind=np.argpartition(Arr,2)
	# 		xi_ind[i]=np.sum(Arr[ind[:2]])
	# 		# xi_ind=2*np.min((xi_a[i],xi_b[i],xi_c[i]))

	# xi_mean=(xi_a+xi_b+xi_c)/3
	# xi_pred=np.sqrt(dif* np.cbrt(k_values)/ mu)
	# ax2.plot(k_array,xi_mean, 'go', k_values, 10*xi_pred, 'b^')
	# ax2.legend((r'$\xi_m$', r'$10*\xi_p$'), loc=legend_loc, fontsize=10)
	# ax2.tick_params(axis='both', which='major', labelsize=10)
	# fig2.savefig(printdir+'xi_comparison'+str(diff_value).replace('.','p')+'.png')
	xi_ind=np.power(np.power(xi_a,1)+np.power(xi_b,1)+np.power(xi_c,1),1)
	fig3, ax3 = plt.subplots()
	om_mean=(om_a+om_b+om_c)/3
	# om_pred=np.power(xi_pred, -2)*dif
	om_xicomp=np.power(xi_ind/3, -2)*dif*(2.0/3)
	ax3.plot(k_values,om_mean, 'go', k_values, om_xicomp, 'r*')
	ax3.legend((r'$\omega_m$',  r'$\omega_c$' ), loc=legend_loc, fontsize=25)
	ax3.tick_params(axis='both', which='major', labelsize=18)
	ax3.set_xlabel(r' k', fontsize=25)
	ax3.set_ylabel(r'$\omega $', fontsize=25)
	ax3.xaxis.set_major_locator(plt.MaxNLocator(6))
	fig3.savefig(printdir+'om_criteria_new'+str(diff_value).replace('.','p')+'.png',bbox_inches='tight')



diff_dict={'2': '0.8','5': '0.1', '6': '5.0','7': '3.0', '8': '4.0', '10': '0.5'}
R=3

os.chdir(printdir)
list=os.listdir(printdir+'/')
diff_values=len(list)

for set_name in list:
	print(set_name)
	set_num=set_name[8:]
	print(float(set_name[8:]))
	diff_value=float(diff_dict[set_num])
	os.chdir(printdir+set_name+'/')
	print(os.getcwd())
	time_array=np.char.strip(np.genfromtxt("timex.csv", dtype="|U30",delimiter=","), chars='"')
	time_array=	time_array.astype(float)
	# print(time_array)
	corr_array=np.char.strip(np.genfromtxt("lengthx.csv", dtype="|U30",delimiter=","), chars='"')
	corr_array=	corr_array.astype(float)
	# print(corr_array)

	std_time_array=np.char.strip(np.genfromtxt("std_timex.csv", dtype="|U30",delimiter=","), chars='"')
	std_time_array=	std_time_array.astype(float)
	# print(time_array)
	std_corr_array=np.char.strip(np.genfromtxt("std_lengthx.csv", dtype="|U30",delimiter=","), chars='"')
	std_corr_array=	std_corr_array.astype(float)
	
	
	k_values, corr_values_a= sort_data(corr_array[:,0], corr_array[:,4])
	k_values, corr_values_b= sort_data(corr_array[:,0], corr_array[:,5])
	k_values, corr_values_c= sort_data(corr_array[:,0], corr_array[:,6])

	k_values, std_corr_values_a= sort_data(std_corr_array[:,0], std_corr_array[:,4])
	k_values, std_corr_values_b= sort_data(std_corr_array[:,0], std_corr_array[:,5])
	k_values, std_corr_values_c= sort_data(std_corr_array[:,0], std_corr_array[:,6])
	
	trois_couleurs_plt(k_values,corr_values_a,corr_values_b,corr_values_c,r'$k \ $',r'$ \xi\  $', R, diff_value, ' ','corr',std_corr_values_a,std_corr_values_b,std_corr_values_c)
	
	# k_values, osc_peak_values_a= sort_data(time_array[:,0], time_array[:,4])
	# k_values, osc_peak_values_b= sort_data(time_array[:,0], time_array[:,5])
	# k_values, osc_peak_values_c= sort_data(time_array[:,0], time_array[:,6])

	# k_values, std_osc_peak_values_a= sort_data(std_time_array[:,0], std_time_array[:,4])
	# k_values, std_osc_peak_values_b= sort_data(std_time_array[:,0], std_time_array[:,5])
	# k_values, std_osc_peak_values_c= sort_data(std_time_array[:,0], std_time_array[:,6])

	# trois_couleurs_plt(k_values,osc_peak_values_a,osc_peak_values_b,osc_peak_values_c,r'$k\ asymmetry\ factor\ $',r'$oscillation\ frequency$', R, diff_value,'Oscillation freq in the asymmetric region','osc_asym',std_osc_peak_values_a,std_osc_peak_values_b,std_osc_peak_values_c)

	k_values, osc_peak_values_a= sort_data(time_array[:,0], time_array[:,1])
	k_values, osc_peak_values_b= sort_data(time_array[:,0], time_array[:,2])
	k_values, osc_peak_values_c= sort_data(time_array[:,0], time_array[:,3])

	k_values, std_osc_peak_values_a= sort_data(std_time_array[:,0], std_time_array[:,1])
	k_values, std_osc_peak_values_b= sort_data(std_time_array[:,0], std_time_array[:,2])
	k_values, std_osc_peak_values_c= sort_data(std_time_array[:,0], std_time_array[:,3])
	trois_couleurs_plt(k_values,osc_peak_values_a,osc_peak_values_b,osc_peak_values_c,r'$k\  $',r'$\omega\ $', R, diff_value,' ','osc',std_osc_peak_values_a,std_osc_peak_values_b,std_osc_peak_values_c)


	# k_values, rel_time_values_a= sort_data(corr_array[:,0], corr_array[:,4])
	# k_values, rel_time_values_b= sort_data(corr_array[:,0], corr_array[:,5])
	# k_values, rel_time_values_c= sort_data(corr_array[:,0], corr_array[:,6])

	# k_values, std_rel_time_values_a= sort_data(std_corr_array[:,0], std_corr_array[:,4])
	# k_values, std_rel_time_values_b= sort_data(std_corr_array[:,0], std_corr_array[:,5])
	# k_values, std_rel_time_values_c= sort_data(std_corr_array[:,0], std_corr_array[:,6])
	# trois_couleurs_plt(k_values,rel_time_values_a,rel_time_values_b,rel_time_values_c,'k','rel_times', R, diff_value,'relaxation time in asymmetric region','rel_time_asym',std_rel_time_values_a,std_rel_time_values_b,std_rel_time_values_c)	
# # 	for name in files:
	k_values, rel_time_values_a_s= sort_data(corr_array[:,0], corr_array[:,1])
	k_values, rel_time_values_b_s= sort_data(corr_array[:,0], corr_array[:,2])
	k_values, rel_time_values_c_s= sort_data(corr_array[:,0], corr_array[:,3])

	k_values, std_rel_time_values_a= sort_data(std_corr_array[:,0], std_corr_array[:,1])
	k_values, std_rel_time_values_b= sort_data(std_corr_array[:,0], std_corr_array[:,2])
	k_values, std_rel_time_values_c= sort_data(std_corr_array[:,0], std_corr_array[:,3])
	trois_couleurs_plt(k_values,rel_time_values_a_s,rel_time_values_b_s,rel_time_values_c_s,'k','rel_times', R, diff_value,'relaxation time','rel_time_symm',std_rel_time_values_a,std_rel_time_values_b,std_rel_time_values_c)	

	k_values, a_osc_peak_values_a= sort_data(time_array[:,0], time_array[:,1])
	k_values, a_osc_peak_values_b= sort_data(time_array[:,0], time_array[:,2])
	k_values, a_osc_peak_values_c= sort_data(time_array[:,0], time_array[:,3])
	k_values, a_rel_time_values_a= sort_data(corr_array[:,0], corr_array[:,4])
	k_values, a_rel_time_values_b= sort_data(corr_array[:,0], corr_array[:,5])
	k_values, a_rel_time_values_c= sort_data(corr_array[:,0], corr_array[:,6])

	plot_criteria(k_values,a_rel_time_values_a,a_rel_time_values_b,a_rel_time_values_c, a_osc_peak_values_a, a_osc_peak_values_b,a_osc_peak_values_c, diff_value, 0.2, 0.2)
