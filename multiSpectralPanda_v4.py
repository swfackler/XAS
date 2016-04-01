#!/usr/bin/python

#
#   Plots MH loops and XAS. XMCD spectra from BL 6.3.1 This is not even an alpha Version (no pun intended...)
#   atndiaye@lbl.gov 
#
import pandas as pd
import matplotlib.pylab as plt
import matplotlib.transforms as mtransforms
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import scipy
from scipy.interpolate import griddata
import sys
import os.path
import numpy as np
import mspScan as msp
import mspAnalyze as mspA
import pickle as pickle
from mpl_toolkits.axes_grid1 import AxesGrid
import math
import itertools


def sort_zy_lines(a,b):
	if hasattr(a, '__iter__') and hasattr(b, '__iter__'):
		if (len(a)>1) & (len(b)>1) :
			if 	a[0]>b[0]: return  1
			elif	a[0]<b[0]: return -1
			elif	a[1]<b[1]: return  1
			elif	a[1]>b[1]: return -1
			else: return 0
	# default:
	if 	a>b: return  1
	elif	a<b: return -1
	else:	return 0


def treat_mh(scan_groups):
	#Start taking the data from the groups
		keys = scan_groups.groups.keys()
		keys.sort()
		elements = [ 	
#				( 635.0, 638.5, "Manganese"),
#				( 721.1, 708.0, "Iron"),
#				( 777.0, 793.2, "Cobalt"),
				( 780.0, 795.0, "Cobalt"),
#				( 852.6, 871.15,"Nickel"),
				]
#		elements = [ 
#				( 779.7, 795.0, "Cobalt"),
#	#			( 853.4, 871.1,"Nickel"),
#				]

		keys = scan_groups.groups.keys()
		keys.sort(sort_zy_lines)
		mh_frame = pd.DataFrame()

		print 'MH Keys:',keys
		for g in keys:
			sg_atloc = scan_groups.get_group(g)
			mhlist = list()
			for E1, E2, element in elements:
				mh_data 	= mspA.get_mh_data(sg_atloc , E1, E2, element)
				mh_data["Index"]= mh_data.index.values

				for c,k in zip(group_columns,g):
					mh_data[c] = k
					print "c,k,",c,k
				mhlist.append(mh_data)
			temp_mh  = pd.concat(mhlist)

			mh_frame = pd.concat([mh_frame,temp_mh ])
			print "###########"
		return(mh_frame)

def average_spec(mh_frame, columns=[]):
#	print "Grouping by ",columns
	mf 		= mh_frame.groupby(columns).agg([np.mean, np.std])

	# Flattening Coloumns and Replacing 'Colounmname mean' by 'Coloumname'
	mf.columns 	= [' '.join([c.replace("mean","") for c in col]).strip() for col in mf.columns.values]

	# Flatten index and restoring order of datapoints
	mf 		= mf.reset_index()
	mf 		= mf.set_index('Index')
	mf 		= mf.sort()
	return mf

def average_mh(mh_frame, group_columns=[]):
	columns		= ["Magnet Field","branch","Element"]
	columns.extend(group_columns)
	return  average_spec(mh_frame, columns=columns)

def average_xmcd(mh_frame, group_columns=[]):
	columns		= ["Energy",]
	columns.extend(group_columns)
	return  average_spec(mh_frame, columns=columns)


def make_ticklabels_invisible(axs1,axs2,x_plotnumber, y_plotnumbner):
	for i, (ax1, ax2) in enumerate(zip(axs1, axs2)):
		x = int(i % y_plotnumber)
		y = int(i / y_plotnumber)
		#or (x ==  x_plotnumber) or (y==0) or (y==y_plotnumber)):
		for tl in ax1.get_xticklabels() + ax1.get_yticklabels()+ ax2.get_xticklabels() + ax2.get_yticklabels():
			tl.set_visible(False)	
# 		if not( (x==0) or (x ==  y_plotnumber-1) or (y==0) or (y==x_plotnumber-1))
		if (x==0): 		[tl.set_visible(True) for tl in ax1.get_yticklabels()]
		if (x==y_plotnumber-1):	[tl.set_visible(True) for tl in ax2.get_yticklabels()]
		if (y==0): 		[tl.set_visible(True) for tl in ax2.get_xticklabels()]
		if (y==x_plotnumber-1):	[tl.set_visible(True) for tl in ax1.get_xticklabels()]
		if (y==0): 		[tl.set_rotation('vertical') for tl in ax2.get_xticklabels()]
		if (y==x_plotnumber-1):	[tl.set_rotation('vertical') for tl in ax1.get_xticklabels()]
		ax1.get_xticklabels()[ 0].set_visible(False)
		ax1.get_xticklabels()[-1].set_visible(False)
		ax2.get_xticklabels()[ 0].set_visible(False)
		ax2.get_xticklabels()[-1].set_visible(False)

		ax2.tick_params(axis='y', colors='red')
		ax1.tick_params(axis='y', colors='blue')

def ax_setup(x_plotnumber, y_plotnumber):
			"""returns a list of axes for the subplots"""
		#	fig, axs = plt.subplots(x_plotnumber, y_plotnumber, figsize=(8,10.5), facecolor='w', edgecolor='k')
#			fig, axs = plt.subplts(x_plotnumber, y_plotnumber, figsize=(28,40), facecolor='w', edgecolor='k')
			fig = plt.figure(figsize=(12,12), facecolor='w', edgecolor='k')
#			fig.subplots_adjust(hspace = .2, wspace=.2)
			if x_plotnumber * y_plotnumber >=0:
				axs = list()
				for i in range( x_plotnumber * y_plotnumber):
					a =  SubplotHost(fig, 1, 1,1)
#					a =  SubplotHost(fig, x_plotnumber, y_plotnumber, i+1)
					axs.append(a)
			else:
				axs = [axs,]

			return fig, axs

		
def plot_mh_frame(mh_frame, group_columns, x_plotnumber, y_plotnumber ):

			fig, axs = ax_setup(x_plotnumber, y_plotnumber)

			mh_frame_group = mh_frame.groupby(group_columns)

			keys = mh_frame_group.groups.keys()
			print "Keys for plotting", keys
			keys.sort(sort_zy_lines)
			for g, ax1 in zip(keys, axs[::-1]):
				mh_data_atloc = mh_frame_group.get_group(g)
				mh_data_atloc = mh_data_atloc.groupby("Element")
				for element in mh_data_atloc.groups.keys():
					mh_data = mh_data_atloc.get_group(element)
					mspA.plot_mh(mh_data, ax1, label = element+" "+str([round(p,2) for p in g]), linewidth = 2)
								
				
#				box = ax1.get_position()
#				ax1.set_position([box.x0, box.y0+ box.height*.2, box.width , box.height*.8])
				# Put a legend to the right of the current axis
#				ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3))
				ax1.legend(loc='lower right')
#				ax1.set_ylim(-0.05,0.05)
#				ax1.set_ylim(-1.05,1.05)
#				ax1.set_title(g)
				fig.add_subplot(ax1)
				print "subplot added"

			fig.suptitle( "Overview "+str(filenames)+".pdf", fontsize=32)

#			fig.subplots_adjust(hspace=0, wspace=0)
#			plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
			plt.show()

def save_txt(mh_data, label):
#	filename = "MH TrajScan30869 " +label + ".txt"
	filename = "MH TrajScan30884 " +label + ".txt"
	mh_data.to_csv(filename)

def save_mh_frame(mh_frame, group_columns):
			mh_frame_group = mh_frame.groupby(group_columns)
			keys = mh_frame_group.groups.keys()
			print "Keys for plotting", keys
			keys.sort(sort_zy_lines)
			for g in keys:
				mh_data_atloc = mh_frame_group.get_group(g)
				mh_data_atloc = mh_data_atloc.groupby("Element")
				for element in mh_data_atloc.groups.keys():
					mh_data = mh_data_atloc.get_group(element)
					save_txt(mh_data, label = "MH_"+element+" "+str([round(p,2) for p in g]))
								


	

def treat_xmcd(scan_groups, element):
		keys = scan_groups.groups.keys()
		keys.sort(sort_zy_lines)
		xmcd_frame = pd.DataFrame()
		print "Keys for analysis", keys
		for g in keys:

			xmcd_data = mspA.get_xmcd(scan_groups.get_group(g)[(scan_groups.get_group(g)["Energy"]>element.E_lower) & (scan_groups.get_group(g)["Energy"]<element.E_upper)])
#			xmcd_data = mspA.removeBG(xmcd_data, pre_edge = element.pre_edge, post_edge= None)
			xmcd_data = mspA.removeV_BG(xmcd_data, pre_edge = element.pre_edge, post_edge= None)

			if not hasattr(g,'__iter__'):
				g =[g,]

			for c,k in zip(group_columns,g):
				xmcd_data[c] = k
			xmcd_data["Index"]= xmcd_data.index.values
			xmcd_frame = pd.concat([xmcd_frame,xmcd_data])

		return(xmcd_frame)


def plot_xmcd_frame(xmcd_frame, group_columns, x_plotnumber, y_plotnumber, mode="XMCD" ):

			keys = xmcd_frame_group.groups.keys()
			print "Keys for plotting", keys
			keys.sort(sort_zy_lines)
			
			fig, axs1 = ax_setup(x_plotnumber, y_plotnumber)

#			axs2  = [ax1.twinx() for ax1 in axs1]
			legendtype = "short"
			once = True

			if mode == "XAS":
				for g in keys:
					ax1,ax2   = mspA.plot_xmcd(xmcd_frame_group.get_group(g) ,ax1 = axs1[0] ,ax2=None, label = g, legend = legendtype)
					print "AX1", ax1
					if once:
						
						fig.add_subplot(ax1)
						once = False
#					ax1.set_ylim(-.65,1.3)		# normalized
#					ax1.set_ylim(0,0.1)		# 
					
				axs2 = axs1

			if (mode == "XAS") & False:
				for g, ax1 in zip(keys, itertools.cycle(axs1[::-1])):
					ax1,ax2   = mspA.plot_xmcd(xmcd_frame_group.get_group(g) ,ax1 = ax1 ,ax2=None, label = g, legend = legendtype)
					print "AX1", ax1
					fig.add_subplot(ax1)
#					ax1.set_ylim(-.65,1.3)		# normalized
#					ax1.set_ylim(0,0.1)		# 
					
				axs2 = axs1
			if mode == "XMCD":
#				aux_trans = mtransforms.Affine2D().scale(1.,2)
				aux_trans = mtransforms.Affine2D().scale(1.,1)
				axs2  = [ax1.twin(aux_trans) for ax1 in axs1]
				for ax2 in axs2: ax2.set_viewlim_mode("transform")
		
				for g, ax1, ax2 in zip(keys, itertools.cycle(axs1[::-1]), itertools.cycle(axs2[::-1])):
					ax1,ax2   = mspA.plot_xmcd(xmcd_frame_group.get_group(g) ,ax1,ax2, label = g ,legend = legendtype)
					fig.add_subplot(ax1)
#					ax1.set_ylim(-.65,1.3)		# normalized
					ax1.set_ylim(-.65,3)

		#			ax1.set_ylim(0.35,0.5)		# Fe (raw)
		#			ax2.set_ylim(-0.005,0.012)	# Fe (raw)
		#			ax1.set_ylim(0.35,0.5)		# Co (raw)
		#			ax2.set_ylim(-0.005,0.015)	# Co (raw)
		#			ax2.set_ylim(-0.06,0.04)	# Cu (raw)

			fig.subplots_adjust(hspace = 0, wspace=0)
			fig.suptitle( "Overview "+element.element+" "+str(filenames)+".pdf", fontsize=32)
			make_ticklabels_invisible(axs1,axs2,x_plotnumber,y_plotnumber)
			plt.show()
 


def get_plot_dimensions(scandata_f, group_columns):
			# Prepare plot and get axes
			# should maybe be based on the keys, rather than on the data
			# should also not be here!

			if len(group_columns) ==0:
				x_plotnumber		= 1
				y_plotnumber		= 1
			elif len(group_columns) == 1:
				x_plotnumber		= int(math.ceil(math.sqrt(len(scandata_f[group_columns[0]].unique()))))
				y_plotnumber		= int(math.ceil(math.sqrt(len(scandata_f[group_columns[0]].unique()))))
			else:
				x_plotnumber		= len(scandata_f[group_columns[0]].unique())
				y_plotnumber		= len(scandata_f[group_columns[1]].unique())
			
			print "Plot dimensions: ", x_plotnumber,y_plotnumber
			return x_plotnumber, y_plotnumber


def plot_map(x,y,z, map_fig, map_ax, title='', method = 'nearest', hide_points = False):
			#### Plot countour ####################################################################
			# define grid.
			xi = np.linspace(min(x),max(x),500)
			yi = np.linspace(min(y),max(y),500)

			# grid the data.
			zi = griddata((x, y), scipy.array(z), (xi[None,:], yi[:,None]), method=method)
			
			# contour the gridded data, plotting dots at the randomly spaced data points.
			CS = map_ax.contourf(xi,yi,zi,500,cmap=plt.cm.CMRmap )
#			CS = map_ax.contourf(xi,yi,zi,500,cmap=plt.cm.CMRmap, vmin=-0.3, vmax=1 )
#			CS = map_ax.contourf(xi,yi,zi,500,cmap=plt.cm.CMRmap, vmin=-0.05, vmax = 0.2)
#			CS = map_ax.contourf(xi,yi,zi,500,cmap=plt.cm.CMRmap, vmin=-0.02, vmax = 0.15)
			map_fig.colorbar(CS, ax=map_ax)
			# plot data points
			if not hide_points:
				map_ax.scatter(x,y,marker='o',c='b',s=5)
			map_ax.set_title(title)
			


class record:
	def __init__(self):
		self.__dict__ = dict()



if __name__ == "__main__":
#if True:
	# Read arguments
	mode 		= sys.argv[1]
	filenames	= sys.argv[2:]
else:
	filenames=['',]

if os.path.isfile(filenames[0]):
	
	# Loading Scan
	scandata_f = msp.read_scans(filenames, datacounter = "Counter 1")

# 	Group 1: 
#	all the groups which get an own subplot

### CHANGE GROUP COLOUMNS TO DEFINE WHICH DATA SHOULD BE GROUPED TOGETHER

#	sg = scandata_f.groupby(["filename_rump",]) 
#	group_columns = ["Z","Theta"]
	group_columns = ["Z", ]
#	group_columns = ["filename_rump"]

#	Z_conv 		= lambda z: z-70.9
#	Y_conv		= lambda y: (y+8.0)*2
#	offset = [Z_conv,Y_conv]
#	scandata_f["Z"] = Z_conv(scandata_f["Z"] )
#	scandata_f["Y"] = Y_conv(scandata_f["Y"] )


	sg = scandata_f.groupby(group_columns) # all the groups which get an own subplot

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .. . . . . . . . . . . . . . . . . . . . 

	if mode == 'MH':
		mh_frame = treat_mh(sg)

		plot_it = True
#		plot_it = False

		mh_frame = average_mh(mh_frame, group_columns)
		if plot_it:

			x_plotnumber, y_plotnumber	=	get_plot_dimensions(scandata_f, group_columns)

			#ARGH!
#			x_plotnumber = int(scipy.ceil(scipy.sqrt(x_plotnumber *y_plotnumber)))
#			y_plotnumber = int(scipy.ceil(scipy.sqrt(x_plotnumber *y_plotnumber)))

			plot_mh_frame(mh_frame, group_columns, x_plotnumber , y_plotnumber)

		save_it = True
		if save_it:
			save_mh_frame(mh_frame, group_columns)
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .. . . . . . . . . . . . . . . . . . . . 

#		make_array = True
		make_array = False
		if make_array:
			# NEEDS WORK : 
			# - generic support for multiple elements
			# - consolidate with the FORC version, which deals better (not great) with NANs

			fig = plt.figure(figsize=(12,12), facecolor='w', edgecolor='k')
	
#			ax1 = fig.add_subplot(121)
			ax2 = fig.add_subplot(111)
#			ax2 = fig.add_subplot(122)

#			aa = mh_frame[mh_frame["Element"] == "Nickel"]
#			ab = aa.apply(pd.Series.reset_index, drop=True).groupby('Z')
#			a  = ab.Asymmetry.apply(pd.Series.reset_index, drop=True).unstack()


			ba = mh_frame[mh_frame["Element"] == "Cobalt"]
			bb = ba.apply(pd.Series.reset_index, drop=True).groupby('Z')
			b  = bb.Asymmetry.apply(pd.Series.reset_index, drop=True).unstack()

#			ma = ax1.imshow(a.values, aspect='auto', interpolation = 'nearest', cmap=plt.cm.CMRmap)
			mb = ax2.imshow(b.values, aspect='auto', interpolation = 'nearest', cmap=plt.cm.CMRmap)
			
#			fig.colorbar(ma, ax=ax1)
			fig.colorbar(mb, ax=ax2)
			plt.show()
				

			

	elif (mode == 'XMCD' or mode == 'XAS'): 
	# XMCD ###################################################################

		element_V 	     = record()
		element_V.element    = "V"
		element_V.pre_edge   = (300,512)
		element_V.E_lower    = 0
		element_V.E_upper    = 600

		element_Fe 	     = record()
		element_Fe.element   = "Fe" 
		element_Fe.pre_edge  = (600,702)
		element_Fe.E_lower   = 600
		element_Fe.E_upper   = 760
		
		element_Co = record()
		element_Co.element   = "Co"			
		element_Co.pre_edge  = (760,774)
		element_Co.E_lower   = 760
		element_Co.E_upper   = 900	

		element_N = record()
		element_N.element   = "N"			
		element_N.pre_edge  = (380,390)
		element_N.E_lower   = 380
		element_N.E_upper   = 410	

		element_O = record()
		element_O.element   = "O"			
		element_O.pre_edge  = (510,525)
		element_O.E_lower   = 510
		element_O.E_upper   = 552	

		elementAll = record()
		elementAll.element   = "all"			
		elementAll.pre_edge  = (0,1000)
		elementAll.E_lower   = 0
		elementAll.E_upper   = 2000		


# 	SET ELEMENT HERE TO DEFINE ENERGY RANGE
		element = element_O

		xmcd_frame = treat_xmcd(sg, element)
		print "xmcd_frame"
		print xmcd_frame
		xmcd_frame_group = xmcd_frame.groupby(group_columns)

		
		plot_it = True
#		plot_it = False
		if plot_it:
	#		xmcd_frame = average_xmcd(xmcd_frame, group_columns)
			print "The frame", xmcd_frame
			x_plotnumber, y_plotnumber	=	get_plot_dimensions(scandata_f, group_columns)
#			plot_xmcd_frame(xmcd_frame, group_columns, x_plotnumber , y_plotnumber, mode=mode)
			plot_xmcd_frame(xmcd_frame, group_columns, x_plotnumber , y_plotnumber, mode=mode)

	
		plot_mapstyle = True
#		plot_mapstyle = False
		if plot_mapstyle:

			map_figure 	= plt.figure()
			map_ax1		= map_figure.add_subplot(111)
#			map_ax2		= map_figure.add_subplot(122)
			xC = xmcd_frame['Energy'] 
			yC = xmcd_frame['Z'] 
			zC = xmcd_frame['XAS'] 
			plot_map(xC,yC,zC, map_figure, map_ax1, title = "XAS intensity map "+element.element+" "+str(filenames)+".pdf", hide_points=True)

#			data_array = get_array(xmcd_frame)
#			plt.imshow(data_array)
			plt.show()


		export_CombiView = True
#		export_CombiView = False
		if export_CombiView:
			keys = xmcd_frame_group.groups.keys()
			keys.sort(sort_zy_lines)
			print keys

#			filenamePOS  	= "Overview "+element.element+" "+str(filenames)+"_POS.txt"
#			filenameXAS  	= "Overview "+element.element+" "+str(filenames)+"_XAS.txt"
#			filenameXMCD 	= "Overview "+element.element+" "+str(filenames)+"_XMCD.txt"
			filenamePOS  	= "Overview "+element.element+"_POS.txt"
			filenameXAS  	= "Overview "+element.element+"_XAS.txt"
			filenameXMCD 	= "Overview "+element.element+"_XMCD.txt"
			filePOS		= open(filenamePOS, "w")
			fileXAS		= open(filenameXAS, "w")
			fileXMCD	= open(filenameXMCD,"w")
			mspA.exportCombiView(xmcd_frame_group, keys, group_columns , filenamePOS, filenameXAS, filenameXMCD)


#		distill_something  = True
		distill_something  = False
		if distill_something:
			

			keys = xmcd_frame_group.groups.keys()
			keys.sort(sort_zy_lines)


#			filter the bad datapoints
			badkeys = []
			badkeys = [
					(71.900000000000006, -7.7999999999999998), # Co
					(77.900000000000006, -6.0), # Fe, Co
#					(83.900000000000006, -7.0), # Fe 
					(75.900000000000006, -6.0)  # Co
					]
			badkeys = [
#					(1, 0.40000000000000036), # Co
#					(7, 4), # Fe, Co
#					(5, 4), # Co (evtl)
					(13, 3), # Fe, Co, V
					(13, 2), # Fe, V
					(11, 2)  # V
					]
			f = lambda x: (x not in badkeys)


#			Composition
			xC,yC,zC = mspA.distill_scalingXAS(xmcd_frame_group, filter(f, keys) , group_columns)
#			pickle.dump((xC,yC,zC), open("XAS_"+element.element+"_Map.pickle", "wb"))

#			XMCD
			xX,yX,zX = mspA.distill_maxXMCD(xmcd_frame_group, filter(f, keys) , group_columns)
#			pickle.dump((xX,yX,zX), open("XMCD_"+element.element+"_Map.pickle", "wb"))

			map_figure 	= plt.figure()
			map_ax1		= map_figure.add_subplot(121)
			map_ax2		= map_figure.add_subplot(122)
			plot_map(xC,yC,zC, map_figure, map_ax1, title = "XAS intensity map "+element.element+" "+str(filenames)+".pdf")
			plot_map(xX,yX,zX, map_figure, map_ax2, title = "XAS intensity map "+element.element+" "+str(filenames)+".pdf")
			plt.show()

	else:
		print "usage: python scriptname mode filename"
		print "modes can be XAS XMCD MH"


#	plt.title( str(filename.split("/")[-2:])+ ' '+element)
#	temperature = str(sg.get_group(k)["Temperature Controller A"].mean().round())
#	fig.suptitle( str(filename)+ ' '+element.element+' '+temperature+"K", fontsize=32)

#	plt.show()

	save_pdf = False
#	save_pdf = True
	if save_pdf:
	#	pdfname	= "_".join(( 'Overview', filename.split('.')[0], mode, element.element+'.pdf' ))
		pdfname = "Overview-"+element.element+"-xmcd_N_dummy.pdf"
		fig.savefig(pdfname, dpi=150, facecolor='w', edgecolor='w',
	        	orientation='portrait', papertype=None, format=None,
	        	transparent=False, bbox_inches=None, pad_inches=0.1,
		        frameon=None)

#	extract = pd.DataFrame( {	"T" : t,
#					"m" : m,
#					"Y" : y,
#					"Z" : z})
#	plt.figure()
#	for i, group in extract.groupby(('Z','Y')):
#		group.plot(x='T', y='m', c=scipy.random.rand(3,1), label=str(i))
#	plt.legend()
#	ax1 = plt.gca()
#	box = ax1.get_position()
#	ax1.set_position([box.x0, box.y0, box.width * 0.66, box.height])
#	# Put a legend to the right of the current axis
#	ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#
#	plt.show()
#	extract.to_csv("mag_temp.csv")
else:
	print "not importing a file"
