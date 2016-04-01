#!/usr/bin/python

#
#   Plots MH loops and XAS. XMCD spectra from BL 6.3.1 This is not even an alpha Version (no pun intended...)
#   atndiaye@lbl.gov 
#
import pandas as pd
import matplotlib.pylab as plt
import scipy
import sys
import os.path
import numpy as np
import mspScan as msp


def get_xmcd(scandata):
	scan_plus	= scandata[scandata["Magnet Field"]<=0] 
	scan_minus 	= scandata[scandata["Magnet Field"]>0]
	xas_plus	= scan_plus ["I_Norm0"].values 
	xas_minus	= scan_minus["I_Norm0"].values
	
	if xas_plus.shape != xas_minus.shape:
		print ("No xas_plus and xas_minus not equal. Using xas_plus. Shapes : {0}, {1}".format(xas_plus.shape, xas_minus.shape) )
		xas_minus = xas_plus

	

	xas 		=  (xas_plus + xas_minus)/2
	xmcd		=  -(xas_plus - xas_minus)
	energy		= scan_plus["Energy"].values

	result_data_f	= pd.DataFrame( {	"Energy":	energy,
						"XAS":		xas,
						"XAS+":		xas_plus,
						"XAS-":		xas_minus,
						"XMCD":		xmcd	}
					)
	return result_data_f

def integrate_xmcd(xmcd_data):
	distances		= xmcd_data['Energy'].values[:-1] - xmcd_data['Energy'].values[1:] 
	areas 			= (xmcd_data['XMCD'][:-1] + xmcd_data['XMCD'][1:])/2 * distances
	xmcd_data['integral']	= areas.cumsum()
	return xmcd_data



def removeBG(xmcd_data, pre_edge = (760,772), post_edge= (797,840)):
	"""Should remove the bg and write everything into the dataframe, then we can plot it nicely later"""

	
	# Preedge to 1
	preedge   = xmcd_data[ (xmcd_data["Energy"]>pre_edge[0]) & (xmcd_data["Energy"]<pre_edge[1]) ]
	xmcd_data["XAS+"] /=    preedge["XAS+"].mean()
	xmcd_data["XAS-"] /=    preedge["XAS-"].mean()

	# Calculate XAS and XMCD
	xmcd_data["XAS"]   =   (xmcd_data["XAS+"] + xmcd_data["XAS-"])/2
	xmcd_data["XMCD"]  =   -(xmcd_data["XAS+"] - xmcd_data["XAS-"])

	# Normalize
	factor 		   = (xmcd_data["XAS"].max()-1)
	xmcd_data["XAS"]   = (xmcd_data["XAS"] - 1) / factor
	xmcd_data["XMCD"] /= factor
	xmcd_data["Factor"]= factor
	print factor

	return (xmcd_data)

def removeV_BG(xmcd_data, pre_edge = (760,772), post_edge= (797,840)):
	"""Should remove the bg and write everything into the dataframe, then we can plot it nicely later"""

	
#	# Preedge to 1
#	preedge   	   = xmcd_data[ (xmcd_data["Energy"] > pre_edge[0]) & (xmcd_data["Energy"]<pre_edge[1]) ]
#	xmcd_data["XAS+"] /=    preedge["XAS+"].mean()
#	xmcd_data["XAS-"] /=    preedge["XAS-"].mean()

	# Update preedge and fit polynomial
	preedge   	   = xmcd_data[ (xmcd_data["Energy"] > pre_edge[0]) & (xmcd_data["Energy"]<pre_edge[1]) ]
	bg_poly_plus_c	   = np.polyfit(preedge["Energy"].values,  preedge["XAS+"].values,  1)
	bg_poly_minus_c	   = np.polyfit(preedge["Energy"],  preedge["XAS-"], deg = 1)
	bg_poly_plus	   = np.poly1d(bg_poly_plus_c)
	bg_poly_minus	   = np.poly1d(bg_poly_minus_c)




	xmcd_data["XAS+"] -= bg_poly_plus( xmcd_data["Energy"] )
	xmcd_data["XAS-"] -= bg_poly_minus(xmcd_data["Energy"] )
	
#	######### Vanadium specific #####################
#	end_of_preedge 	   = ((xmcd_data["Energy"] >= 508) & (xmcd_data["Energy"] <= 509))
#	L3_energies	   = ((xmcd_data["Energy"] >= 510) & (xmcd_data["Energy"] <= 517.5))
#	bg_const	   = (bg_poly_plus(xmcd_data["Energy"][end_of_preedge] ).mean() + bg_poly_minus(xmcd_data["Energy"][end_of_preedge] ).mean() ) / 2
#
	bg_const	   = (bg_poly_plus(xmcd_data["Energy"]).mean() + bg_poly_minus(xmcd_data["Energy"] ).mean() ) / 2


	# set preedge to 1 based on subtracted BG
	xmcd_data["XAS+"]  = (xmcd_data["XAS+"] / bg_const) 
	xmcd_data["XAS-"]  = (xmcd_data["XAS+"] / bg_const)
	

	# Calculate XAS and XMCD
	xmcd_data["XAS"]   =   (xmcd_data["XAS+"] + xmcd_data["XAS-"])/2
	xmcd_data["XMCD"]  =   (xmcd_data["XAS+"] - xmcd_data["XAS-"])

	# Normalize
#	factor 		   = xmcd_data["XAS"][L3_energies].max()
	factor 		   = xmcd_data["XAS"].max() # *0+1
	xmcd_data["XAS"]  /= factor
	xmcd_data["XMCD"] /= factor
	xmcd_data["Factor"]= factor
#	print factor
		
	return (xmcd_data)



def distill_maxXMCD(xmcd_frame, keys , groupColumns):
	x,y,z = list(), list(), list()
	for g in keys:
		if np.isfinite( abs(xmcd_frame.get_group(g)["XMCD"]).max()):
# 			print g
			x.append(g[1])
			y.append(g[0])
			z.append(abs(xmcd_frame.get_group(g)["XMCD"]).max())
	return x,y,z

def distill_scalingXAS(xmcd_frame, keys , groupColumns):
	x,y,z = list(), list(), list()
	for g in keys:
		if np.isfinite( abs(xmcd_frame.get_group(g)["Factor"]).max()):
# 			print g
			x.append(g[1])
			y.append(g[0])
			z.append(abs(xmcd_frame.get_group(g)["Factor"]).max())
	return x,y,z



def distill_maxXAS(xmcd_frame, keys , groupColumns):
	x,y,z = list(), list(), list()
	for g in keys:
		scan_data = xmcd_frame.get_group(g)
	
#		if np.isfinite( abs(scan_data.where["XAS"]).max()):
		if np.isfinite( abs(scan_data[scan_data["Energy"]<520]["XAS"]).max()):
			x.append(g[1])
			y.append(g[0])
			z.append( abs(scan_data[scan_data["Energy"]<519.]["XAS"]).max())
			print g,'\t', z[-1]
		else:
			x.append(g[1])
			y.append(g[0])
			z.append(-1)
	return x,y,z



def exportCombiView(xmcd_frame_group, keys ,groupColumns, filenamePOS, filenameXAS, filenameXMCD):
	filePOS		  = open(filenamePOS, "w")
	fileXAS		  = open(filenameXAS, "w")
	fileXMCD	  = open(filenameXMCD,"w")	
	delim		  = '\t'
	EOL		  = '\r\n'

	#Header:
	filePOS.write(  delim.join(groupColumns)		+ EOL)
	fileXAS.write(  delim.join(xmcd_frame_group.get_group(keys[0])['Energy'].astype(str))	+ EOL)
	fileXMCD.write( delim.join(xmcd_frame_group.get_group(keys[0])['Energy'].astype(str))	+ EOL)
		
	for g in keys:
		print g
		scan_data = xmcd_frame_group.get_group(g)
		
		# Only if keys are groups: 
		# filePOS.write( delim.join([str(i) for i in g])		+ EOL)
		# Else if keys are singe things (e.g. a string or number):
		filePOS.write( str(g) + EOL)
		fileXAS.write( delim.join(scan_data['XAS' ].astype(str))	+ EOL)
		fileXMCD.write(delim.join(scan_data['XMCD'].astype(str))	+ EOL)
		
	filePOS.close()
	fileXAS.close()
	fileXMCD.close()


def plot_xmcd(xmcd_data, ax1=None, ax2=None, label = None, legend='full'):

	if legend == 'short':
		legend_label = str(label) # '%5.1f'%label[0] +' '+ '%5.1f'%label[1]
		if not ax2 is None: ax2.plot(xmcd_data["Energy"], xmcd_data["XMCD"], 	linewidth = 2, color = 'red')
		if not ax1 is None: ax1.plot(xmcd_data["Energy"], xmcd_data["XAS" ] , 	linewidth = 2 ,color = 'blue')

		bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.7)
		if not ax2 is None: ax2.text(0.9, 0.9, legend_label, transform=ax2.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='right', bbox=bbox_props)
		elif not ax1 is None: ax1.text(0.9, 0.9, legend_label, transform=ax1.transAxes, fontsize=8, verticalalignment='top', horizontalalignment='right', bbox=bbox_props)

		
	elif legend == 'full' or True:
		legend_label = '%5.3f'%label[0] +' '+ '%5.3f'%label[1]
		if not ax2 is None: ax2.plot(xmcd_data["Energy"], xmcd_data["XMCD"], 		linewidth = 2, color = 'red',  label = "XMCD " + str(legend_label))
		if not ax2 is None: ax2.plot([], [] , 						linewidth = 2 ,color = 'blue', label = "XAS  " + str(legend_label))
		if not ax1 is None: ax1.plot(xmcd_data["Energy"], xmcd_data["XAS" ] , 		linewidth = 2 ,color = 'blue', label = "XAS  " + str(legend_label))
		if not ax1 is None: ax1.plot([], [] , 						linewidth = 2 ,color = 'red',  label = "XCMD " + str(legend_label))
		if not ax2 is None: ax2.legend(fancybox=True, framealpha=0.7, fontsize=8)
		

	if not ax1 is None: 
		for tl in ax1.get_yticklabels():
			tl.set_color('blue')
	if not ax2 is None: 
		for tl in ax2.get_yticklabels():
			tl.set_color('red')
	return ax1,ax2








def get_mh_data(scandata, E1, E2, element="", ax1=None, label = None, digits=1):
	print "The Energies I seek:", round(E1), round(E2)
	print "The energies I got: ", (scandata["Energy"]).round().unique()

	scan_E1		= scandata[(scandata["Energy"]).round() == round(E1)] 
	scan_E2		= scandata[(scandata["Energy"]).round() == round(E2)] 
	
	# print "The scans which came out of E1 and E1", round(E1), round(E2)
	# print scan_E1.count(), scan_E2.count() 
	
	I_E1		=  scan_E1["I_Norm0"].values
	I_E2		=  scan_E2["I_Norm0"].values
	
	if I_E1.shape == I_E2.shape :
		loop 		= ((I_E1 - I_E2)/ (I_E1 + I_E2))#[1604:]
		field		= scan_E1["Magnet Field"].values#[1604:]
#		field		= scan_E1["Magnet Power"].values#[1604:]
		print "Loop points",loop.shape
	else:
		print "No equal shapes", I_E1.shape, I_E2.shape 
		loop = []
		field = []

	
#	#Cosmetics:
#	F = len(field)
#	if F>1:
#		drift = np.arange(-F/2,F/2)*(loop[-1] - loop[1])/(F-1)
#		loop -= drift
#	field = np.array(field)[range(2,F/2)+range(F/2+2,F)]
#	loop  = np.array(loop)[range(2,F/2)+range(F/2+2,F)]

#	#Normalize 
#	temp 		= abs(loop)
#	temp.sort()
#	print "Tempshape:", temp.shape
#	if int(temp.shape[0]*0.9) < temp.shape[0] :
#		loop /= 	scipy.median(temp[int(temp.shape[0]*0.9)])
#	loop -= 	scipy.median(loop)
	result_data_f	= pd.DataFrame( {	"Magnet Field":	field,
						"Asymmetry":	loop ,
						"Element":	element}
					)

	# Subtract previous
	result_data_f["branch"]				  	= result_data_f['Magnet Field'] - result_data_f['Magnet Field'].shift(1) 
	# for first point of if field has not changed, it's the next branch already: Sot to same brach value as next
	result_data_f["branch"].values[0] 			= result_data_f["branch"].values[1]
	result_data_f["branch"][result_data_f["branch"]==0] 	= result_data_f["branch"][result_data_f["branch"].shift(-1)==0]
	
	print "MH LOOP"
	print result_data_f #['branch']

	return result_data_f

def process_mh_data(mh_data):
#	factor = 1
	# normalization (This is not a good spot)
# 	factor 		= 2/(mh_data['Asymmetry'].max()-mh_data['Asymmetry'].min())
	factor		= 1/mh_data['Asymmetry'].max()
#	mh_data['Asymmetry'] -  mh_data['Asymmetry'].median()) * facto
#	mh_data['Asymmetry'] - mh_data['Asymmetry'].max()
	mh_data['Asymmetry'] = mh_data['Asymmetry'] * factor					
	return(mh_data)

def plot_mh(mh_data, ax1=None, label = None, linewidth = 1):
#	print "MH DATA"
#	print mh_data

	if (mh_data['Magnet Field'].shape) != (mh_data['Asymmetry'].shape):
		print "Different shape of field and asymmetry!",
		print  mh_data['Magnet Field'].shape, mh_data['Asymmetry'].shape


#	legend_label = '%5.3f'%label[0] +' '+ '%5.3f'%label[1] + ' '+ label[2]
	legend_label = label
	ax1.plot(mh_data['Magnet Field'],  mh_data['Asymmetry'], linewidth = linewidth,alpha=1, label = legend_label)
#	ax1.plot(mh_data['Magnet Field'],  mh_data['branch'], linewidth = linewidth,alpha=1, label = legend_label)
	return()


