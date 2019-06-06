###### SAO Preprocessing code
###### Developed by G.Lim (lim9gu@gmail.com)
### run /data3/SAO1m/code/SAO_process.py
### or import SAO_process as SAO
#=================================================================
def CCDinfo(camera='STX16803'):
	"""
	###### SAO Preprocessing code Ver 0.6
	###### developed by G.Lim (lim9gu@gmail.com)

	*** Prerequisites ***

	- python 2.7 based
	- numpy
	- astropy
	- pyraf
	- iraf (> 2.14)
	- SExtractor
	- Astrometry.net
	- hotpants
	- astroquery

	1. Description 
	: This function prints SAO 1-m information. 
	
	2. Usage 
	>>> CCDinfo()

	3. History
	2018.12.20 G. Lim created. 
	2019.02.08 Rdnoise information is added.
	2019.04.18 FLI Kepler CCD information is added by G.Lim
	"""

	if camera == 'STX16803':
		print(""" =======     SAO 1m SBIG STX-16803     =======\n Gain 		= 1.3600000143051147\n Rdnoise 	= 9.\n FOV 		= 21.2 x 21.2 arcmin^2\n Pixel scale	= 0.311 arcsec/pixel\n ==============================================\n""")
	elif camera == 'Kepler' :
		print(""" =======     SAO 1m FLI Kepler    ======\n Gain 		= 1.5\n Rdnoise    = ? \n FOV     = ? x ? arcmin^2\n Pixel scale	= ? arcsec/pixel\n==============================================\n""")
#=================================================================
def fileset(camera = 'STX16803'):
	"""
	1. Description
	: This function classifies calibration frames (bias, dark, skyflat, domeflat), science frames making each directory. Using other calibration frames will be added soon. If no calibration frames today, copy recent master calibration frames on current directory.
	
	2. Usage
	>>> fileset()

	3. History
	2018. 12. 28 G. Lim created.
	2019. 01. 23 Finding recent master calibration frame is added by G. Lim
	2019. 01. 24 Image size cut section is added by G. Lim
				 : Classify binned data (1x1, 2x2) with their size. Use 1x1 binned images which have 33.5MB size. 2x2 binned images are moved to bin2 directory. 
	2019. 02. 07 Loading recent calibration images is added by G. Lim
	2019. 03. 11 Image name change code is added. 
				 Domeflat, skyflat is standard naming.
				 DOMEFLAT, domeflat, SKYFLAT, Skyflat is not compatible.
	2019. 06. 06 bias, dark, flat, sci images are classified as header information 'IMAGETYP'. and if there are no calibrations, the code calls recent master images. This function had an issue that cannot call exptime and filter which don't exist. 'hdrcheck' function is added to calculate MJD and to check if object name is entered correctly by comparing CRVAL1, CRVAL2 with IMSNG galaxy catalog. 
	"""	
	import glob
	import os, sys
	import numpy as np
	from astropy.io import fits
	from astropy.table import Table
	# Image size cut
	allimage = glob.glob('*-00*.fit')
	allimage.sort()
	# Binning classification
	xbin = []
	ybin = []
	for im in allimage :
		hdr = fits.getheader(im)
		xbin_dum = hdr['xbinning']
		ybin_dum = hdr['ybinning']
		xbin.append(xbin_dum)
		ybin.append(ybin_dum)
	cat = Table({'allimage' : allimage, 'xbinning': xbin, 'ybinning' : ybin}, names = ['allimage', 'xbinning', 'ybinning'])
	bin1 = np.where((cat['xbinning'] == 1) & (cat['ybinning'] == 1))[0]
	bin2 = np.where((cat['xbinning'] == 2) & (cat['ybinning'] == 2))[0]
	if len(bin2) !=0 :
		print('There are 2x2 binned data...')
		os.system('mkdir bin2')
		os.system('mv '+" ".join(cat['allimage'][bin2])+' ./bin2')	
		print('Find 2x2 binned data at bin2 folder.')
	# Name change
	os.system('rename zero bias *tion*zero.fit')
	os.system('rename Bias bias *tion*Bias.fit')
	os.system('rename Cal cal Cal*bias.fit')
	os.system('rename D dk *tion*D*.fit')
	#os.system('rename d dk *tion*?d*.fit')
	os.system('rename dark dk *tion*dark*.fit')
	os.system('rename Cal cal Cal*dk*.fit')
	#calibrations = glob.glob('*.fit')
	#calibrations.sort()
	imbin1 = cat[bin1]['allimage']
	bias, dark, flat, sci = [], [], [], []
	for i in xrange(len(imbin1)):
		hdr      = fits.getheader(imbin1[i])
		IMAGETYP = hdr['IMAGETYP']
		if IMAGETYP == 'Bias Frame' :
			bias.append(imbin1[i])
		elif IMAGETYP == 'Dark Frame' :
			dark.append(imbin1[i])
		elif IMAGETYP == 'Flat Field' :
			flat.append(imbin1[i])
		elif IMAGETYP == 'Light Frame' :
			sci.append(imbin1[i])
	print('All images are classified.')
	bias.sort()
	dark.sort()
	flat.sort()
	sci.sort()
	
	# Bias classification
	if  bias == [] :
		print('No bias today or bias directory is already created!')
		masterbias_loc = '/data3/SAO1m/red/'+camera+'/masterbias/'
		os.system('ls '+masterbias_loc)
		masterbias_list = glob.glob(masterbias_loc+'*.fits')
		masterbias_list.sort()
		recent_masterbias = masterbias_list[-1]
		print('Recent master bias of '+recent_masterbias+' is found.')
		os.system('cp '+recent_masterbias+' ./')
		print(recent_masterbias+' is copied.')
	elif bias != [] :
		os.system('mkdir bias')
		biasjoin = " ".join(bias)
		os.system('mv '+biasjoin+' ./bias')

	# Dark exp time classification
	masterdark_loc = '/data3/SAO1m/red/'+camera+'/masterdark/'
	os.system('ls '+masterdark_loc)

	# Exptime of sci images
	sciexptime = []
	for i in xrange(len(sci)) :
		scihdr = fits.getheader(sci[i])
		sciexptime.append(scihdr['exptime'])
		sciexpset = set(sciexptime)
		sciexptime = list(sorted(sciexpset))	

	for exp in sciexptime:
			darkexp = glob.glob('cal*dk'+str(int(exp))+'.fit')
			if len(darkexp) != 0:
				os.system('mkdir dark')
				darkjoin = " ".join(darkexp)
				os.system('mv '+darkjoin+' ./dark')	
			elif len(darkexp) == 0 :
				dark_dum_name = '2*_dark'+str(int(exp))+'.fits'
				print('Find master dark of '+dark_dum_name+'...')
				dark_dum = glob.glob(dark_dum_name)
				if len(dark_dum) == 1 :
					print('I found '+dark_dum_name+'!')
					print('Work with this masterdark. STOP.')
				elif len(dark_dum) == 0 :
					print('No dark images of '+str(int(exp)))
					masterdark_list = glob.glob(masterdark_loc+'2*dark'+str(int(exp))+'.fits')
					masterdark_list.sort()
					print(masterdark_list)
					recent_masterdark = masterdark_list[-1]
					if len([recent_masterdark]) == 1 :
						print('Recent master dark of '+recent_masterdark+' is found.')
						os.system('cp '+recent_masterdark+' ')
						print(recent_masterdark+' is copied.')
					elif len([recent_masterdark]) == 0 :
						print('No masterdark. Stop here :(')	

	# Flat classification
	os.system('rename Skyflat skyflat Skyflat*.fit')
	os.system('rename SKYFLAT skyflat SKYFLAT*.fit')
	os.system('rename SkyFlat skyflat SkyFlat*.fit')
	os.system('rename domeflat Domeflat domeflat*.fit')
	os.system('rename DomeFlat Domeflat DomeFlat*.fit')
	os.system('rename DOMEFLAT Domeflat domeflat*.fit')
	os.system('rename t-000 t_000 Domeflat*.fit')
	os.system('rename flat_ flat- *flat*.fit')

	# Filter of science images

	scifilter = []
	i=0
	for i in range(len(sci)):
		scihdr = fits.getheader(sci[i])
		scifilter.append(scihdr['filter'])
	scifilterset = set(scifilter)
	sciinfilter = list(sorted(scifilterset))
	print(sciinfilter, 'filter set today.')

	for band in sciinfilter :
		masterflat_loc = '/data3/SAO1m/red/'+camera+'/masterflat_'+band+'/'
		os.system('ls '+masterflat_loc)
		flatband = glob.glob('*flat*'+band+'.fit')
		if len(flatband) != 0 :
			if 'skyflat' in flatband[0].split('-') :
				skyflat = flatband
				os.system('mkdir skyflat')
				skyflatjoin = " ".join(skyflat)
				os.system('mv '+skyflatjoin+' ./skyflat')
			elif 'Domeflat' in flatband[0].split('-') :
				domeflat = flatband
				os.system('mkdir domeflat')
				domeflatjoin = " ".join(domeflat)
				os.system('mv '+domeflatjoin+' ./domeflat')			
		elif len(flatband) == 0 :
			print('No flat frames today.')
			masterflat_list = glob.glob(masterflat_loc+'2*_n'+band+'flat.*.fits')
			masterflat_list.sort()
			print(masterflat_list)
			recent_masterflat = masterflat_list[-1]
			if len([recent_masterflat]) == 1:
				print('Recent master flat of '+recent_masterflat+' is found.')
				os.system('cp '+recent_masterflat+' ./')
			elif len([recent_masterflat]) == 0:
				print('No masterflat. Stop here :(')
	print('File setting is done.')
	return bias, dark, flat, sci
#=================================================================
def biascom(camera='STX16803'):
	"""
	1. Description 
	: This function makes a master bias image of SAO 1-m using Pyraf. Put bias images to 201XXXXX/bias/ directory. Run this code on 201XXXXX directory, then pyraf chdir task enter bias directory and makes process. Output image is zero.fits and it will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 
	
	2. Usage 
	: Start on 2018XXXX directory. Make bias directory which contains each bias frame. Naming of each bias images should be cal*bias.fit. Then just use SAO_biascom()
	
	>>> SAO_biascom()

	3. History
	2018.03    Created by G.Lim.
	2018.12.17 Edited by G.Lim. Define SAO_biascom function. 
	2019.02.07 Assign archive of masterbias in each date by G. Lim
	"""
	import glob
	import os, sys
	import numpy as np
	from pyraf import iraf
	from astropy.io import fits
	curdir = os.getcwd()
	curdate = curdir.split('/')[-1]
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	iraf.chdir('./bias')
	input_name = 'bias.list'
	output_name = curdate+'_zero.fits'
	#os.system('ls cal*bias.fit > '+input_name)
	calibrations = glob.glob('cal*.fit')
	f    = open(input_name,'w+')
	for i in xrange(len(calibrations)):
		hdr      = fits.getheader(calibrations[i])
		IMAGETYP = hdr['IMAGETYP']
		if IMAGETYP == 'Bias Frame' :	
			f.write(calibrations[i]+'\n')
	f.close()
	print('Zerocombine is running...')
	iraf.zerocombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
	print('Output master '+output_name+' is created.')
	os.system('cp '+output_name+' ../')
	os.system('cp '+output_name+' /data3/SAO1m/red/'+camera+'/masterbias/')
	iraf.chdir('../')
	iraf.dir('.')
#=================================================================
def darkcom(camera='STX16803') :
	"""
	1. Description 
	: This function makes a master dark image of SAO 1-m using Pyraf. Put dark images to 201XXXXX/dark/ directory. Run this code on 201XXXXX directory, then pyraf chdir task enter dark directory and makes process. Output image is darkXXX.fits (XXX is exposure time. This function will classify each exposure of dark frames!) and it will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 

	2. Usage
	: Start on 2018XXXX directory. Make dark directory which contains each dark frame. Naming of each dark image should be cal*dk*.fit. Then just use SAO_darkcom() 

	3. History
	2018.03    Created by G.Lim.
	2018.12.20 Edited by G.Lim. Define SAO_darkcom function.
	2019.02.07 Assign archive of masterdark in each date by G. Lim
	"""
	import glob
	import os, sys
	import numpy as np	
	from pyraf import iraf
	from astropy.io import fits
	curdir = os.getcwd()
	curdate = curdir.split('/')[-1]
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	iraf.chdir('./dark')
	dark = glob.glob('cal*dk*.fit')
	allexptime = []
	for i in xrange(len(dark)) :
		hdr = fits.getheader(dark[i])
		allexptime.append(hdr['exptime'])
	expset = set(allexptime)
	exptime = list(sorted(expset))
	i=0
	for i in xrange(len(exptime)) :
		print('Find images with exptime of '+str(exptime[i]))
		imlist = []
		for j in xrange(len(dark)) :
			hdr = fits.getheader(dark[j])
			if hdr['exptime'] == exptime[i] :
				imlist.append(dark[j])
			else :
				pass
		print(imlist)
		output_name = curdate+'_dark'+str(int(exptime[i]))+'.fits'
		input_name = output_name[:-5]+'.list'
		f=open(input_name,'w+')
		for k in xrange(len(imlist)) : 
			f.write(imlist[k]+'\n')
		f.close()
		print('Darkcombine is running...')
		iraf.imstat('@'+input_name)
		iraf.darkcombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='none', ccdtype='' )
		os.system('cp '+output_name+' ../')
		os.system('cp '+output_name+' /data3/SAO1m/red/'+camera+'masterdark/')
	os.system('rm d*.list')
	iraf.chdir('../')
	iraf.dir('.')
	print('Output master '+output_name+' is created.')
#=================================================================
def flatcom(flattype='sky', process='bias', camera='STX16803') :
	"""
	1. Description 
	: This function makes master-normalised images of SAO 1-m using Pyraf. Put flat images to 201XXXXX/flat/ directory. Run this code on 201XXXXX directory, then pyraf chdir task enter skyflat, domeflat directory and makes process. Keyword 'process' will decide if bias or dark subtraction is required or not. If process = bias, bias frame will be subtracted. If process = dark, dark frame will be subtracted (This function is for reducing domeflats which have long exposure time.). And then flatcombine and normalizing will be performed. Output image is nflatX.YYY.fits (X is filter and YYY is sky or dome. This function will classify each exposure of frames!) and they will be copied on upper directory. Due to iraf.chdir task, you should reset python when this code is finished in order to make confusion of current directory between iraf and python! 

	2. Usage
	: Start on 2018XXXX directory. Make skyflat or domeflat directory which contains each flat frame. Naming of each flat image should be *flat*.fit. And domeflat naming is Domeflat*.fit. Then just use SAO_flatcom(). 
	>>> SAO_flatcom('sky') --> Use skyflat
	>>> SAO_flatcom('dome') --> Use domeflat
	* default configuration is sky using bias processing.

	3. History
	2018.03    Created by G. Lim.
	2018.12.20 Edited by G. Lim. Define SAO_flatcom function.
	2018.12.28 Edited by G. Lim. Add bias or dark keyword. Join function is used when performing imarith, combine tasks.
	2019.02.07 Assign archive of masterflat in each date by G. Lim
	"""
	import glob
	import os, sys
	import itertools
	import numpy as np
	from pyraf import iraf
	from astropy.io import fits
	from astropy.io import ascii
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	curdir = os.getcwd()
	curdate = curdir.split('/')[-1]
	#flattype = sys.argv[1] # dome or sky
	#if flattype == 'dome_Hal' :
	#	iraf.chdir('./domeflat_Hal')
	#elif flattype == 'dome_LED' :
	#	iraf.chdir('./domeflat_LED')
	if flattype == 'dome' :
		iraf.chdir('./domeflat')
	elif flattype == 'sky' :
		iraf.chdir('./skyflat')
	flat = glob.glob('*flat*.fit')
	flat.sort()
	if process == 'bias' : # Bias subtraction : mainly skyflat
		input_name = 'flat.list'
		k=0
		f=open(input_name,'w+')
		for k in xrange(len(flat)) : 
			f.write(flat[k]+'\n')
		f.close()
		#flatjoin = ",".join(flat)
		#outflatjoin = ",z".join(flat)
		#outflatjoin = 'z'+outflatjoin
		print('zero subtraction with '+flattype+'flat images...')
		iraf.imarith(operand1='@'+input_name, op='-', operand2='../*_zero.fits', result='z@'+input_name)	
		#iraf.imarith(operand1=flatjoin, op='-', operand2='../zero.fits', result=outflatjoin)	
	elif process == 'dark' : # Dark subtraction : mainly domeflat
		allexptime = []
		i=0
		for i in xrange(len(flat)) :
			hdr = fits.getheader(flat[i])
			allexptime.append(hdr['exptime'])
		expset = set(allexptime)
		exptime = list(sorted(expset))
		i=0
		for i in xrange(len(exptime)) :
			print('Find images with exptime of '+str(int(exptime[i])))
			imlist = []
			j=0
			for j in xrange(len(flat)) :
				hdr = fits.getheader(flat[j])
				if hdr['exptime'] == exptime[i] :
					imlist.append(flat[j])
				else :
					pass
			print(imlist)
			imlist.sort()
			input_name = 'flat.list'
			k=0
			f=open(input_name,'w+')
			for k in xrange(len(imlist)) : 
				f.write(imlist[k]+'\n')
			f.close()		
			# use join
			#darkjoin = ",".join(imlist)
			#outdarkjoin = ",f".join(imlist)
			#outdarkjoin = 'f'+outdarkjoin
			iraf.imarith(operand1='@'+input_name, op='-', operand2='../*_dark'+str(int(exptime[i]))+'.fits', result='d@'+input_name)
			#iraf.imarith(operand1=darkjoin, op='-', operand2='../dark'+str(int(exptime[i]))+'.fits', result=outdarkjoin)
			print('Dark subtracted flat images are created.')
	# Flat combine
	if process == 'bias' :
		calflat = glob.glob('z*.fit')
	elif process == 'dark' :
		calflat = glob.glob('d*.fit')
	allfilter = []
	i=0
	for i in xrange(len(calflat)) :
		hdr = fits.getheader(calflat[i])
		allfilter.append(hdr['filter'])
	filterset = set(allfilter)
	infilter = list(sorted(filterset))
	i=0
	for i in xrange(len(infilter)) :
		print('Find images with filter of '+str(infilter[i]))
		imlist = []
		for j in xrange(len(calflat)) :
			hdr = fits.getheader(calflat[j])
			if hdr['filter'] == infilter[i] :
				imlist.append(calflat[j])
			else :
				pass
		print(imlist)
		imlist.sort()
		input_name = str(infilter[i])+'flat.list'
		k=0
		f=open(input_name,'w+')
		for k in xrange(len(imlist)) : 
			f.write(imlist[k]+'\n')
		f.close()
		#flatjoin = ",".join(imlist)
		output_name = input_name[:-5]+'.fits'
		iraf.flatcombine(input='@'+input_name, output=output_name, combine='median', reject='minmax', process='no', scale='mode', ccdtype='')
		#iraf.flatcombine(input=flatjoin, output=output_name, combine='median', reject='minmax', process='no', scale='mode', ccdtype='')
		print(output_name+' is created. Normalizing...')
		data, newhdr = fits.getdata(output_name, header=True)
		x = np.mean(data)
		nimage = data/x
		newflat_name = curdate+'_n'+str(infilter[i])+'flat.'+flattype+'.fits'
		fits.writeto(newflat_name, nimage, header=newhdr, overwrite=True)
		os.system('cp '+newflat_name+' ../')
		os.system('cp '+newflat_name+' /data3/SAO1m/red/'+camera+'/masterflat_'+infilter[i]+'/')
	print('Normalised master flats are created.')
	iraf.imstat(images='*n?flat.'+flattype+'.fits')
	os.system('rm *.list ?flat.fits')
	iraf.chdir('../')
	iraf.dir('./')
#=================================================================
def objpre(sci_list) :
	"""
	1. Description
	: This function applies master calibration images to science frames, including dark subtraction, flat fielding. STX16803 CCD has an issue of calibration frame when image type is set to bias and dark. So the CCD takes calibration images with H alpha light frame with exposure 0s as bias and >0s as dark frame.
 
	2. Usage
	: objpre('sky')

	3. History
	2018.03    Created by G. Lim
	2019.02.07 change name of master calibration frames in each date by G. Lim
	"""
	import glob
	import os, sys
	import itertools
	import numpy as np
	from pyraf import iraf
	from astropy.io import fits
	from astropy.io import ascii
	curdir = os.getcwd()
	curdate = curdir.split('/')[-1]
	iraf.noao()
	iraf.imred()
	iraf.ccdred()
	iraf.ccdred.setinst(instrume='camera', directo='/iraf/iraf/noao/imred/ccdred/ccddb/', query='q', review='no')
	#obj_list = '*-00*.fit'
	#obj = glob.glob(obj_list)
	obj = sci_list
	# dark subtraction
	allexptime = []
	i=0
	for i in xrange(len(obj)) :
		hdr = fits.getheader(obj[i])
		allexptime.append(hdr['exptime'])
	expset = set(allexptime)
	exptime = list(sorted(expset))
	i, j, k = 0, 0, 0
	for i in xrange(len(exptime)) :
		print('Find images with exptime of '+str(exptime[i]))
		imlist = []
		for j in xrange(len(obj)) :
			hdr = fits.getheader(obj[j])
			if hdr['exptime'] == exptime[i] :
				imlist.append(obj[j])
			else :
				pass
		print(imlist)
		imlist.sort()
		print('Creating object list for dark subtraction...')
		f = open("obj"+str(int(exptime[i]))+".list", 'w+')
		for im in xrange(len(imlist)) :
			f.write(imlist[im]+"\n")
		f.close()
		#imjoin = ",".join(imlist)
		#outimjoin = ",d".join(imlist)
		#outimjoin = 'd'+outimjoin
		print('dark subtraction with dark'+str(int(exptime[i]))+'.fits')
		input_dark = glob.glob('20*dark'+str(int(exptime[i]))+'.fits')[0]
		iraf.imarith(operand1 = '@obj'+str(int(exptime[i]))+'.list', op = '-', operand2 = input_dark, result = 'd@obj'+str(int(exptime[i]))+'.list')
		#iraf.imarith(operand1 = imjoin, op = '-', operand2 = './dark'+str(int(exptime[i]))+'.fits', result = outimjoin)
		print('dark subtracted object images are created.')
	# Flat fielding
	#dobj = glob.glob('d'+obj_list)
	dobj = ['d'+x for x in obj]
	dobj.sort()
	allfilter = []
	i=0
	for i in xrange(len(dobj)) :
		hdr = fits.getheader(dobj[i])
		allfilter.append(hdr['filter'])
	filterset = set(allfilter)
	infilter = list(sorted(filterset))
	i, j, k = 0, 0, 0
	for i in xrange(len(infilter)) :
		print('Find images with filter of '+str(infilter[i]))
		imlist = []
		imlist.sort()
		for j in xrange(len(dobj)) :
			hdr = fits.getheader(dobj[j])
			if hdr['filter'] == infilter[i] :
				imlist.append(dobj[j])
			else :
				pass
		print(imlist)

		g = open("obj"+infilter[i]+".list", 'w+')
		for im in xrange(len(imlist)) :
			g.write(imlist[im]+"\n")
		g.close()
		#dimjoin = ",".join(imlist)
		#doutimjoin = ",f".join(imlist)
		#doutimjoin = 'f'+doutimjoin
		print('Performing flat fielding...')
		#nflats = glob.glob('2*n'+str(infilter[i])+'flat.'+flattype+'.fits')[0]
		nflats = glob.glob('2*n'+str(infilter[i])+'flat.*.fits')[0]
		flattype = nflats.split('.')[1]
		print(nflats)
		#nflats = 'n'+str(infilter[i])+'flat.'+flattype+'.fits'
		iraf.imarith(operand1='@obj'+infilter[i]+'.list', op='/', operand2=nflats, result='f@obj'+infilter[i]+'.list')
		#iraf.imarith(operand1=dimjoin, op='/', operand2=nflats, result=doutimjoin)
	print('Flat fielding is finished. Check the images.')
#=================================================================
def astrometry() :
	"""
	1. Description 
	: Solving WCS coordinates using Astrometry.net software. For better performance in especially B band images, --use-sextractor mode is added. This mode needs SExtractor configuration files. So please posit configuration files for your working directory. cpulimit 300 is also added to prevent too long processing time for bad images.  

	2. Usage
	>>> SAO_astrometry()

	3. History
	2018.03    Created by G.Lim.
	2018.12.18 Edited by G.Lim. SExtractor mode is added.
	2018.12.21 Edited by G.Lim. Define SAO_astrometry function.
	"""
	import os,sys
	import glob
	import subprocess
	import numpy as np
	addlist = glob.glob('fd*.fit')
	addlist.sort()
	sexconfig = '/data3/SAO1m/code/astrom.config/astrometry.net.sex'
	print('Solving WCS using Astrometry.net...')
	for n in range(len(addlist)):
		com='solve-field '+addlist[n]+' --cpulimit 300 --overwrite --use-sextractor  --sextractor-config '+sexconfig+' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low 0.3 --scale-high 0.4 --no-remove-lines --uniformize 0 --no-plots  --new-fits a'+addlist[n]+' --temp-dir .\n'
		print(com)
		print(str(n)+' th of '+str(len(addlist)))
		os.system(com) 
	orinum = subprocess.check_output('ls fd*.fit | wc -l', shell=True)
	resnum = subprocess.check_output('ls afd*.fit | wc -l', shell=True)
	print("from "+str(orinum[:-1])+" files , "+str(resnum[:-1])+" files are solved.")
	print("All done.")
	os.system('rm tmp*')
	os.system('rm *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match ')
	print('Astrometry process is complete.')
#=================================================================
def hdrcheck(imlist_name):
	"""
	Check header and put galaxy name comparing to IMSNG catalog. 
	"""
	import os
	import sys
	import glob
	import astropy.units as u
	from astropy.io import ascii
	from astropy.io import fits
	from astropy.coordinates import SkyCoord

	qso_data_path = '/data3/SAO1m/obsdata/STX16803/'
	qso_red_path  = '/data3/SAO1m/red/STX16803/'

	all_catname = '/data3/IMSNG/alltarget.dat'
	all_cat     = ascii.read(all_catname) 
	ra, dec = all_cat['ra'], all_cat['dec']
	radeg, decdeg = [], []
	for i in range(len(all_cat)) :
		c   = SkyCoord(str(ra[i])+' '+str(dec[i]), unit=(u.hourangle, u.deg))
		radeg.append(c.ra.deg)
		decdeg.append(c.dec.deg)
	all_cat['radeg'] = radeg
	all_cat['decdeg'] = decdeg
	coo_all = SkyCoord(radeg, decdeg, unit=(u.deg,u.deg))

	imlist_name = 'a*.fit'
	imlist = glob.glob(imlist_name)
	imlist.sort()
	for i in range(len(imlist)) :
		inim            = imlist[i]
		print(inim)
		data, hdr       = fits.getdata(inim, header=True)
		CRVAL1          = hdr['CRVAL1']
		CRVAL2          = hdr['CRVAL2']
		jd              = hdr['jd']
		coo_target      = SkyCoord(CRVAL1, CRVAL2, unit=(u.deg, u.deg))
		indx, d2d, d3d  = coo_target.match_to_catalog_sky(coo_all)
		if indx.size == 0 :
			print('No matching. Maybe you obtained wrong field. OR Non IMSNG target.' )
			#os.system('mkdir bad')
			#os.system('mv '+inim+' .//')
		elif indx.size == 1 :
				obj  = all_cat[indx]['obj']
				mjd0 = 2400000.5
				mjd  = jd - mjd0
				print('======================================')
				print(obj)+ ' is matched.'
				print(str(round(d2d.arcmin[0],3))+ ' arcmin apart')
				print('======================================')
				hdr['object']   = obj
				hdr['mjd']      = mjd
				fits.writeto(inim, data, header=hdr, overwrite=True)
		#elif len(idx) > 1 :
		#		print('2 targets in 1 field.')
	print('Header info inspection is finished.')
#=================================================================
def fnamechange(camera='STX16803') :
	"""
	1. Description 
	: Change file name of WCS solving images (afd*.fits) using naming sequence following: Calib-SAO-CCD_name-OBJECT-UTDATE-UTSTART-FILTER-EXPTIME.fits. Then the code copies changed files to IMSNGgalaxies directory.

	2. Usage
	>>> SAO_fnamechange()

	3. History
	2018.03    Created by G.Lim.
	2018.12.21 Edited by G.Lim. Define SAO_fnamechange function.
	"""
	import glob
	import os, sys
	import subprocess
	import numpy as np 
	from astropy.io import fits
	lists = glob.glob('afd*.fit')
	lists.sort()
	for i in xrange(len(lists)):
		hdr = fits.getheader(lists[i])
		EXPTIME = int(hdr['exptime'])
		if camera == 'Kepler' :
			FILTER = 'R'
		elif camera == 'STX16803' :
			FILTER = hdr['filter']
		OBJECT = hdr['object']
		if (OBJECT == 'M51') | (OBJECT == 'M51a') :
			print('Object name is corrected.')
			OBJECT = 'M51A'
		UTDATE = hdr['date-obs'][0:10]
		UTSTART = hdr['date-obs'][11:19]
		newimage = 'Calib-SAO_'+camera+'-'+OBJECT+'-'+str(UTDATE[0:4])+str(UTDATE[5:7])+str(UTDATE[8:10])+'-'+str(UTSTART[0:2])+str(UTSTART[3:5])+str(UTSTART[6:8])+'-'+FILTER+'-'+str(EXPTIME)+'.fits'
		os.system('cp '+lists[i]+' '+newimage)
		print('Copy '+lists[i]+' to '+newimage+'.')
	#os.system('cp Cal*.fits /data3/IMSNG/IMSNGgalaxies/')
	print("Basic preprocessing of SAO 1-m is finished.")
#=================================================================
def filemove(camera='STX16803'):
	"""
	1. Description 
	: This code distributes all the calibrated images of SAO 1m data to IMSNGgalaxies/SAO directory, based on C. Choi's code.
	
	2. Usage
	: Run this code on '/data3/IMSNG/IMSNGgalaxies' location.
	>>> SAO_filemove() 

	3. History
	2018.12    Created by G.Lim 
	2018.01.24 Docstring is added by G. Lim
	"""
	import os
	import numpy as np
	import astropy.io.fits as fits
	import astropy.io.ascii as ascii
	loc = '/data3/IMSNG/IMSNGgalaxies/'
	os.system('ls Cal*.fits -l > currentdir.list')
	curdirlist = ascii.read('currentdir.list',data_start=0)
	dirpar     = curdirlist['col2']
	name       = curdirlist['col9']
	#dirnames   = name[np.where(dirpar >  1)]
	filenames  = name[np.where(dirpar == 1)]
	calframes  = []
	for i in filenames : 
		if (i[:9] =='Calib-SAO') & (i[-4:] == 'fits') : calframes.append(i)
	print(str(len(calframes))+' files exist.')
	for n in calframes:
		galname  = n.split('-')[2]
		print(galname)
		if galname == 'M51a' :
			galname = 'M51A'
		hdr = fits.getheader(n)
		if camera == 'STX16803' :
			infilter = hdr['filter']
		elif camera == 'Kepler' :
			infilter = 'R'
		makedir    = 'mkdir '+loc+galname
		os.system(makedir)
		makeobsdir = 'mkdir '+loc+galname+'/SAO/' # If IMSNG target
		os.system(makeobsdir)
		makecamdir = 'mkdir '+loc+galname+'/SAO/'+camera
		os.system(makecamdir)
		makefildir = 'mkdir '+loc+galname+'/SAO/'+camera+'/'+infilter+'/'
		os.system(makefildir)
		mvcommand='cp '+ n +' '+loc+galname+'/SAO/'+camera+'/'+infilter+'/'
		mvsubtract = 'cp hd'+n+' hc'+n+' '+loc+galname+'/SAO/'+camera+'/'+infilter+'/'
		os.system(mvcommand)
		os.system(mvsubtract)
		chmodcom = 'chmod 777 -R '+loc+galname+'/SAO/'
		os.system(chmodcom)
#=================================================================
def makedir(camera='STX16803'):
	"""
	1. Description
	: Make directory of camera name (STX16803 or Kepler)

	2. Usage
	>>> 

	3. History
	2019.04.21 : G. Lim created.
	"""
	import glob
	import os, sys
	workdir = '/data3/IMSNG/IMSNGgalaxies'
	curdir = os.getcwd()
	if curdir != workdir :
		print('Please run this code in '+workdir)
	elif curdir == workdir :
		folder = glob.glob('*')
		folder.sort()
		for i in range(len(folder)):
			os.chdir(folder[i])
			obs = glob.glob('*')
			obs.sort()
			if 'SAO' in obs :
				print('No camera dir. Make '+camera)
				os.chdir('SAO')
				os.system('mkdir '+'./'+ camera)
				os.system('mv Cal*.fits '+camera)
				#os.system('rsync -a B '+camera+'/B')
				os.chdir('../../')
			else :
				print(folder[i]+' has no SAO data.')
				os.chdir('../')
				pass
		print('Done.')
#=================================================================
def gregister(imlist_name) :
	"""
	1. Description 
	: 
	
	2. Usage
	>>> SAO_gregister() 

	3. History
	2019.04.17  Created by G. Lim
	"""
	import glob
	import alipy
	import os,sys
	import numpy as np
	from astropy.io import fits

	### Object classification
	# imlist_name = 'Cal*skysub.fits'
	lists      = glob.glob(imlist_name)
	lists.sort()
	split = lists[0].split('-')
	objects    = []
	for i in xrange(len(lists)) :
		hdr    = fits.getheader(lists[i])
		objects.append(hdr['object'])
	objectset  = set(objects)
	objectlist = list(sorted(objectset))

	def gregister(images):
		id     = alipy.ident.run(ref_image, images, visu=False)
		print("%20s : %20s, flux ratio %.2f" % (id[0].ukn.name, id[0].trans, id[0].medfluxratio))
		alipy.align.irafalign(id[0].ukn.filepath, id[0].uknmatchstars, id[0].refmatchstars, shape=outputshape, makepng=False)

	### Filter classification
	obj = 0
	for obj in xrange(len(objectlist)):
		object_name = objectlist[obj]
		if 'gregister' in split :
			image_list  = glob.glob(split[0]+'*'+object_name+'*gre*.fits')
		elif 'skysub' in split :
			image_list  = glob.glob(split[0]+'*'+object_name+'*skysub.fits')
		elif ('gregister' in split) & ('skysub' in split) :
			image_list  = glob.glob(split[0]+'*'+object_name+'*skysub*gre*.fits')
		elif ('gregister' not in split) & ('skysub' not in split) :
			image_list  = glob.glob(split[0]+'*'+object_name+'*0.fits')
		elif len(objectlist) == 1 :
			image_list  = objectlist[0]
		#image_list  = glob.glob('Cal*'+object_name+'*0.fits')
		k           = 0
		allfilter   = []
		for k in xrange(len(image_list)) :
			hdr     = fits.getheader(image_list[k])
			allfilter.append(hdr['filter'])
		filterset   = set(allfilter)
		infilter    = list(sorted(filterset))	
		band        = 0
		for band in xrange(len(infilter)):
			if 'gregister' in split :
				image_list_filter     = glob.glob(split[0]+'*'+object_name+'*'+infilter[band]+'*gre*.fits')
			elif 'skysub' in split :
				image_list_filter     = glob.glob(split[0]+'*'+object_name+'*'+infilter[band]+'*skysub.fits')
			elif ('gregister' in split) & ('skysub' in split) :
				image_list_filter     = glob.glob(split[0]+'*'+object_name+'*'+infilter[band]+'*skysub*gre*.fits')
			elif ('gregister' not in split) & ('skysub' not in split) :
				image_list_filter     = glob.glob(split[0]+'*'+object_name+'*'+infilter[band]+'*0.fits')
			if len(image_list_filter) > 1 : 
				images_to_align   = image_list_filter
				ref_image         = image_list_filter[1]
				outputshape       = alipy.align.shape(ref_image)
				images_to_align_1 = []
				n = 0
				for n in xrange(len(images_to_align)) :
					images        = images_to_align[n:n+1][0]
					print('\n',images,'\n')
					gregister([images])
					images_to_align_1.append(images)
			else :
				print(str(image_list_filter[0])+' is the only element. Pass.')
				pass	
	os.system('mv alipy_out/*.fits .')
	os.system('rm -r alipy_out')
	print('Done. \a')
#=================================================================
def wcsremap(imlist_name, tmpim):
	'''
	wcsremap by Andrew Becker script
	wcsremap -template template.fits -source input.fits -outIm input_remapped.fits
	'''
	import glob
	import os

	imlist  = glob.glob(imlist_name)
	imlist.sort()
	for i in xrange(len(imlist)) :
		inim = imlist[i]
		parts	= inim.split('-')
		parts[0]= 'Remap'
		outim	= '-'.join(parts)
		os.system('rm '+outim)
		com		= 'wcsremap -template '+tmpim+' -source '+inim+' -outIm '+outim
		os.system(com)
		print(outim)
	#return outim
#=================================================================
def skysex(input_image, mode='single'):
	"""
	1. Description 
	: Sky background subtraction code using SExtractor sky estimation algorithm. When using only one image, mode should be 'single'. When using more than two images, mode should be 'multi'.
	
	2. Usage
	>>> SAO_skysex() 

	3. History
	2019.12     Created by G. Lim
	2019.01.25  Added to SAO_process.py by G. Lim
	"""
	import glob
	import os,sys
	import numpy as np
	from astropy.io import ascii
	from astropy.io import fits

	pixscale = 0.311 # SAO SBIG stx16803
	config = '/data3/SAO1m/code/sex.config/'

	DETECT_MINAREA  = '5'
	DETECT_THRESH   = '5'
	ANALYSIS_THRESH = '5'
	DEBLEND_NTHRESH = '32'
	DEBLEND_MINCONT = '0.005'
	BACK_SIZE       = '128'
	BACK_FILTERSIZE = '5'
	BACKPHOTO_TYPE  = 'GLOBAL'
	BACKPHOTO_THICK = '24' # default = 24

	print('DETECT_MINAREA    : '+DETECT_MINAREA)
	print('DETECT_THRESH     : '+DETECT_THRESH)
	print('ANALYSIS_THRESH   : '+ANALYSIS_THRESH)
	print('DEBLEND_NTHRESH   : '+DEBLEND_NTHRESH)
	print('DEBLEND_MINCONT   : '+DEBLEND_MINCONT)
	print('BACK_SIZE         : '+BACK_SIZE)
	print('BACK_FILTERSIZE   : '+BACK_FILTERSIZE) 
	print('BACKPHOTO_TYPE    : '+BACKPHOTO_TYPE)
	print('BACKPHOTO_THICK   : '+BACKPHOTO_THICK)

	### Single
	if mode == 'single' :
		#input_image = input('Input image? : ')
		inim = input_image
		print(inim+' is entered.')
		infilter = inim.split('-')[-2]
		data, hdr = fits.getdata(inim, header=True)
		gain = hdr['EGAIN']
		check = 'BACKGROUND,-BACKGROUND'
		checkim = inim[:-5]+'-sky.fits,'+inim[:-5]+'-skysub.fits'
		sexcom='sex -c '+config+'skysub.sex '+inim+' -CATALOG_NAME '+inim[:-5]+'.skysub.cat -PARAMETERS_NAME '+config+'skysub.param -STARNNW_NAME '+config+'default.nnw -FILTER_NAME '+config+'default.conv -GAIN '+str(gain)+' -DETECT_THRESH '+DETECT_THRESH+' -ANALYSIS_THRESH '+ANALYSIS_THRESH+' -DETECT_MINAREA '+DETECT_MINAREA+ ' -PIXEL_SCALE ' +str(pixscale)+ ' -BACK_SIZE '+BACK_SIZE+' -BACK_FILTERSIZE '+BACK_FILTERSIZE+' -BACKPHOTO_TYPE '+BACKPHOTO_TYPE + ' -BACKPHOTO_THICK ' +BACKPHOTO_THICK+' -CHECKIMAGE_TYPE '+check+' -CHECKIMAGE_NAME '+checkim
		print(sexcom)
		os.system(sexcom)
		txt=open(inim[:-5]+'.skysub.txt','w+')
		txt.write('#image DETECT_MINAREA DETECT_THRESH ANALYSIS_THRESH DEBLEND_NTHRESH DEBLEND_MINCONT BACK_SIZE BACK_FILTERSIZE BACKPHOTO_TYPE BACKPHOTO_THICK'+'\n')
		txt.write(inim+' '+DETECT_MINAREA+' '+DETECT_THRESH+' '+ANALYSIS_THRESH+' '+DEBLEND_NTHRESH+' '+DEBLEND_MINCONT+' '+BACK_SIZE+' '+BACK_FILTERSIZE+' '+BACKPHOTO_TYPE+' '+BACKPHOTO_THICK+'\n')
		txt.close()
		print(checkim, ' are created.')

	### multi
	elif mode == 'multi' : 
		image = glob.glob('Cal*.fits')
		os.system('rm ./Cal*sky*')
		for i in range(len(image)):
			inim = image[i]
			infilter = inim.split('-')[-2]
			data, hdr = fits.getdata(inim, header=True)
			gain = hdr['EGAIN']
			check = 'BACKGROUND,-BACKGROUND'
			checkim = inim[:-5]+'-sky.fits,'+inim[:-5]+'-skysub.fits'

			sexcom='sex -c '+config+'skysub.sex '+inim+' -CATALOG_NAME '+inim[:-5]+'.skysub.cat -PARAMETERS_NAME '+config+'skysub.param -STARNNW_NAME '+config+'default.nnw -FILTER_NAME '+config+'default.conv -GAIN '+str(gain)+' -DETECT_THRESH '+DETECT_THRESH+' -ANALYSIS_THRESH '+ANALYSIS_THRESH+' -DETECT_MINAREA '+DETECT_MINAREA+ ' -PIXEL_SCALE ' +str(pixscale)+ ' -BACK_SIZE '+BACK_SIZE+' -BACK_FILTERSIZE '+BACK_FILTERSIZE+' -BACKPHOTO_TYPE '+BACKPHOTO_TYPE + ' -BACKPHOTO_THICK ' +BACKPHOTO_THICK+' -CHECKIMAGE_TYPE '+check+' -CHECKIMAGE_NAME '+checkim
			print(sexcom)
			os.system(sexcom)

			txt=open(inim[:-5]+'.skysub.txt','w+')
			txt.write('#image DETECT_MINAREA DETECT_THRESH ANALYSIS_THRESH DEBLEND_NTHRESH DEBLEND_MINCONT BACK_SIZE BACK_FILTERSIZE BACKPHOTO_TYPE BACKPHOTO_THICK'+'\n')
			txt.write(inim+' '+DETECT_MINAREA+' '+DETECT_THRESH+' '+ANALYSIS_THRESH+' '+DEBLEND_NTHRESH+' '+DEBLEND_MINCONT+' '+BACK_SIZE+' '+BACK_FILTERSIZE+' '+BACKPHOTO_TYPE+' '+BACKPHOTO_THICK+'\n')
			txt.close()
	print(checkim, ' are created.')
#=================================================================
def imcombine(imlist_name, wcsremap = False):
	"""
	1. Description
	: Stacking images in each object.
	(1) classify objects.
	(2) classify filters.
	(3) 
	2. Usage
	>>> SAO_imcombine()
	
	3. History
	2019.04.28 G. Lim added this code from separated code.
	"""
	import glob
	import os, sys
	import subprocess
	import numpy as np 
	from pyraf import iraf
	from astropy.io import fits

	lists= glob.glob(imlist_name)
	lists.sort()

	# (1) Object classification

	#lists = glob.glob('Cal*gre*.fits')
	objects = []
	for i in xrange(len(lists)) :
		inim = lists[i]
		split = inim.split('-')
		hdr = fits.getheader(inim)
		objects.append(hdr['object'])
	objectset = set(objects)
	objectlist = list(sorted(objectset))
	objectlist.sort()
	obj = 0
	for obj in xrange(len(objectlist)):
		object_name = objectlist[obj]
		if wcsremap == False :
			image_list = glob.glob('Cal*'+object_name+'*gre*.fits')
		elif wcsremap == True :
			image_list = glob.glob('Remap*'+object_name+'*0.fits')
		image_list.sort()
		k=0
		allfilter = []
		for k in xrange(len(image_list)) :
			hdr = fits.getheader(image_list[k])
			allfilter.append(hdr['filter'])
		filterset = set(allfilter)
		infilter = list(sorted(filterset))
		band = 0
		for band in xrange(len(infilter)):
			if wcsremap == False :
				list_name = 'Cal*'+object_name+'*'+infilter[band]+'*gre*.fits'
			elif wcsremap == True :
				list_name = 'Remap*'+object_name+'*'+infilter[band]+'*0.fits'
			image_list_gre = glob.glob(list_name)
			image_list_gre.sort()
			os.system('ls '+list_name+' > calibrated.list')
			hdr_gre = fits.getheader(image_list_gre[0])
			UTDATE = hdr_gre['date-obs'][0:10]
			UTSTART = hdr_gre['date-obs'][11:19]
		
			newimage = 'Calib-SAO_STX16803-'+object_name+'-'+str(UTDATE[0:4])+str(UTDATE[5:7])+str(UTDATE[8:10])+'-'+str(UTSTART[0:2])+str(UTSTART[3:5])+str(UTSTART[6:8])+'-'+infilter[band]+'-'+str(int(hdr['exptime']*len(image_list_gre)))+'.imcomb.fits'
			iraf.imcombine('@'+list_name, newimage, combine = 'median', reject='ccdclip', scale='none', zero='mode')

	print('Done.')
#=================================================================
#------------------------------------------------------------
def getimages(ra,dec,size=240,filters="grizy"):
	"""Query ps1filenames.py service to get a list of images

	ra, dec = position in degrees
	size = image size in pixels (0.25 arcsec/pixel)
	filters = string with filters to include
	Returns a table with the results
	"""	
	from astropy.table import Table
	service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
	url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
		   "&filters={filters}").format(**locals())
	table = Table.read(url, format='ascii')
	return table
#------------------------------------------------------------
def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
	"""Get URL for images in the table


	ra, dec = position in degrees
	size = extracted image size in pixels (0.25 arcsec/pixel)
	output_size = output (display) image size in pixels (default = size).
				  output_size has no effect for fits format images.

	filters = string with filters to include
	format = data format (options are "jpg", "png" or "fits")
	color = if True, creates a color image (only for jpg or png format).
			Default is return a list of URLs for single-filter grayscale images.
	Returns a string with the URL
	"""
	import numpy as np
	import requests 
	from astropy.io.votable import parse_single_table 

	if color and format == "fits":
		raise ValueError("color images are available only for jpg or png formats")
	if format not in ("jpg","png","fits"):
		raise ValueError("format must be one of jpg, png, fits")

	table	= getimages(ra, dec, size=size, filters=filters)
	url		= (	"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
				"ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
	if output_size:
		url = url + "&output_size={}".format(output_size)
	# sort filters from red to blue
	flist = ["yzirg".find(x) for x in table['filter']]
	table = table[np.argsort(flist)]
	if color:
		if len(table) > 3:
			# pick 3 filters
			table = table[[0,len(table)//2,len(table)-1]]
		for i, param in enumerate(["red","green","blue"]):
			url = url + "&{}={}".format(param,table['filename'][i])
	else:
		urlbase = url + "&red="
		url = []
		for filename in table['filename']:
			url.append(urlbase+filename)
	return url
#------------------------------------------------------------
def ps1cut(imlist_name, infilter='r'):
	"""
	1. Description 
	: Sky background subtraction code using SExtractor sky estimation algorithm. When using only one image, mode should be 'single'. When using more than two images, mode should be 'multi'.
	
	2. Usage
	>>> SAO_skysex() 

	3. History
	2019.12     Created by G. Lim
	2019.01.25  Added to SAO_process.py by G. Lim
	"""
	import glob
	import os, sys
	import numpy as np
	from astropy.io import fits

	SAO_fov    = 21.231
	pixscale   = 0.25    # PS1
	ext_factor = 1.2

	imlist     = glob.glob(imlist_name)
	imlist.sort()

	for i in xrange(len(imlist)):
		inim           = imlist[i]		
		hdr            = fits.getheader(inim)
		obj            = hdr['object'] 
		CRVAL1, CRVAL2 = hdr['CRVAL1'], hdr['CRVAL2']
		fov_pix        = ((SAO_fov*60.)/pixscale)*ext_factor
		if fov_pix > 6000 : 
			fov_pix = 6000
		param_geturl   = dict(  ra          = CRVAL1,
								dec         = CRVAL2,
								size        = int(fov_pix),
								output_size = None,
								filters     = infilter,
								format      = "fits")
		#infilter       = list(param_geturl.values())[3]
		try :
			url            = geturl(**param_geturl)
		except :
			try :
				url        = geturl(**param_geturl)
			except :
				try : 
					url        = geturl(**param_geturl)
				except :
					url        = geturl(**param_geturl)	
		save_dir   = '/data3/IMSNG/IMSNGgalaxies/refimg/SAO/'+infilter
		os.system('mkdir '+save_dir)
		fh             = fits.open(url[0])
		newname        = 'Ref-PS1-'+obj+'-'+infilter+'.fits'
		fh.writeto(save_dir+'/'+newname, overwrite=True)
		pan, panhd     = fits.getdata(save_dir+'/'+newname, header=True)
		pan0           = np.nan_to_num(pan)
		fits.writeto(save_dir+'/'+newname, pan0, panhd, overwrite=True)
	print('Please check the image at '+save_dir+'/'+newname)
#=================================================================
def gregister_inv(imlist_name, infilter='r'):
	"""
	Alipy gregister for subtraction using PS1
	from https://obswww.unige.ch/~tewes/alipy/index.html and its tutorial
	usage : run /data3/SAO1m/code/
	small pix --> large pix

	1. Description 
	: Alipy gregister for subtraction using archive data. (Ex. PS1). Reference images will be aligned to input_images. 

	from https://obswww.unige.ch/~tewes/alipy/index.html 
	
	2. Usage
	>>> SAO_gregister_inv('Cal*imcomb.fits', infilter='r') 

	3. History
	2019.12     Created by G. Lim
	2019.01.25  Added to SAO_process.py by G. Lim
	"""
	import alipy
	import glob
	import os,sys
	import numpy as np
	from astropy.io import fits

	def gregister(ref_image, images):
		id = alipy.ident.run(ref_image, images, visu=False)
		print("%20s : %20s, flux ratio %.2f" % (id[0].ukn.name, id[0].trans, id[0].medfluxratio))
		alipy.align.irafalign(id[0].ukn.filepath, id[0].uknmatchstars, id[0].refmatchstars, shape=outputshape, makepng=False)

	imlist = glob.glob(imlist_name)
	imlist.sort()
	os.system('mkdir bad_align')
	for i in xrange(len(imlist)) :
		inim = imlist[i]
		hdr  = fits.getheader(inim)
		obj = hdr['object']
		ref = 'Ref-PS1-'+obj+'-'+infilter+'.fits'
		ref_path = '/data3/IMSNG/IMSNGgalaxies/refimg/SAO/'+infilter+'/'
		images_to_align = ref_path+ref #sys.argv[1] # Public data

		images = inim[:-5]+'.ref.fits'
		print('==================================')
		print('For '+inim+'...')
		print(ref+' is copied to '+images)
		print('==================================')
		os.system('cp '+images_to_align+' '+images)
		outputshape = alipy.align.shape(inim)
		try :
			gregister(inim, [images])
		except :
			os.system('mv '+inim+' '+images+' bad_align/')
			print(inim+' is not gregistered.')
			pass
	os.system('mv alipy_out/*.fits .')
	os.system('rm -r alipy_out')
	print('Done. \a')
#=================================================================
def hotpants_public(imlist_name) :
	"""
	1. Description 
	: This code aims to convolution and subtraction of already-gregistered reference images which are downloaded from other public archives. ex) PanStarrs DR1. These public images should be gregistered in advance by running gregister code using Alipy-based Python code.
	(1) Read Calibrated images.
	(2) Read reference images of each Calibrated images.
	(3) Running hotpants code (Becker's code) setting upper and lower limit of counts, which is investigated from original public image, and force convolution on template images with -c keyword. For upper & lower count limit of PS1 image, maximum count is more than 1.e+06 and minimum count is larger than -10000. 
	
	2. Usage
	>>> SAO_hotpants_public('Calib*.fits') :

	3. History
	2018.12.12 Created by G.Lim for Maidanak image and PS1 public data.
	2018.12.14 Test for SAO image and PS1 public data (on going)
	"""
	import glob, os
	# (1)
	objlist = glob.glob(imlist_name)
	objlist.sort()

	# (2)
	ref = []
	for i in xrange(len(objlist)):
		ref_name = objlist[i][:-5]+'.ref_gregister.fits'
		ref.append(ref_name)
	ref.sort()

	# (3)
	infile=objlist
	for n in xrange(len(infile)):
		outfile='hd'+infile[n]
		convfile='hc'+infile[n]
		com = 'hotpants -c t -n i -iu 2000000 -tu 2000000 -il -10000 -tl -10000 -v 0 -inim '+infile[n]+' -tmplim '+ref[n]+' -outim '+outfile+' -oci '+convfile
		print(infile[n])
		os.system(com)
	print('All done, check it out!')
#=================================================================
def hotpants(imlist_name, refim) :
	"""
	1. Description 
	: This code aims to convolution and subtraction of already-gregistered reference images which are downloaded from other public archives. ex) PanStarrs DR1. These public images should be gregistered in advance by running gregister code using Alipy-based Python code.
	(1) Read Calibrated images.
	(2) Read reference images of each Calibrated images.
	(3) Running hotpants code (Becker's code) setting upper and lower limit of counts, which is investigated from original public image, and force convolution on template images with -c keyword. For upper & lower count limit of PS1 image, maximum count is more than 1.e+06 and minimum count is larger than -10000. 
	
	2. Usage
	>>> SAO_hotpants('Calib*.fits') :

	3. History
	2018.12.12 Created by G.Lim for Maidanak image and PS1 public data.
	2018.12.14 Test for SAO image and PS1 public data (on going)
	"""
	import glob, os
	# (1)
	objlist = glob.glob(imlist_name)
	objlist.sort()

	# (2)
	ref = refim

	# (3)
	for n in xrange(len(objlist)):
		infile = objlist[n]
		outfile='hd'+infile
		convfile='hc'+infile
		com = 'hotpants -c t -n i -iu 2000000 -tu 2000000 -il -10000 -tl -10000 -v 0 -inim '+infile+' -tmplim '+ref+' -outim '+outfile+' -oci '+convfile
		print(infile[n])
		os.system(com)
	print('All done, check it out!')
#=================================================================
