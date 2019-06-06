### run /data3/SAO1m/code/SAO_command.py
def process(today=True, subtract=False):
	"""
	1. Description
	: Perform data processing of SAO observatory.

	2. Usage
	>>> process()
	
	3. Histroy

	"""
	import glob
	import os, sys
	from datetime import datetime
	sys.path.append('/data3/SAO1m/code')
	import SAO.SAO_process as SAO

	nas_data_path = '/volume1/1mobs/STX16803/'
	qso_data_path = '/data3/SAO1m/obsdata/STX16803/'
	qso_red_path  = '/data3/SAO1m/red/STX16803/'

	if today == True :
		now = datetime.now()
		curdate = str(now.year)+str(format(now.month,'02'))+str(format((now.day-1), '02'))
	elif today == False :
		curdate = str(input('Specify date : '))

	print('=== === === === SAO pipeline starts ! === === === ===')
	#os.system('cp -r '+qso_data_path+curdate+' '+qso_red_path)
	os.chdir(qso_red_path+(curdate))
	SAO.CCDinfo(camera = 'STX16803')
	bias, dark, flat, sci = SAO.fileset(camera = 'STX16803')

	# Is there master bias?
	mzero = glob.glob('2*zero.fits')
	if len(mzero) == 0 : 
		SAO.biascom(camera='STX16803')

	# Are there master dark frames?
	mdark = glob.glob('2*dark*.fits')
	if len(mdark) == 0 :	
		SAO.darkcom(camera='STX16803')
	elif len(mdark) != 0 :
		darkframes = glob.glob('./dark/cal*dk*.fit')
		if len(darkframes) == 0:
			print('Already have master dark.')
			pass
		elif len(darkframes) != 0:
			SAO.darkcom(camera='STX16803')

	# Are there master flat in each filter?
	try :
		mskyflat = glob.glob('2*n*.sky.fits')
		if len(mskyflat) == 0 :
			SAO.flatcom(flattype='sky', process='bias', camera='STX16803')
		elif len(mskyflat) != 0 :
			skyflatframes = glob.glob('./skyflat/sky*.fit')
			if len(skyflatframes)  == 0 :
				print('Already have master skyflat.')
				pass
			elif len(skyflatframes) != 0 :
				SAO.flatcom(flattype='sky', process='bias', camera='STX16803')
	except :
		mdomeflat = glob.glob('2*n*.dome.fits')
		if len(mdomeflat) == 0 :
			SAO.flatcom(flattype='dome', process='dark', camera='STX16803')
		elif len(mdomeflat) != 0 :
			domeflatframes = glob.glob('./domeflat/dome*.fit')
			if len(domeflatframes) == 0:
				print('Already have master domeflat.')
				pass
			elif len(domeflatframes) != 0:
				SAO.flatcom(flattype='dome', process='dark', camera='STX16803')

	# Reduction 
	SAO.objpre(sci)
	# Enter WCS
	SAO.astrometry()
	# Put MJD & Check object name
	SAO.hdrcheck('a*.fit')
	# Change file name
	SAO.fnamechange()
	### Image subtraction
	if subtract == True :
		print('subtract = True.')
		SAO.gregister('Cal*.fits')
		SAO.hotpants('Cal*gre*.fits', refim)
	elif subtract == False :
		print('subtract = False')
	SAO.filemove()
	os.chdir('../')
	print 'Done.'

#process(today=True, subtract=False)
process(today=False, subtract=False)
