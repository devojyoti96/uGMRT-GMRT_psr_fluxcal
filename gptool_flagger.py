import numpy as np,matplotlib.pyplot as plt,sys,os,glob

# Code is written by Devojyoti Kansabanik, NCRA-TIFR, Pune, 23/06/2020
if len(sys.argv)<9:
	print ('##############################\nCode to do RFI flagging on calibrators and create bandshapes. Note that gptool.in file must be in the same directory where the code is and all parameters were put in gptool.in correctly.\n#############################\n')
	print ('Usage : python gptool_flagger.py datadir fluxcal_on fluxcal_off1 fluxcal_off2 phasecal_on phasecal_off1 phasecal_off2 outdir\n############################################')
	print ('datadir : Directory where raw data is stored\n')
	print ('fluxcal_on: Raw fluxcal on scan file\n')
	print ('fluxcal_off1: Raw fluxcal off first scan\n')
	print ('fluxcal_off2: Raw fluxcal off second scan (If not available put the same name as fluxcal_off1)\n')
	print ('fluxcal_on: Raw phasecal on scan file\n')
	print ('fluxcal_off1: Raw phasecal off first scan\n')
	print ('fluxcal_off2: Raw phasecal off second scan (If not available put the same name as phasecal_off1)\n') 
	print ('outdir : Output directory to store bandshapes\n##########################################\n')
else:
	datadir=sys.argv[1]

	fluxcal_on=sys.argv[2]
	fluxcal_off1=sys.argv[3]
	fluxcal_off2=sys.argv[4]

	phasecal_on=sys.argv[5]
	phasecal_off1=sys.argv[6]
	phasecal_off2=sys.argv[7]

	outdir=sys.argv[8]

	print ('RFI flagging for fluxcal on scan.....................\n##############################\n')
	os.system('gptool -f '+datadir+'/'+fluxcal_on+' -o '+datadir+' -m 64 -nodedisp')
	os.system('mv bandshape.gpt '+fluxcal_on+'.bandshape')

	print ('RFI flagging for fluxcal off scan ...................\n###############################\n')
	os.system('gptool -f '+datadir+'/'+fluxcal_off1+' -o '+datadir+' -m 64 -nodedisp')
	os.system('mv bandshape.gpt '+fluxcal_off1+'.bandshape')

	if fluxcal_off1!=fluxcal_off2:
		os.system('gptool -f '+datadir+'/'+fluxcal_off2+' -o '+datadir+' -m 64 -nodedisp')
		os.system('mv bandshape.gpt '+fluxcal_off2+'.bandshape')


	print ('RFI flagging for phasecal on scan...........................\n##############################\n')
	os.system('gptool -f '+datadir+'/'+phasecal_on+' -o '+datadir+' -m 64 -nodedisp')
	os.system('mv bandshape.gpt '+phasecal_on+'.bandshape')

	print ('RFI flagging for phasecal off scan...........................\n###############################\n')
	os.system('gptool -f '+datadir+'/'+phasecal_off1+' -o '+datadir+' -m 64 -nodedisp')
	os.system('mv bandshape.gpt '+phasecal_off1+'.bandshape')

	if phasecal_off1!=phasecal_off2:
		os.system('gptool -f '+datadir+'/'+phasecal_off2+' -o '+datadir+' -m 64 -nodedisp')
		os.system('mv bandshape.gpt '+phasecal_off2+'.bandshape')

	try:
		os.system('mkdir '+outdir)
	except OSError:
		pass
	os.system('mv *.bandshape '+outdir)
	os.system('rm -rf *.gpt *.inf')
	print ('\n###############################\n Bandshapes saved at : '+outdir+'\n##########################################\n')






















