import psrchive,numpy as np,matplotlib.pyplot as plt,sys
from scipy.signal import savgol_filter as sv
from scipy.ndimage import gaussian_filter1d as gf

# Code is written by Devojyoti Kansabanik, NCRA.-TIFR, Pune, 23/06/2020

if len(sys.argv)<8:
	print ('Code for absolute flux calibration of  GMRT/uGMRT pulsar observation in PA mode from off-on scan on calibrator source\n\n############################################\n')
	print('Usage : python calibrate.py raw_fil caltable bad_freq gptool_bandshape fscrunch nant SB\n\n####################################\n')
	print ('raw_fil: Uncalibrated archive file\n')
	print ('caltable: g_by_tsys table\n') 
	print ('bad_freq: Bad frequency lists (for e.g.[250-270,360-380])\n')
	print ('gptool_bandhape: GPTOOL bandshape of pulsar data\n')
	print ('fscrunch: Frequency in MHz to average\n')
	print ('nant: Number of antennas used in PA mode\n')
	print ('SB: Side band used (USB or LSB)\n \n#########################################\n')
else:
	######## Inputs
	raw_fil=sys.argv[1]
	caltable=sys.argv[2]
	bad_freq=sys.argv[3]
	gptool_bandshape=sys.argv[4]
	bscrunch=float(sys.argv[5])
	nant=int(sys.argv[6])
	sb=sys.argv[7]

	bad_freq=bad_freq.split(',')
	bandshape=np.loadtxt(gptool_bandshape,unpack=True,usecols=(1))
	cal=np.load(caltable)*nant
	if sb=='USB':		# If LSB is used gptool bandshape have to be fliped
		cal=np.flip(cal)
		bandshape=np.flip(bandshape)
	elif sb=='LSB':
		cal=np.flip(cal)
	else:
		print ('Choose proper sideband....\n')
		sys.exit()
	print ('Loading raw data..................\n')
	raw=psrchive.Archive_load(raw_fil)
	nchan=raw.get_nchan()
	nsub=raw.get_nsubint()
	frequencies=raw.get_frequencies()
	start_freq=frequencies[0]
	bw=raw.get_bandwidth()
	chanres=bw/nchan
	bad_chan=[]
	nbin=raw.get_nbin()
	for i in bad_freq:			# Making bad channels list
		a=float(i.split('-')[0])
		b=float(i.split('-')[1])
		c=int((a-start_freq)/chanres)
		d=int((b-start_freq)/chanres)
		if d>c:
			for i in range(c,d):
				bad_chan.append(i)
		elif c>d:
			for i in range(d,c):
				bad_chan.append(i)
		else:
			bad_chan.append(c)

	raw.remove_baseline()

	weight=[]
	print ('Calculating weights...............\n')
	for ichan in range(nchan):		# Caliculatinh weight, i.e. unflagged fraction in each channel
		if ichan in bad_chan:
			weight.append(0)
		else:
			count=0
			for isub in range(nsub):
				subint=raw[isub]
				pdata=subint.get_Profile(0,ichan).get_amps()
				if np.sum(np.count_nonzero(pdata))!=0:
					count+=1
			frac=count/nsub
			weight.append(frac)

	print ('Set weight to data........\n')
	for isub in range(nsub):		# Put weight in data
		subint=raw[isub]
		for ichan in range(nchan):
			if ichan in bad_chan:
				subint.set_weight(ichan,0)
				subint.get_Profile(0,ichan).scale(0)
			else:
				subint.set_weight(ichan,weight[ichan])
				
	#Un-normalise the gptool normalisation
	for ichan in range(nchan):
		for isub in range(nsub):
			subint=raw[isub]
			prof=subint.get_Profile(0,ichan)
			if ichan in bad_chan:
				prof.scale(0)
			else:
				prof.scale(bandshape[ichan])

	if not raw.get_dedispersed():
		raw.dedisperse()
	raw.centre_max_bin()
	data=raw.get_data()
	data=np.nanmean(data,axis=0)
	data=data.reshape(data.shape[-2],data.shape[-1])
	folded_data=np.mean(data,axis=0)
	folded_data=folded_data/max(folded_data)
	plt.plot(folded_data)
	plt.show()
	onleft=int(input('On pulse start window:'))
	onright=int(input('Off pulse end window:'))
	print ('On pulse window:',onleft,'-',onright,'.........\n')
	offscale=[]
	print ('Calibrating......\n') # Calculating the Jy scaling factor
	for ichan in range(nchan):
		g_by_tsys=cal[ichan]
		w=weight[ichan]
		if w!=0 and g_by_tsys!=0:
			prof=data[ichan,:]
			offpulse=np.append(prof[:onleft],prof[onright:])
			sigma=(np.std(offpulse))		#Off pulse standard deviation
			del_s=1/(g_by_tsys*np.sqrt(2*w*(nsub/len(prof))*abs(chanres)*(10**6)))
			offscale.append(del_s/(sigma))				# Required Jy scaling 
		else:
			offscale.append(0)
	nbin=raw.get_nbin()
	for isub in range(nsub):						# Scale the data in Jy
		subint=raw[isub]
		for ichan in range(nchan):
			prof=subint.get_Profile(0,ichan)
			pdata=prof.get_amps()
			c=offscale[ichan]
			if np.sum(np.count_nonzero(pdata))!=0:
				prof.scale(float(c))
			else:
				prof.scale(0)
		
	print ('Unload calibrated archive.....\n')		# Saving calibrated data
	raw.unload(raw_fil[:-3]+'_fluxcalib.ar')
	if not raw.get_dedispersed():
		raw.dedisperse()
	raw.tscrunch()
	raw.fscrunch(int(bscrunch/abs(chanres)))
	raw.centre_max_bin()
	nchan=raw.get_nchan()
	freq=raw.get_frequencies()
	freqlist=[]
	peak_flux=[]
	mean_flux=[]
	fluxerr=[]
	subint=raw[0]
	print ('Frequency crunching and profile saving...\n') # Average calibrated archive making
	for ichan in range(nchan):
		prof=subint.get_Profile(0,ichan)
		pdata=prof.get_amps()
		if np.sum(np.count_nonzero(pdata))>=0:
			freqlist.append(freq[ichan])
			off_pulse=np.append(pdata[:onleft],pdata[onright:])
			mean=np.mean(off_pulse)
			pdata-=mean
			off_pulse=np.append(pdata[:onleft],pdata[onright:])
			std=np.std(off_pulse)
			if std>0:
				snr=np.max(pdata)/std
				if snr>4:
					peak_flux.append(np.max(pdata))
					mean_flux.append(np.nanmean(pdata))
					fluxerr.append(std)
					np.save(raw_fil[:-3]+'_'+str(int(freq[ichan])),pdata)
				else:
					prof.scale(0)
			else:
				prof.scale(0)
		else:
			prof.scale(0)
	print ('Unloading frequency crunched archive....\n')
	raw.unload(raw_fil[:-3]+'_fluxcalib.fTp')
	np.save(raw_fil[:-3]+'_calibrated_flux',np.array([mean_flux,peak_flux,fluxerr]))   # Save the averaged mean flux density calculated from Sum of bins
	print ('Calibration finished for : '+raw_fil+'\n#################################################')

















































