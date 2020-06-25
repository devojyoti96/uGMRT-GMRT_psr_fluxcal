import numpy as np,matplotlib.pyplot as plt,sys,psrchive
from scipy.ndimage import gaussian_filter1d as gf
from scipy.signal import medfilt,savgol_filter as sv
from scipy import interpolate 
from scipy.ndimage import gaussian_filter1d as gf
from scipy.optimize import curve_fit

# Code is written by Devojyoti Kansabanik, NCRA-TIFR, Pune, 23/06/2020

def fill_nan(A): 
	inds = np.arange(A.shape[0]) 
	good = np.where(np.isfinite(A)) 
	f = interpolate.interp1d(inds[good], A[good],bounds_error=False) 
	B = np.where(np.isfinite(A),A,f(inds)) 
	return B 

def flux_cal_spectra(f,a0,a1,a2,a3): #Function calculating spectra of calibrator, a are the coefficients in polynominal power law from VLA calibrator manual and f is frequency.
        logf=a0+a1*np.log10(f)+a2*(np.log10(f))**2+a3*(np.log10(f))**3
        return 10**logf

def bad_chan_func(chan1_freq,chanres):
	bad_chan=[]
	for i in bad_freq:
		a=float(i.split('-')[0])
		b=float(i.split('-')[1])
		c=int((a-chan1_freq)/chanres)
		d=int((b-chan1_freq)/chanres)
		if d>c:
			for i in range(c,d):
				bad_chan.append(i)
		elif c>d:
			for i in range(d,c):
				bad_chan.append(i)
		else:
			bad_chan.append(c)
	return bad_chan

def resolution_match(x,outchanres):
	if len(x)<outchanres:
		print ('Number of output channel can not be less than')
	avg_factor=int(len(x)/outchanres)
	y=np.mean(x.reshape(-1,avg_factor),axis=1)
	return y

if len(sys.argv)<18:
	print ('\n#######################################\nCode to calculate G/T_sys value using flux calibrator and phase calibrator observation\n#########################################\n')
	print ('Usage : python bandshape_cal.py fluxcal_on fluxcal_off1 fluxcal_off2 phasecal_on phasecal_off1 phasecal_off2 pulsar_off start_freq bandwidth sideband caltable bad_freq nant out_chanres smooth_window cal_beam psr_beam\n######################################')
	print ('fluxcal_on : Fluxcal on scan gptool bandshape\n')
	print ('fluxcal_off1 : Fluxcal off first scan gptool bandshape\n')
	print ('fluxcal_off2 : Fluxcal off second scan gptool bandshape (If not available put the first off scan)\n')
	print ('phasecal_on : Phasecal on scan gptool bandshape\n')
	print ('phasecal_off1 : Phasecal off first scan gptool bandshape\n')
	print ('phasecal_off2 : Phasecal off second scan gptool bandshape\n')
	print ('pulsar_off : Pulsar data gptool bandshape\n')
	print ('start_freq : Lowest frequency of the band in MHz (for e.g. 250-450 MHz it is 250)\n')
	print ('bandwidth : Bnadwidth in MHz\n')
	print ('sideband : USB or LSB\n')
	print ('caltable: Output caltable name\n')
	print ('bad_freq : Bad frequency list in MHz (for e.g. \'250-270,360-380\')\n')
	print ('nant : Number of antenna used in PA\n')
	print ('out_chanres : Output number of channels\n')
	print ('smooth_window : NUmber of channels for smoothing\n')
	print ('cal_beam : Calibrator beam mode (\'PA\' or \'CDP\')\n')
	print ('psr_beam : Pulsar beam mode (\'PA\' or \'CDP\')\n\n##################################################################\n')
else:
	# Inputs
	on_fluxcal=sys.argv[1]
	off1_fluxcal=sys.argv[2]
	off2_fluxcal=sys.argv[3]
	on_phasecal=sys.argv[4]
	off1_phasecal=sys.argv[5]
	off2_phasecal=sys.argv[6]
	psr_bandshape=sys.argv[7]
	chan1_freq=float(sys.argv[8])
	bw=float(sys.argv[9])
	sideband=sys.argv[10]
	caltable=sys.argv[11]
	bad_freq=sys.argv[12]
	nant=int(sys.argv[13])
	outchanres=float(sys.argv[14])
	smooth_window=int(sys.argv[15])
	cal_mode=sys.argv[16]
	psr_mode=sys.argv[17]	
	
	if smooth_window%2==0:
		smooth_window+=1

	bad_freq=bad_freq.split(',')
	y_on_fluxcal,y_on_fluxcal_rms=np.loadtxt(on_fluxcal,unpack=True,usecols=(1,3))
	y_off1_fluxcal,y_off1_fluxcal_rms=np.loadtxt(off1_fluxcal,unpack=True,usecols=(1,3))
	y_off2_fluxcal,y_off2_fluxcal_rms=np.loadtxt(off2_fluxcal,unpack=True,usecols=(1,3))
	y_off_fluxcal=np.median(np.concatenate((y_off1_fluxcal.reshape(1,-1),y_off2_fluxcal.reshape(1,-1)),axis=0),axis=0)
	y_off_fluxcal_rms=np.median(np.concatenate((y_off1_fluxcal_rms.reshape(1,-1),y_off2_fluxcal_rms.reshape(1,-1)),axis=0),axis=0)

	y_on_phasecal,y_on_phasecal_rms=np.loadtxt(on_phasecal,unpack=True,usecols=(1,3))
	y_off1_phasecal,y_off1_phasecal_rms=np.loadtxt(off1_phasecal,unpack=True,usecols=(1,3))
	y_off2_phasecal,y_off2_phasecal_rms=np.loadtxt(off2_phasecal,unpack=True,usecols=(1,3))
	y_off_phasecal=np.median(np.concatenate((y_off1_phasecal.reshape(1,-1),y_off2_phasecal.reshape(1,-1)),axis=0),axis=0)
	y_off_phasecal_rms=np.median(np.concatenate((y_off1_phasecal_rms.reshape(1,-1),y_off2_phasecal_rms.reshape(1,-1)),axis=0),axis=0)
	psr_off,psr_off_rms=np.loadtxt(psr_bandshape,unpack=True,usecols=(1,3))

	if sideband=='LSB':
		y_on_fluxcal=np.flip(y_on_fluxcal)
		y_on_fluxcal_rms=np.flip(y_on_fluxcal_rms)
		y_off_fluxcal=np.flip(y_off_fluxcal)
		y_off_fluxcal_rms=np.flip(y_off_fluxcal_rms)
		y_on_phasecal=np.flip(y_on_phasecal)
		y_on_phasecal_rms=np.flip(y_on_phasecal_rms)
		y_off_phasecal=np.flip(y_off_phasecal)
		y_off_phasecal_rms=np.flip(y_off_phasecal_rms)
		psr_off=np.flip(psr_off)
		psr_off_rms=np.flip(psr_off_rms)
	if cal_mode=='CDP':
		y_on_fluxcal*=np.sqrt(2)
		y_on_fluxcal_rms*=np.sqrt(2)
		y_off_fluxcal*=np.sqrt(2)
		y_off_fluxcal_rms*=np.sqrt(2)
		y_on_phasecal*=np.sqrt(2)
		y_on_phasecal_rms*=np.sqrt(2)
		y_off_phasecal*=np.sqrt(2)
		y_off_phasecal_rms*=np.sqrt(2)
	if psr_mode=='CDP':
		psr_off*=np.sqrt(2)
		psr_off_rms*=np.sqrt(2)
		
	if len(y_on_fluxcal)!=outchanres:
		if len(y_on_fluxcal)<outchanres:
			print ('Input channels can not be smaller than output number of channels\n')
		else:
			y_on_fluxcal=resolution_match(y_on_fluxcal,outchanres)	
	if len(y_off_fluxcal)!=outchanres:
                if len(y_off_fluxcal)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        y_off_fluxcal=resolution_match(y_off_fluxcal,outchanres)
	if len(y_off_fluxcal)!=outchanres:
                if len(y_off_fluxcal)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        y_off_fluxcal=resolution_match(y_off_fluxcal,outchanres)
	if len(y_on_phasecal)!=outchanres:
                if len(y_on_phasecal)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        y_on_phasecal=resolution_match(y_on_phasecal,outchanres)
	if len(y_on_phasecal)!=outchanres:
                if len(y_on_phasecal)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        y_on_phasecal=resolution_match(y_on_phasecal,outchanres)
	if len(y_off_phasecal)!=outchanres:
                if len(y_off_phasecal)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        y_off_phasecal=resolution_match(y_off_phasecal,outchanres)
	if len(psr_off)!=outchanres:
                if len(psr_off)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        psr_off=resolution_match(psr_off,outchanres)
	if len(y_on_phasecal_rms)!=outchanres:
                if len(y_on_phasecal_rms)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        y_on_phasecal_rms=resolution_match(y_on_phasecal_rms,outchanres)
	if len(y_off_phasecal_rms)!=outchanres:
                if len(y_off_phasecal_rms)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        y_off_phasecal_rms=resolution_match(y_off_phasecal_rms,outchanres)
	if len(psr_off_rms)!=outchanres:
                if len(psr_off_rms)<outchanres:
                        print ('Input channels can not be smaller than output number of channels\n')
                else:
                        psr_off_rms=resolution_match(psr_off_rms,outchanres)

	phasecal_off_snr=sv(((y_off_phasecal)/y_off_phasecal_rms),smooth_window,3)
	bad_chan=bad_chan_func(chan1_freq,bw/len(phasecal_off_snr))
	for ichan in bad_chan:
		phasecal_off_snr[ichan]=0
	phasecal_off_snr_norm=phasecal_off_snr/max(phasecal_off_snr)
	phasecal_on_snr=sv(((y_on_phasecal)/y_on_phasecal_rms),smooth_window,3)
	bad_chan=bad_chan_func(chan1_freq,bw/len(phasecal_on_snr))
	for ichan in bad_chan:
	        phasecal_on_snr[ichan]=0
	phasecal_on_snr_norm=phasecal_on_snr/max(phasecal_on_snr)
	freq=np.linspace(chan1_freq,chan1_freq+bw,len(y_on_fluxcal))
	psr_snr=psr_off/psr_off_rms
	bad_chan=bad_chan_func(chan1_freq,bw/len(psr_snr))
	for ichan in bad_chan:
		psr_snr[ichan]=0
	psr_snr_norm=psr_snr/max(psr_snr)

	fluxcal_spectra=flux_cal_spectra(freq/1000,1.2481,-0.4507,-0.1798,0.0357)	#3C286
	count_fluxcal=(y_on_fluxcal-y_off_fluxcal)/(fluxcal_spectra)	# counts/Jy for fluxcal
	bad_chan=bad_chan_func(chan1_freq,bw/len(count_fluxcal))
	for ichan in bad_chan:
		count_fluxcal[ichan]=np.nan
	f1=(y_on_phasecal-y_off_phasecal)/count_fluxcal
	f1[f1<0]=np.nan
	valid = ~(np.isnan(f1))
	f1=gf(f1,20)
	try:
		popt,pcov=curve_fit(flux_cal_spectra,freq[valid],f1[valid])
		f1_fitted=flux_cal_spectra(freq,popt[0],popt[1],popt[2],popt[3])
		plt.plot(freq,f1_fitted)
		g_by_tsys=((y_on_phasecal*phasecal_on_snr_norm)-(y_off_phasecal*phasecal_off_snr_norm))/(psr_off*f1_fitted)
	except:
		f1=(flux_cal_spectra(freq/1000,0.98941947, -0.70830277, -0.15327315,  0.14418446))
		g_by_tsys=((y_on_phasecal*phasecal_on_snr_norm)-(y_off_phasecal*phasecal_off_snr_norm))/(psr_off*f1)
	for ichan in bad_chan:
		g_by_tsys[ichan]=np.nan
	g_by_tsys=fill_nan(g_by_tsys)
	g_by_tsys=np.nan_to_num(g_by_tsys/nant)
	g_by_tsys[g_by_tsys<0]=0
	g_by_tsys=medfilt(g_by_tsys,smooth_window)
	g_by_tsys=sv(g_by_tsys,smooth_window,2)
	plt.plot(freq,g_by_tsys*1000,label=str(int(min(freq)))+'-'+str(int(max(freq)))+'MHz')
	plt.xlabel('Frequency (MHz)')
	plt.ylabel(r'$\frac{G}{T_{sys}}\times1000 (Jy^{-1})$')
	pos=np.where(g_by_tsys>0)
	x1=int(pos[0][0])
	x2=int(pos[0][-1])
	plt.xlim(freq[x1],freq[x2])
	plt.ylim(bottom=0)
	np.save(caltable,g_by_tsys)
	plt.savefig(caltable+'.png')
	plt.show()
