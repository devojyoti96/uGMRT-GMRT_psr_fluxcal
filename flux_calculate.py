import os,glob,sys,psrchive,numpy as np,matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Code is written by Devojyoti Kansabanik, NCRA-TIFR, Pune, 23/06/2020

def std_model(p,scale):
	return p*scale

def chisq_func(data,model,sigma):
	chi=np.sum(((data-model)/sigma)**2)
	nu=len(data)-1
	red_chi=chi/nu
	return chi,red_chi

if len(sys.argv)<9:
	print ('Code to calculate mean flux density of pulsar using template fitting method. Before running the code save the time crunched and frequency slices pulse profiles in a directory. Also save standard profiles in a separate directory. Make sure there are no other files in those directories\n\n########################################\n')
	print ('Useage : python flux_calculate.py prof_dir_path std_prof_dir_path fluxtable minsnr nbin avg_prof\n\n#########################################\n')
	print ('prof_dir_path : Give full path of profile directory\n')
	print ('std_prof_dir_path : Give full path of standard profile directory\n')
	print ('fluxtable : Name of output fluxtable to save result\n')
	print('minsnr : Minimum signal-to-noise of profile to save in fluxtable\n')
	print('bin : Number of phase bin average\n')
	print('avg_prof : Frequency crunched average profile\n')
	print('bad_freq : List of bad frequencies \n')
	print ('fluxmode : Method to calculate flux density. Three methods are available: \'peak\',\'mean_width\',\'mean_prof\'\n####################################################\n')
else:
	path=sys.argv[1]
	std_path=sys.argv[2]
	fluxtable=sys.argv[3]
	minsnr=float(sys.argv[4])
	nbin=int(sys.argv[5])
	avg_prof=sys.argv[6]
	pwd=os.getcwd()
	os.chdir(path)
	bad_freq=sys.argv[7]
	fluxmode=sys.argv[8]

	bad_freq=bad_freq.split(',')
	bad_freq_list=[]
	for i in bad_freq:
		a=int(i.split('-')[0])
		b=int(i.split('-')[-1])
		bad_freq_list.append(np.array([x for x in range(a,b)]))
	bad_freq_list=np.array(bad_freq_list)
	bad_freq=np.append(bad_freq_list[0],bad_freq_list[1]).tolist()
	file_list=glob.glob('*')
	if len(file_list)==0:
		print ('No files found, check the path\n')
		sys.exit()
	freq_list=[]
	flux_list=[]
	fluxerr_list=[]
	red_chi_list=[]
	print (avg_prof)
	os.system('pam -FTp '+avg_prof+' -m')
	do_fit=input('Want to do fitting for average profile? y/n: ')
	while do_fit=='y':
		centering=input('Want center the profile? y/n: ')
		if centering=='n':
			os.system('paas -i -s '+avg_prof+'.std'+' '+avg_prof)
		else:
			os.system('paas -i -C -s '+avg_prof+'.std'+' '+avg_prof)
		repeat=input('Want to repat the fiiting? y/n: ')
		if repeat=='n':
			break
	if do_fit=='y':
		want_save=input('Want to keep this standard profile? y/n: ')
		if want_save=='y':
			print ('Standard profile unloaded : '+avg_prof+'.std'+'\n###########################################\n')
		elif want_save=='n':
			os.system('rm -rf '+avg_prof+'.std')
	os.system('rm -rf paas*')
	std_avg_prof=psrchive.Archive_load(avg_prof+'.std')
	pdata=std_avg_prof[0].get_Profile(0,0).get_amps()
	w_avg=int(np.sum(pdata)/max(pdata))
	print ('Effective average pulse width W: '+str(w_avg)+'\n############################\n')
	plt.plot(pdata)
	plt.show()
	onleft=int(input('On pulse start bin:'))
	onright=int(input('On pulse end bin:'))
	for fil in file_list:
		raw=psrchive.Archive_load(fil)
		freq=(int(raw.get_frequencies()[0]))
		if int(freq) in bad_freq:
			continue
		else:
			std_file_name=std_path+'/freq_'+str(freq)+'.std'
			if os.path.isfile(std_file_name)==False:
				print ('Standard profiles are not found at ',std_file_name,' check the path...\n')
				continue
			try:
				std_prof=psrchive.Archive_load(std_file_name)
				if nbin>=2:
					raw.bscrunch(nbin)
					std_prof.bscrunch(nbin)
			
				pdata=raw[0].get_Profile(0,0).get_amps()
				pstd=std_prof[0].get_Profile(0,0).get_amps()
				pstd_norm=pstd/max(pstd)
				max_pos_pdata=np.where(pdata==max(pdata))[0][0]
				max_pos_std=np.where(pstd_norm==max(pstd_norm))[0][0]
				if abs(max_pos_std-max_pos_pdata)>5:
					pdata=np.roll(pdata,max_pos_std-max_pos_pdata)
				off_pulse=np.append(pdata[:onleft],pdata[onright:])
				sigma=np.ones(len(pdata))*np.std(off_pulse)
				popt,pcov=curve_fit(std_model,pstd_norm,pdata,sigma=sigma,p0=[max(pdata)],absolute_sigma=True)
				chisq,redchi=chisq_func(pdata,std_model(pstd_norm,popt[0]),sigma)
				fitted_profile=std_model(pstd_norm,popt[0])
				snr=popt[0]/np.std(off_pulse)
				weff=int(np.sum(fitted_profile)/max(fitted_profile))
				if fluxmode=='peak':
					flux=popt[0]
					fluxerr=np.std(off_pulse)
				elif fluxmode=='mean_width':
					flux=(popt[0]*w_avg)/len(pdata)
					flux_err=np.std(off_pulse)/(np.sqrt(len(pdata)))
				elif fluxmode=='mean_prof':
					flux=(popt[0]*weff)/len(pdata)
					flux_err=np.std(off_pulse)/(np.sqrt(len(pdata)))
				if snr>minsnr:
					if fluxmode=='peak':
						print ('\n##########################\nFlux density at '+str(int(freq))+' MHz : '+str(flux*1000)+','+str(flux_err*1000)+' mJy\n######################################\n\n')
					elif fluxmode=='mean_width' or fluxmode=='mean_prof':
						print ('\n##########################\nFlux density at '+str(int(freq))+' MHz : '+str(flux*1000)+','+str(flux_err*1000)+' mJy , width : '+ste(weff)+'\n######################################\n\n')
					else:
						print ('Choose proper flux calculation method...\n')
						sys.exit()
					freq_list.append(freq)
					flux_list.append(flux)
					fluxerr_list.append(flux_err)
					red_chi_list.append(redchi)
					wlist.append(weff)
			except RuntimeError:
				pass
	os.chdir(pwd)
	plt.errorbar(freq_list,[a*1000 for a in flux_list],[b*1000 for b in fluxerr_list],fmt='bo',markersize=4,mfc='white',capsize=4,zorder=1)
	if fluxmode='peak':
		plt.ylabel('Peak flux density (mJy)')
	elif fluxmode=='mean_width' or fluxmode=='mean_prof':
		plt.ylabel('Mean flux density (mJy)')
	plt.xlabel('Frequency (MHz)')
	plt.savefig(fluxtable.split('.')[0]+'.png')
	plt.show()
	result=np.array([freq_list,flux_list,fluxerr_list,red_chi_list])
	np.save(fluxtable,result)
	d=[(i/len(final_data))*100 for i in w_list]
	from scipy.signal import savgol_filter as sv
	if len(d)<=5:
		plt.plot(freq_list,d)
	else:
	plt.plot(freq_list,sv(d,5,1))
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Duty cycle (%)')
	plt.savefig(fluxtable+'_width.png')
	plt.show()










