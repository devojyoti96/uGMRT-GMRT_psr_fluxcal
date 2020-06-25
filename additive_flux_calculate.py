import psrchive,numpy as np,matplotlib.pyplot as plt,sys,os,glob
from scipy.optimize import curve_fit
# Code is written by Devojyoti Kansabanik , NCRA-TIFR, Pune, 23/06/2020

def model(pstd_norm,scale):
	return pstd_norm*scale

def chifun(data,model,sigma):
	chi=np.sum(((data-model)/sigma)**2)
	redchi=chi/(len(data)-1)
	return chi,redchi

if len(sys.argv)<11:
	print ('\n#########################################################\n\nCode to calculate mean flux density by template fitting with interactively adding frequency slices uniformly/non-uniformly based on user requirement.\nOutput will be a numpy table and spectra plot.\nInitial profiles and standard profiles must have same bandwidth\n##############################################\n\n')
	print ('Usage : python additive_flux_calculate.py prof_path std_prof_path fluxtable minsnr nbin prof_out_dir avg_prof bad_freq flux_mode\n#############################################\n')
	print ('prof_path : Full path of input profile directory\n')
	print ('std_prof_path : Full path of input standard profile directory\n')
	print ('fluxtable : Output fluxtable name\n')
	print ('minsnr : Minimum signal-to-noise of the profile to save in result\n')
	print ('nbin : Number of phase bin to average\n')
	print ('prof_out_dir : Output directory to save average profiles\n')
	print ('avg_prof : Frequency scrunched profile \n')
	print ('bad_freq : Bad frequency ranges (for e.g. 250-170,360-380)\n')
	print ('bypass : Bypass bad frequencies\n')
	print ('flux_mode: Flux density calculation method. Three modes are available, \'peak\',\'mean_width\',\mean_prof\'\n\n#################################################################\n')
else:
	prof_path=sys.argv[1]
	std_prof_path=sys.argv[2]
	fluxtable=sys.argv[3]
	minsnr=float(sys.argv[4])
	nbin=int(sys.argv[5])
	prof_out_dir=sys.argv[6]
	avg_prof=sys.argv[7]
	bad_freq=sys.argv[8]
	bypass=sys.argv[9]
	fluxmode=sys.argv[10]

	bad_freq=bad_freq.split(',')
	bad_freq_list=[]
	for i in bad_freq:
		a=int(i.split('-')[0])
		b=int(i.split('-')[-1])
		bad_freq_list.append(np.array([x for x in range(a,b)]))
	bad_freq_list=np.array(bad_freq_list)
	bad_freq=np.append(bad_freq_list[0],bad_freq_list[1]).tolist()
	profiles=glob.glob(prof_path+'/*')
	std_profiles=glob.glob(std_prof_path+'/*')
	std_profiles=sorted(std_profiles)
	std_freques=np.array([int(i.split('/')[-1].split('.')[0].split('_')[-1]) for i in std_profiles])
	data=[]
	sdata=[]
	flist=[]
	flux_list=[]
	freq_list=[]
	flux_err_list=[]
	chi_list=[]
	red_chi_list=[]
	c=0
	os.system('pam -FTp '+avg_prof+' -m')
	do_fit=input('Want to do fitting for average profile? y/n: ')
	while do_fit=='y':
		centering=input('Want center the profile? y/n: ')
		if centering=='n':
			os.system('paas -i -s '+avg_prof+'.std'+' '+avg_prof+' -d /xs')
		else:
			os.system('paas -i -C -s '+avg_prof+'.std'+' '+avg_prof+' -d /xs')
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
	print (bypass)
	print ('Effective average pulse width W: '+str(w_avg)+'\n############################\n')
	plt.plot(pdata)
	plt.show()
	w_list=[]
	onleft=int(input('On pulse start bin:'))
	onright=int(input('On pulse end bin:'))
	for i in range(len(profiles)):
		prof=profiles[i]
		raw=psrchive.Archive_load(prof)
		freq=int(raw.get_frequencies()[0])
		if freq in bad_freq and bypass=='True':
			continue	
		else:
			index= (np.argmin(abs(std_freques-freq)))
			std_file_name=std_profiles[index]
			std=psrchive.Archive_load(std_file_name)
			if nbin>=2:
				raw.bscrunch(nbin)
				std.bscrunch(nbin)
			raw.tscrunch()
			raw.fscrunch()
			pdata=raw[0].get_Profile(0,0).get_amps()
			std.tscrunch()
			std.fscrunch()
			std_data=std[0].get_Profile(0,0).get_amps()
			want_to_add=input('\nWant add '+str(int(freq))+' MHz frequency slice? y/n :')
			if want_to_add=='y':
				data.append(pdata.tolist())
				sdata.append(std_data.tolist())
				flist.append(int(freq))
				c+=1
			else:
				if c==0:
					final_data=pdata
					final_sdata=std_data
					flist=[int(freq)]		
				else:
					data.append(pdata.tolist())
					sdata.append(std_data.tolist())
					final_data=np.array(data)
					final_sdata=np.array(sdata)
					final_data=np.mean(final_data,axis=0)
					final_sdata=np.mean(final_sdata,axis=0)
					flist.append(int(freq))
					c=0
				max_pos_prof=np.argmax(final_data)
				max_pos_std=np.argmax(final_sdata)
				if abs(max_pos_std-max_pos_prof)>5: 
					final_data=np.roll(final_data,max_pos_std-max_pos_prof)
				if onleft<onright:
					off_pulse=np.append(final_data[:onleft],final_data[onright:])
				else:
					off_pulse=final_data[onright:onleft]
				off_sigma=np.std(off_pulse)
				sigma=np.ones(len(final_data))*off_sigma
				final_sdata_norm=final_sdata/max(final_sdata)
				popt,pcov=curve_fit(model,final_sdata_norm,final_data,sigma=sigma,absolute_sigma=True)
				chi,redchi=chifun(final_data,model(final_sdata_norm,popt[0]),sigma)
				fitted_profile=model(final_sdata_norm,popt[0])
				weff=int(np.sum(fitted_profile)/max(fitted_profile))
				if fluxmode=='peak':
					flux=popt[0]
					fluxerr=np.sqrt(np.std(off_pulse)**2+(0.1*popt[0])**2)
				elif fluxmode=='mean_width':
					flux=popt[0]*w_avg/len(final_data)
					fluxerr=(np.sqrt(np.std(off_pulse)**2+(0.1*flux)**2))/np.sqrt(len(final_data))
				elif fluxmode=='mean_prof':
					flux=popt[0]*weff/len(final_data)
					fluxerr=(np.sqrt(np.std(off_pulse)**2+(0.1*flux)**2))/np.sqrt(len(final_data))
				else:
					print ('Choose correct flux density calculation method....\n')
					sys.exit()
				print ('\n#####################\nFrequencies added :',flist,'\n############################################\n')
				snr=float(popt[0]/np.std(off_pulse))
				if snr>minsnr and np.std(off_pulse)!=0 and snr!=np.inf:
					if os.path.isdir(prof_out_dir)==False:
						os.system('mkdir '+prof_out_dir)
					np.save(prof_out_dir+'/freq_'+str(np.mean(np.array(flist))),final_data)
					if fluxmode=='peak' or fluxmode=='mean_width':
						print ('\n#########################################\nMean flux density at '+str(np.mean(np.array(flist)))+' MHz: '+str(flux*1000)+' mJy , Flux error :'+str((fluxerr)*1000)+' mJy, Effective width :'+str(int(weff))+'\n############################################')
					elif fluxmode=='mean_prof':
						print ('\n#########################################\nMean flux density at '+str(np.mean(np.array(flist)))+' MHz: '+str(flux*1000)+' mJy , Flux error :'+str((fluxerr)*1000)+' mJy, Effective width :'+str(int(w_avg))+'\n############################################')
					else:
						print ('Please choose correct flux density method...\n')
						sys.exit()
					freq_list.append(np.mean(np.array(flist)))
					w_list.append(weff)
					flux_list.append(flux)
					flux_err_list.append(fluxerr)
					chi_list.append(chi)
					red_chi_list.append(redchi)
				del data
				del sdata
				del flist
				data=[]
				sdata=[]
				flist=[]
		
	result=np.array([freq_list,flux_list,flux_err_list,red_chi_list])
	np.save(fluxtable,result)
	plt.errorbar(freq_list,[a*1000 for a in flux_list],[b*1000 for b in flux_err_list],fmt='ro',mfc='white',zorder=1,capsize=3)
	plt.xlabel('Frequency (MHz)')
	if fluxmode=='peak':
		plt.ylabel('Peak flux density (mJy)')
	else:	
		plt.ylabel('Mean flux density (mJy)')
	plt.savefig(fluxtable+'.png')
	plt.show()
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
