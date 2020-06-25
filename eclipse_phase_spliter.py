import psrchive,os,glob,sys,numpy as np,matplotlib.pyplot as plt
from scipy.signal import savgol_filter as sv
from scipy.ndimage import gaussian_filter1d as gf
# Code is written by Devojyoti Kansabanik, NCRA-TIFR, Pune, 23/06/2020

if len(sys.argv)<3:
	print ('\n#####################################\n\nEclipse phase spliting code\n\n###############################\n')
	print ('Usage: python eclipse_phase_spliter.py raw_file time_window_smooth\n\n############################\n')
	print ('raw_file : Name of the uncalibrated archive file\n')
	print ('time_window_smooth : Time window in number of sun-integration to smooth the tome profile\n\n###############################\n')
else:
	arch_name=sys.argv[1]
	smooth_window=int(sys.argv[2])

	if smooth_window%2==0:
		smooth_window+=1
	
	print('Loading data....\n')
	arch=psrchive.Archive_load(arch_name)
	print ('####################################\n')
	if not arch.get_dedispersed():
		arch.dedisperse()
	arch.remove_baseline()
	arch.fscrunch()
	arch.centre_max_bin()
	data=arch.get_data()
	data=data.reshape(data.shape[0],data.shape[-1])
	data_prof=np.mean(data,axis=0).flatten()
	width=int(np.sum(data_prof)/max(data_prof))
	onleft=np.argmax(data_prof)-int(width/2)
	onright=np.argmax(data_prof)+int(width/2)
	plt.plot(data_prof/max(data_prof))
	plt.show()
	print ('\n##################\nOn pulse start bin :', onleft,'\n')
	print ('On pulse end bin : ',onright,'\n######################\n')
	data_split=data[:,onleft:onright]
	data_split_avg=sv(np.mean(data_split,axis=1).flatten(),smooth_window,1)
	data_split_avg=gf(data_split_avg,10)
	data_split_avg=sv(np.mean(data_split,axis=1).flatten(),smooth_window,1)
	plt.plot(data_split_avg)
	plt.show()
	start_nsub=int(input('Enter start sub-integration of eclipse phase:'))
	end_nsub=int(input('ENter end sub-integration of eclipse phase:'))
	nsub=arch.get_nsubint()
	pwd=os.getcwd()

	print ('Spliting sub integrations')
	try:
		os.mkdir('split')
	except :
		pass
	os.system('cp -r '+arch_name+' split')
	os.chdir('split')
	print(os.getcwd())
	rawfil=glob.glob('*.ar')[0]
	os.system('psrsplit -q -n 1 '+rawfil)
	os.system('rm -rf '+rawfil)

	file_list=[]
	for i in range(start_nsub,end_nsub):
		file_list.append(glob.glob(rawfil[:3]+'*'+str(i)+'*')[0])

	file_string=''
	for i in file_list:
		file_string+=i+' '
	print ('Adding '+str(len(file_list))+' sub integrations in eclipse phase........\n')
	os.system('psradd '+file_string+' -o ../'+arch_name[:-3]+'_eclipse_phase.ar')
	print ('Eclipse phase archive made...\n')

	for fil in file_list:
		os.system('rm -rf '+fil)
	os.chdir(pwd)
	os.system('rm -rf split')
	print ('Finished eclipse phase spliting .....\n#################################################')





