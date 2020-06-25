import os,psrchive,glob,sys

# Code is written by Devojyoti Kansabanik, NCRA-TIFR, Pune, 23/06/2020

if len(sys.argv)<3:
	print('Code to create standard profile for each time crunched frequency slices. Fitting will be done interactively. Keep all time crunched files in a separate directory without any other files. This is an interactive code. You need to select the full phasebin and profile height manually during the fiiting.\n\n#####################################################\n')
	print ('Usage : python std_profile_maker.py input_path output_path\n#############################################\n')
	print ('input_path : Path of input directory\n')
	print ('output_path : Path of output directory\n#########################################\n')

else:
	path=sys.argv[1]
	outpath=sys.argv[2]
	if os.path.isdir(path)==False:
		raise('Input file path not found')
		sys.exit()
	if os.path.isdir(outpath)==False:
		os.system('mkdir '+outpath)
	os.chdir(path)
	os.system('rm -rf paas*')
	file_list=glob.glob('*')

	for fil in file_list:
		raw=psrchive.Archive_load(fil)
		freq=int(raw.get_frequencies()[0])
		std_file_name='freq_'+str(freq)+'.std'
		do_fit=input('Frequency :'+str(freq)+' Want to do fitting? y/n: ')
		while do_fit=='y':
			centering=input('Want center the profile? y/n: ')
			if centering=='n':
				os.system('paas -i -s '+outpath+'/'+std_file_name+' '+fil)
			else:
				 os.system('paas -i -C -s '+outpath+'/'+std_file_name+' '+fil) 
			repeat=input('Want to repat the fiiting? y/n: ')
			if repeat=='n':
				break
		if do_fit=='y':
			want_save=input('Want to keep this standard profile? y/n: ')
			if want_save=='y':
				print ('Standard profile unloaded : '+outpath+'/'+std_file_name+'\n###########################################\n')
			elif want_save=='n':
				os.system('rm -rf '+outpath+'/'+std_file_name)
		os.system('rm -rf paas*')
	print ('Total '+str(len(file_list))+' standard profiles made and saved at :'+outpath+'\n####################################')




