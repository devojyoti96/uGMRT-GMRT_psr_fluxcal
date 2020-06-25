import numpy as np,matplotlib.pyplot as plt,sys
from scipy.optimize import curve_fit
from astropy.modeling import models,fitting

def power_law(x,a,b):
	return a*(x**b)

nonec_file_list=sys.argv[1]
ec_file_list=sys.argv[2]
outname=sys.argv[3]

nonec_list=nonec_file_list.split(',')
ec_list=ec_file_list.split(',')
combined_nonec=[]
combined_nonec_freq=[]
combined_nonec_err=[]
for nonec in nonec_list:
	print (nonec)
	try:
		nonec_data=np.load(nonec,allow_pickle=True,encoding='bytes')
	except:
		nonec_data=np.load(nonec)

	# Noneclipse plot
	freq=nonec_data[0].astype('float')
	flux=nonec_data[1].astype('float')*1000
	flux_err=nonec_data[2].astype('float')*1000
	for i in range(len(freq)):
		combined_nonec_freq.append(freq[i])
		combined_nonec.append(flux[i])
		combined_nonec_err.append(flux_err[i])
	plt.errorbar(freq,flux,flux_err,fmt='o',mfc='white',zorder=1,capsize=5,label=nonec.split('/')[-1][:-4])
	
combined_nonec_freq=np.array(combined_nonec_freq)
combined_nonec=np.array(combined_nonec)
combined_nonec_err=np.array(combined_nonec_err)

combined_ec_freq=[]
combined_ec_flux=[]
combined_ec_fluxerr=[]

for ec in ec_list:
	# Eclipse spectra plot
	try:
		ec_data=np.load(ec,allow_pickle=True,encoding='bytes')
	except:
		ec_data=np.load(ec)
	freq=ec_data[0]
	flux=ec_data[1]*1000
	flux_err=ec_data[2]*1000
	for i in range(len(freq)):
		combined_ec_freq.append(freq[i])
		combined_ec_flux.append(flux[i])
		combined_ec_fluxerr.append(flux_err[i])
	plt.errorbar(freq,flux,flux_err,fmt='go',mfc='white',zorder=1,capsize=5,label=ec.split('/')[-1][:-4])

# Ploting labels and legends
plt.legend(loc='lower left')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Mean flux density (mJy)')
plt.savefig(outname+'_spectra.png')
plt.show()


for ec in ec_list:
	try:
		ec_data=np.load(ec,allow_pickle=True,encoding='bytes')
	except:
		ec_data=np.load(ec)
	freq=ec_data[0]
	flux=ec_data[1]*1000
	flux_err=ec_data[2]*1000
	logfreq_ec=[np.log10(i) for i in freq]
	logflux_ec=[np.log10(i) for i in flux]
	logfluxerr_ec=(flux_err/flux)/2.7182818
	plt.errorbar(logfreq_ec,logflux_ec,logfluxerr_ec,fmt='o',mfc='white',zorder=1,capsize=5,label=ec.split('/')[-1][:-4])

for nonec in nonec_list:
	try:
		nonec_data=np.load(nonec,allow_pickle=True,encoding='bytes')
	except:
		nonec_data=np.load(nonec)
	freq=nonec_data[0].astype('float')
	flux=nonec_data[1].astype('float')*1000
	flux_err=nonec_data[2].astype('float')*1000	
	logfreq_nonec=[np.log10(i) for i in freq]
	logflux_nonec=[np.log10(i) for i in flux]
	logfluxerr_nonec=(flux_err/flux)/2.7182818
	plt.errorbar(logfreq_nonec,logflux_nonec,logfluxerr_nonec,fmt='o',mfc='white',zorder=1,capsize=5,label=nonec.split('/')[-1][:-4])

plt.xlabel('log(Frequency)')
plt.ylabel('log(Mean flux density)')
plt.legend(loc='lower left')
plt.savefig(outname+'_logspectra.png')
plt.show()

from scipy.optimize import curve_fit
def power_law(x,b,alpha,beta,nu_c):
	y=[]
	for i in x:
		if i <=nu_c:
			y.append(float(i**alpha))
		else:
			y.append(float((i**beta)*(nu_c**(alpha-beta))))
	y=np.array(y)
	y.astype('float')
	return y
#	return b*(freq**alpha)*np.exp((alpha/beta)*freq**(-beta))
'''
#plt.plot(np.linspace(250,450,1024),power_law(np.linspace(250,450),1,)
popt,pcov=curve_fit(power_law,combined_nonec_freq,combined_nonec,absolute_sigma=True)
print (popt,np.sqrt(np.diag(pcov)))

x=np.linspace(min(combined_nonec_freq),max(combined_nonec_freq))
y=power_law(x,popt[0],popt[1],popt[2],popt[3])
plt.plot(np.log10(x),np.log10(y))
plt.errorbar(np.log10(combined_nonec_freq),np.log10(combined_nonec),np.divide(np.array(combined_nonec_err),np.array(combined_nonec))/2.7182818,fmt='o',capsize=5)
plt.show()


'''
f = models.SmoothlyBrokenPowerLaw1D(amplitude=12, x_break=290,alpha_1=-1.5, alpha_2=2.6,delta=0.01)
fit= fitting.LevMarLSQFitter()
fitted_line = fit(f, combined_nonec_freq,combined_nonec)
plt.errorbar((combined_nonec_freq),(combined_nonec),(combined_nonec_err),fmt='o',capsize=5)
x=np.linspace(min(combined_nonec_freq),max(combined_nonec_freq))
alpha1=fitted_line.alpha_1.value
alpha2=fitted_line.alpha_2.value
error=(np.diag(fit.fit_info['param_cov']))
print(fitted_line,error)
alpha1_err=error[2]
alpha2_err=error[3]
plt.plot((x),(fitted_line(x)))
plt.ylim(bottom=0)
plt.xlabel('(Frequency)')
plt.ylabel('(Mean flux density)')
plt.text(250,50,r'$\alpha_1=$'+str(-1*alpha1)[:4]+r'$\pm$'+str(alpha1_err)[:4])
plt.text(357,50,r'$\alpha_2=$'+str(-1*alpha2)[:4]+r'$\pm$'+str(alpha2_err)[:4])
plt.savefig(outname+'_logspectra_fit.png')
plt.show()



