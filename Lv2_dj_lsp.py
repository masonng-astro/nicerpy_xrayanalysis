import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
import Lv2_swift_lc
from scipy import stats
from scipy.optimize import curve_fit
import glob

def rebin_lc(corr_lc_files,corr_bg_files,bg_scale,tbin,cmpltness):
	"""
	Rebinning the original light curve that was corrected through xrtlccorr

	lc_files - list of corrected light curve files
	tbin - size of new time bins
	cmpltness - level of completeness for the bins
	"""

	"""
	#### time order the lc files!!!
	start_times = [fits.open(lc_files[i])[1].header['TIMEZERO'] for i in range(len(lc_files))]
	time_ordered = np.argsort(start_times)
	ordered_eventlist = np.array(lc_files)[time_ordered]

	times = []
	rates = []
	errors = []
	for i in tqdm(range(len(ordered_eventlist))): #getting every corrected light curve FITS file
		data_arrays = fits.open(ordered_eventlist[i])[1].data
		tstart = fits.open(ordered_eventlist[i])[1].header['TIMEZERO'] #get TIMEZERO
		times += list(tstart+data_arrays['TIME']) #add to an overall list, the 'TIME' array from the corrected light curve FITS file, plus timezero
		rates += list(data_arrays['RATE']) #get the corresponding column for the rates
		errors += list(data_arrays['ERROR']) #get the corresponding column for the errors

	times = np.array(times)
	trunc_times = times-times[0]
	rates = np.array(rates)
	errors = np.array(errors)
	"""
	times,rates,errors,fracexp = Lv2_swift_lc.get_bgsub(corr_lc_files,corr_bg_files,bg_scale)
	trunc_times = times-times[0]

	rebinned_time = []
	rebinned_rate = []
	rebinned_errs = []
	rebinned_fracexp = []
	completeness = []

	time_bins = np.arange(0,trunc_times[-1]+tbin,tbin)
	print('Rebinning...')
	for i in tqdm(range(len(time_bins)-1)):
		time_interval = trunc_times[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
		rate_interval = rates[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
		error_interval = errors[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
		fracexp_interval = fracexp[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]

		comp = len(time_interval)/(tbin/10)

		if len(time_interval) != 0 and comp >= cmpltness:
			mean_time = np.mean(time_interval)
			mean_rate = np.mean(rate_interval)
			mean_error = np.sqrt(np.sum(error_interval**2))/np.size(error_interval)
			sum_fracexp = sum(fracexp_interval)

			rebinned_time.append(mean_time)
			rebinned_rate.append(mean_rate)
			rebinned_errs.append(mean_error)
			rebinned_fracexp.append(sum_fracexp)
			completeness.append(len(time_interval)/(tbin/10))

	#plt.hist(completeness,bins=10)
	#plt.xlabel('Completeness (fraction)',fontsize=12)
	#plt.ylabel('Number',fontsize=12)
	#plt.show()

	#print(len(rebinned_time))

	return np.array(rebinned_time),np.array(rebinned_rate),np.array(rebinned_errs),np.array(rebinned_fracexp)

def psd_error(times,rates,errors):
	"""
	obtain errors for the best frequency estimate of the signal
	"""

	"""
	print(len(times),len(rates),len(errors))
	newdatachoice = np.random.choice(len(times),size=int(0.1*len(times)))

	newtimes = list(np.array([times[0]])) + list(np.array([times[-1]])) + list(times[np.array(list(set(newdatachoice)))])
	newrates = list(np.array([rates[0]])) + list(np.array([rates[-1]])) + list(rates[np.array(list(set(newdatachoice)))])
	newerrs = list(np.array([errors[0]])) + list(np.array([errors[-1]])) + list(errors[np.array(list(set(newdatachoice)))])

	times = newtimes
	rates = newrates
	errors = newerrs
	print(len(times),len(rates),len(errors))
	"""
	freqs_list = []
	psds_list = []
	for j in tqdm(range(1000)):
		new_rates = np.zeros(len(rates))
		for i in range(len(rates)):
			if rates[i] != 0:
				new_rates[i] = np.random.normal(loc=rates[i],scale=errors[i])

		trunc_times = times-times[0]

		newchoice = np.random.choice(len(trunc_times),size=len(trunc_times))
		rand_times = trunc_times[np.array(list(set(newchoice)))]
		rand_rates = new_rates[np.array(list(set(newchoice)))]
		omega,psd,prob3,prob4,prob5 = lsp(rand_times,rand_rates)

		nu_reg = omega/(2.0*np.pi)
		freq = omega/(2*np.pi)
		psds_list.append( np.max(psd[(freq>=8.2e-6)&(freq<=8.4e-6)]) )
		freqs_list.append( freq[psd==psds_list[-1]][0])

	#plt.figure()
	#plt.plot(freq,psd,'rx-')
	#plt.show()

	return freqs_list,psds_list

def lsp(times, rates):
	"""
	cast times and rates as numpy arrays
	"""
	times = np.array(times)
	rates = np.array(rates)
	"""
	Set the initial time to 0
	"""
	tmin = min(times)
	times = times - tmin
	"""
	calculate the number of independent frequencies
	"""
	n0 = len(times)
	Ni = int(-6.362 + 1.193*n0 + 0.00098*n0**2)

	"""
	estimate the minimum and maximum frequencies to be sampled
	"""
	fmin = 1/np.max(times)
	fmax = n0/(2.0*np.max(times))
	#print('Minimum frequency: ' + str(fmin))
	#print('Maximum frequency: ' + str(fmax))
	#print(fmax)
	fmax = 1e-5

	"""
	The frequencies that will be sampled
	"""
	omega = 2*np.pi *(fmin+(fmax-fmin)*np.arange(Ni)/(Ni-1.))

	"""
	subtract the mean from data
	"""
	cn = rates - np.mean(rates)
	"""
	compute the lsp
	"""
	psd = np.zeros(Ni)

	print('Calculating the PSD...')
	for i in tqdm(range(Ni)):
		tau = np.arctan(np.sum(np.sin(2.*omega[i]*times))/np.sum(np.cos(2.*omega[i]*times)))
		tau = tau/(2*omega[i])

		co = np.cos(omega[i]*(times-tau))
		si = np.sin(omega[i]*(times-tau))
		psd[i] = 0.5*(np.sum(cn*co)**2/np.sum(co**2)+np.sum(cn*si)**2/np.sum(si**2))
	"""
	correct normalization
	"""
	var = np.std(cn)
	psd = psd/var**2

	"""
	calculate the false alarm probabilities
	"""

	sigmathree=2.7e-03
	sigmafour=6.33e-05
	sigmafive=5.73e-7

	probthree = -np.log(1.0 - ((1.0-sigmathree)**(1./Ni))  )
	probfour = -np.log(1.0 - ((1.0-sigmafour)**(1./Ni))  )
	probfive = -np.log(1.0 - ((1.0-sigmafive)**(1./Ni))  )

	return omega, psd, probthree, probfour, probfive

def rebin(array, rebinfac):
	if len(array) % rebinfac == 0:
		Nrebin = int(len(array) / rebinfac)
		newarray = [None]*Nrebin
		for i in range(Nrebin):
			startinx = i*rebinfac
			stopinx = startinx + rebinfac
			newarray[i] = np.mean(array[startinx:stopinx])
		return newarray
	else:
		raise TypeError("rebinfac MUST be an integer factor of array size")

def powlaw(x,a,b,c):
	return a + b*x**(-c)

if __name__ == '__main__':
	#eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/lightcurve/ngc300x1_100s.lc'
	#lc_files = glob.glob('/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/lightcurve/sw*bgsub*corr.lc')

	bary_outputfolder = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/lightcurve/'
	obsids = [str(i) for i in range(49834027,49834042)] + [str(i) for i in range(49834043,49834062)] + [str(i) for i in range(49834063,49834066)] + ['88810002'] + [str(i) for i in range(49834066,49834069)] + [str(i) for i in range(49834070,49834079)] + [str(i) for i in range(49834080,49834088)]
	corr_lc_files = [bary_outputfolder + 'sw000' + obsids[i] + '_corr.lc' for i in range(len(obsids))]
	corr_ulx1_files = [bary_outputfolder + 'sw000' + obsids[i] + '_ulx1_corr.lc' for i in range(len(obsids))]
	corr_bg_files = [bary_outputfolder + 'sw000' + obsids[i] + '_bg_corr.lc' for i in range(len(obsids))]
	bg_scale_x1 = (30/120)**2
	#bg_scale_ulx1 = (35/120)**2
	rebinned_t, rebinned_rate, rebinned_err, rebinned_fracexp = rebin_lc(corr_lc_files,corr_bg_files,bg_scale_x1,100,0)
	#print(rebinned_fracexp[:50])
	#rebinned_t_ulx1, rebinned_rate_ulx1, rebinned_err_ulx1 = rebin_lc(corr_ulx1_files,corr_bg_files,bg_scale_ulx1,3600)

	#x1_text = open(bary_outputfolder + 'ngc300x1_bg_exp_corr_lc_3600s.txt','w')
	#ulx1_text = open(bary_outputfolder + 'ngc300ulx1_bg_exp_corr_lc_3600s.txt','w')

	tstart_49834027 = 546830295.758713

	#for i in range(len(rebinned_t)):
	#	x1_text.write(str(51910 + 7.428703700000000E-04+(rebinned_t[i]+tstart_49834027)/86400) + ' ' + str(rebinned_rate[i]) + ' ' + str(rebinned_err[i]) + '\n')
	#x1_text.close()

	#for i in range(len(rebinned_t_ulx1)):
	#	ulx1_text.write(str(51910 + 7.428703700000000E-04 + (rebinned_t_ulx1[i]+tstart_49834027)/86400) + ' ' + str(rebinned_rate_ulx1[i]) + ' ' + str(rebinned_err_ulx1[i]) + '\n')
	#ulx1_text.close()

	#print(51910 + 7.428703700000000E-04+(tstart_49834027+rebinned_t[rebinned_err>=0.06])/86400,rebinned_rate[rebinned_err>=0.06],rebinned_err[rebinned_err>=0.06])
	#print(51910 + 7.428703700000000E-04+(tstart_49834027+rebinned_t_ulx1[rebinned_err_ulx1>=0.06])/86400,rebinned_rate_ulx1[rebinned_err_ulx1>=0.06],rebinned_err_ulx1[rebinned_err_ulx1>=0.06])

	#mjd = 51910 + 7.428703700000000E-04 + (tstart_49834027+rebinned_t)/86400
	#mjd_ulx1 = 51910 + 7.428703700000000E-04 + (tstart_49834027+rebinned_t_ulx1)/86400
	#plt.errorbar(x=mjd[rebinned_err<=0.06],y=rebinned_rate[rebinned_err<=0.06],yerr=rebinned_err[rebinned_err<=0.06],fmt='rx')
	#plt.errorbar(x=mjd_ulx1[rebinned_err_ulx1<=0.06],y=rebinned_rate_ulx1[rebinned_err_ulx1<=0.06],yerr=rebinned_err_ulx1[rebinned_err_ulx1<=0.06],fmt='bx')
	#plt.legend(('X-1','ULX-1'),fontsize=12)
	#plt.xlabel('Time (MJD)',fontsize=12)
	#plt.ylabel('[Exposure-corrected] Count rate (c/s)',fontsize=12)
	#plt.axhline(y=0,color='k',lw=0.5,alpha=0.5)
	#plt.show()

	#psd_error(times,rates,errors)
	#freqs_list, psd_list = psd_error(rebinned_t,rebinned_rate,rebinned_err)

	#print('Median frequency: ' + str(np.median(freqs_list)))
	#print('Error in frequency: ' + str(np.std(freqs_list)))
	#print('Powers: ' + str(psd_list))

	mjd_x1 = 51910 + 7.428703700000000E-04 + (tstart_49834027+rebinned_t)/86400

	xmm_lc1 = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_lccorr.lc'
	rebinned_t_xmm1 = fits.open(xmm_lc1)[1].data['TIME']
	rebinned_rate_xmm1 = fits.open(xmm_lc1)[1].data['RATE']
	xmm_lc2 = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_lccorr.lc'
	rebinned_t_xmm2 = fits.open(xmm_lc2)[1].data['TIME']
	rebinned_rate_xmm2 = fits.open(xmm_lc2)[1].data['RATE']

	mjd_x1_xmm = fits.open(xmm_lc1)[1].header['MJDREF'] + np.array(list(rebinned_t_xmm1) + list(rebinned_t_xmm2))/86400
	rebinned_rate_xmm = np.array(list(rebinned_rate_xmm1) + list(rebinned_rate_xmm2))

	total_t = np.array(list(mjd_x1_xmm) + list(mjd_x1))
	total_t_s = (total_t - total_t[0])*86400
	total_rate = np.array(list(rebinned_rate_xmm) + list(rebinned_rate))

	### For the light curve
	plt.figure()
	plt.plot(total_t_s,total_rate,'r-')
	plt.plot(total_t_s,0.2*np.sin(2*np.pi*8.4712e-6 * total_t_s + 0.2)+0.2,'b-')

	### For the PSD
	omega,psd,prob3,prob4,prob5 = lsp(total_t_s,total_rate)
	#omega,psd,prob3,prob4,prob5 = lsp(times,rates)
	freq = omega/(2*np.pi)

	plt.figure()
	plt.plot(freq,psd,'rx-')
	#plt.yscale('log')
	#plt.xscale('log')
	plt.xlabel('Frequency (Hz)',fontsize=12)
	plt.ylabel('Normalized Power',fontsize=12)
	plt.axhline(y=prob3,lw=0.5,alpha=0.5)
	plt.axhline(y=prob4,lw=0.5,alpha=0.5)
	plt.axhline(y=prob5,lw=0.5,alpha=0.5)
	print(prob3,prob4,prob5)

	plt.show()
















	"""
	df = pd.read_csv("lc0853980201.dat",delim_whitespace=True, header=None)
	times = df.iloc[:,0]
	rates = df.iloc[:,1]

	np.random.shuffle(rates)

	omega, psd,_,_,_ = lsp(times, rates)
	nu_reg = omega/(2.0*np.pi)

	refac = 2
	if refac > 1:
		hot = rebin(nu_reg,refac)
		pot = rebin(psd, refac)
		zot = pot/np.sqrt(refac)
	else:
		hot = nu_reg
		pot = psd
		zot = pot
	print(len(psd))
	print(np.mean(psd))
	plt.step(hot, pot,where='mid',color='orange')
	plt.errorbar(hot,pot,yerr=zot, ls='none', color='orange')
	plt.xscale('log')
	plt.yscale('log')
	popt, pcov = curve_fit(powlaw, hot, pot)
	plt.plot(nu_reg, powlaw(nu_reg, *popt), 'g--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
	plt.legend(loc='best')
	#plt.show()
	"""
