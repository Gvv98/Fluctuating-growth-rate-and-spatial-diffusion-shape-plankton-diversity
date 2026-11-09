import numpy as np
import scipy as sp
import pandas as pd
import warnings
import copy

# expression for the SAD (Generalised Inverse Gaussian distribution)
def logGIG(n,b,c,l):

	return(sp.stats.geninvgauss.logpdf(n,1.0-l,np.sqrt(4.0*b/c),scale=np.sqrt(b*c)))

# loglikelihood
def loglik(fit_args,n_array):

	b = fit_args[0]
	c = fit_args[1]
	l = fit_args[2]

	return(-np.sum(logGIG(n_array,b,c,l)))

# maximise the loglikelihood
def minimise_loglik(fit_init,n_array):

	return(sp.optimize.minimize(loglik,fit_init,args=n_array,bounds=[(0,None),(0,None),(0,None)],method="Nelder-Mead",options={'maxiter':1000}).x)

# calculate the Kolmogorov-Smirnov statistic for the GIG
def KS_GIG(n_array,fit_args):

	counts = np.asarray(pd.Series(n_array).value_counts().sort_index(ascending=True))
	d = {'x':np.unique(n_array),'S':np.cumsum(counts)/np.sum(counts)}
	# create a dataframe
	df_cumul = pd.DataFrame(data=d)
	# calculate the theoretical CDF
	df_cumul['P'] = sp.stats.geninvgauss.cdf(df_cumul['x'].values,1.0-fit_args[2],np.sqrt(4.0*fit_args[0]/fit_args[1]),scale=np.sqrt(fit_args[0]*fit_args[1]))
	# return the KS statistic = maximum difference 
	# between the theoretical and empirical CDF
	return(np.max(abs(df_cumul['S'].values-df_cumul['P'].values)))

# calculate the Anderson-Darling statistic for the GIG
def AD_GIG(n_array,fit_args):

	b = fit_args[0]
	c = fit_args[1]
	l = fit_args[2]
	# initialise the S sum
	S_stat = 0
	# sort the abundances
	n_array = np.asarray(sorted(n_array))
	# calculate the S sum
	for i in range(1,len(n_array)+1):
		S_stat += (2*i-1)*(sp.stats.geninvgauss.logcdf(n_array[i-1],1.0-l,np.sqrt(4.0*b/c),scale=np.sqrt(b*c)) + 
					np.log(1.0-sp.stats.geninvgauss.cdf(n_array[len(n_array)-i],1.0-l,np.sqrt(4.0*b/c),scale=np.sqrt(b*c))))
	# return the AD statistic by correcting the S sum
	return(np.sqrt(-len(n_array)-S_stat/len(n_array)))	

# calculate the Cramer-von Mises statistic for the GIG
def CvM_GIG(n_array,fit_args):

	# initialise the S sum
	S_stat = 0
	# sort the abundances
	n_array = np.asarray(sorted(n_array))
	# calculate the S sum
	for i in range(1,len(n_array)+1):
		S_stat += ((2.0*i-1.0)/2.0/len(n_array) - 
					sp.stats.geninvgauss.cdf(n_array[i-1],1.0-fit_args[2],np.sqrt(4.0*fit_args[0]/fit_args[1]),scale=np.sqrt(fit_args[0]*fit_args[1])))**2
	# return the AD statistic by correcting the S sum
	return(1.0/12.0/len(n_array) + S_stat)

# calculate the Kolmogorov-Smirnov statistic for the lognormal
def KS_LNM(n_array,fit_args):

	counts = np.asarray(pd.Series(n_array).value_counts().sort_index(ascending=True))
	d = {'x':np.unique(n_array),'S':np.cumsum(counts)/np.sum(counts)}
	# create a dataframe
	df_cumul = pd.DataFrame(data=d)
	# calculate the theoretical CDF
	df_cumul['P'] = sp.stats.lognorm.cdf(df_cumul['x'].values,s=fit_args[0],scale=fit_args[1])
	# return the KS statistic = maximum difference 
	# between the theoretical and empirical CDF
	return(np.max(abs(df_cumul['S'].values-df_cumul['P'].values)))

# calculate the Anderson-Darling statistic for the lognormal
def AD_LNM(n_array,fit_args):

	# initialise the S sum
	S_stat = 0
	# sort the abundances
	n_array = np.asarray(sorted(n_array))
	# calculate the S sum
	for i in range(1,len(n_array)+1):
		S_stat += (2*i-1)*(sp.stats.lognorm.logcdf(n_array[i-1],s=fit_args[0],scale=fit_args[1]) + 
					np.log(1.0-sp.stats.lognorm.cdf(n_array[len(n_array)-i],s=fit_args[0],scale=fit_args[1])))
	# return the AD statistic by correcting the S sum
	return(np.sqrt(-len(n_array)-S_stat/len(n_array)))	

# calculate the Cramer-von Mises statistic for the lognormal
def CvM_LNM(n_array,fit_args):

	# initialise the S sum
	S_stat = 0
	# sort the abundances
	n_array = np.asarray(sorted(n_array))
	# calculate the S sum
	for i in range(1,len(n_array)+1):
		S_stat += ((2.0*i-1.0)/2.0/len(n_array) - 
					sp.stats.lognorm.cdf(n_array[i-1],s=fit_args[0],scale=fit_args[1]))**2
	# return the AD statistic by correcting the S sum
	return(1.0/12.0/len(n_array) + S_stat)

# sample from the fitted GIG distribution and refit
def replicate_GIG(n_array,fit_args):

	sub_array = sp.stats.geninvgauss.rvs(1.0-fit_args[2],np.sqrt(4.0*fit_args[0]/fit_args[1]),scale=np.sqrt(fit_args[0]*fit_args[1]),size=len(n_array))
	res = minimise_loglik(fit_args,sub_array)
	
	return([np.sum(sub_array),np.mean(sub_array),np.var(sub_array),np.min(sub_array),res[0],res[1],res[2],KS_GIG(sub_array,res),AD_GIG(sub_array,res),CvM_GIG(sub_array,res)])

# sample from the fitted lognormal distribution and refit
def replicate_LNM(n_array,fit_args):

	sub_array = sp.stats.lognorm.rvs(s=fit_args[0],scale=fit_args[1],size=len(n_array))
	res = [np.std(np.log(sub_array)),np.exp(np.mean(np.log(sub_array)))]
	
	return([np.sum(sub_array),np.mean(sub_array),np.var(sub_array),np.min(sub_array),0.0,0.0,0.0,KS_LNM(sub_array,res),AD_LNM(sub_array,res),CvM_LNM(sub_array,res)])

# gather the statistics for one sample in the GRUMP data
def do_sample(n_array,Nrep):

	# initial parameter values
	fit_init = [np.min(n_array),np.var(n_array)/np.mean(n_array),1.5]

	replicates = np.zeros((Nrep,10))

	# fit the GIG distribution
	res_GIG = minimise_loglik(fit_init,n_array)
	# fit the lognormal distribution (method of the moments)
	res_LNM = [np.std(np.log(n_array)),np.exp(np.mean(np.log(n_array)))]
	# get the Bayesian Information Criterion for both fits
	BIC_GIG = 3.0*np.log(len(n_array))+2.0*loglik(res_GIG,n_array)
	BIC_LNM = 2.0*np.log(len(n_array))-2.0*np.sum(np.log(sp.stats.lognorm.pdf(n_array,s=res_LNM[0],scale=res_LNM[1])))
	
	# if the GIG is a better fit (BIC), 
	# do the goodness-of-fit sampling with the GIG
	if(BIC_GIG < BIC_LNM):
		replicates[0,:] = [np.sum(n_array),np.mean(n_array),np.var(n_array),np.min(n_array),res_GIG[0],res_GIG[1],res_GIG[2],KS_GIG(n_array,res_GIG),AD_GIG(n_array,res_GIG),CvM_GIG(n_array,res_GIG)]

		for i in range(1,Nrep):

			replicates[i,:] = replicate_GIG(n_array,res_GIG)

		return([replicates[0,0],len(n_array),replicates[0,3],BIC_GIG,BIC_LNM,
				replicates[0,1],np.quantile(replicates[:,1],0.025),np.quantile(replicates[:,1],0.975),
				replicates[0,2],np.quantile(replicates[:,2],0.025),np.quantile(replicates[:,2],0.975),
				np.mean(replicates[1:,0]),np.quantile(replicates[1:,0],0.025),np.quantile(replicates[1:,0],0.975),
				np.mean(replicates[1:,3]),np.quantile(replicates[1:,3],0.025),np.quantile(replicates[1:,3],0.975),
				replicates[0,4],np.quantile(replicates[:,4],0.025),np.quantile(replicates[:,4],0.975),
				replicates[0,5],np.quantile(replicates[:,5],0.025),np.quantile(replicates[:,5],0.975),
				replicates[0,6],np.quantile(replicates[:,6],0.025),np.quantile(replicates[:,6],0.975),
				replicates[0,7],len(replicates[:,7][replicates[:,7] >= replicates[0,7]])/Nrep,
				replicates[0,8],len(replicates[:,8][replicates[:,8] >= replicates[0,8]])/Nrep,
				replicates[0,9],len(replicates[:,9][replicates[:,9] >= replicates[0,9]])/Nrep])

	# if the lognormal is a better fit (BIC), 
	# do the goodness-of-fit sampling with the lognormal
	else:
		replicates[0,:] = [np.sum(n_array),np.mean(n_array),np.var(n_array),np.min(n_array),res_GIG[0],res_GIG[1],res_GIG[2],KS_LNM(n_array,res_LNM),AD_LNM(n_array,res_LNM),CvM_LNM(n_array,res_LNM)]

		for i in range(1,Nrep):

			replicates[i,:] = replicate_LNM(n_array,res_LNM)


		return([replicates[0,0],len(n_array),replicates[0,3],BIC_GIG,BIC_LNM,
				replicates[0,1],np.quantile(replicates[:,1],0.025),np.quantile(replicates[:,1],0.975),
				replicates[0,2],np.quantile(replicates[:,2],0.025),np.quantile(replicates[:,2],0.975),
				np.mean(replicates[1:,0]),np.quantile(replicates[1:,0],0.025),np.quantile(replicates[1:,0],0.975),
				np.mean(replicates[1:,3]),np.quantile(replicates[1:,3],0.025),np.quantile(replicates[1:,3],0.975),
				replicates[0,4],0.0,0.0,
				replicates[0,5],0.0,0.0,
				replicates[0,6],0.0,0.0,
				replicates[0,7],len(replicates[:,7][replicates[:,7] >= replicates[0,7]])/Nrep,
				replicates[0,8],len(replicates[:,8][replicates[:,8] >= replicates[0,8]])/Nrep,
				replicates[0,9],len(replicates[:,9][replicates[:,9] >= replicates[0,9]])/Nrep])

# do all samples in the data for a certain type (16S / 18S)
def do_type(type,abund,Nrep):

	abund_type = copy.deepcopy(abund[abund["Sequence_Type"]==type])
	abund_type.drop(columns=["Sequence_Type"],inplace=True)
	
	vec_cruise = np.unique(abund_type["cruise"].values)

	Nsample = len(abund_type[["SampleID","cruise"]].drop_duplicates().index)

	out = np.zeros((Nsample,32))
	vec_cruise_out = [""]*Nsample
	vec_sample_out = [""]*Nsample

	k = 0
	for i in range(len(vec_cruise)):

		abund_sub = abund_type[abund_type["cruise"]==vec_cruise[i]]
		vec_sample = np.unique(abund_sub["SampleID"].values)

		for j in range(len(vec_sample)):

			print(k)

			n_array = abund_sub[abund_sub["SampleID"]==vec_sample[j]]["Raw_Sequence_Counts"].values
			
			if(len(n_array) < 50):
				out[k,:] = np.concatenate(([np.sum(n_array),len(n_array),np.min(n_array)],np.repeat(-999.0,29)))
			else:
				with warnings.catch_warnings():
					warnings.simplefilter("ignore")
					out[k,:] = do_sample(n_array,Nrep)

			vec_cruise_out[k] = vec_cruise[i]
			vec_sample_out[k] = vec_sample[j]
			k += 1

	df_out = pd.DataFrame(out,columns=["N","S","Nmin","BIC_GIG","BIC_LNM",
								"Nmean","Nmean_lo","Nmean_up",
								"Nvar","Nvar_lo","Nvar_up",
								"N_mn","N_lo","N_up",
								"Nmin_mn","Nmin_lo","Nmin_up",
								"b","b_lo","b_up","c","c_lo","c_up",
								"l","l_lo","l_up","KS","p_KS","AD","p_AD",
								"CvM","p_CvM"])
	df_out["sample"] = vec_sample_out
	df_out["cruise"] = vec_cruise_out
	df_out["type"] = type[(len(type)-3):]

	df_out = df_out.astype({"N":np.int64,"S":np.int64,"Nmin":np.int64})
	return(df_out[["type","cruise","sample","N","S","Nmin","BIC_GIG","BIC_LNM",
			"Nmean","Nmean_lo","Nmean_up",
			"Nvar","Nvar_lo","Nvar_up",
			"N_mn","N_lo","N_up",
			"Nmin_mn","Nmin_lo","Nmin_up",
			"b","b_lo","b_up","c","c_lo","c_up",
			"l","l_lo","l_up","KS","p_KS","AD","p_AD",
			"CvM","p_CvM"]])

def main():

	# number of replicates in the goodness-of-fit tests
	# minimal p-value = 1/Nrep
	Nrep = 100

	abund = pd.read_csv("GRUMP_data.csv",usecols=["cruise","SampleID","Raw_Sequence_Counts","Sequence_Type","depth"],
							dtype={"SampleID":str})
	# only keep surface stations (depth < 20 metres)
	abund = abund[abund["depth"] < 20.0].drop(columns=["depth"])
	# do the fitting and goodness-of-fit tests for 16S
	df_16S = do_type("Prokaryotic_16S",abund,Nrep)
	# do the fitting and goodness-of-fit tests for 18S
	df_18S = do_type("Eukaryote_18S",abund,Nrep)

	# combine the results and write to csv file
	pd.concat([df_16S,df_18S]).to_csv("fit_goodness_GRUMP.csv",index=False)
	

if __name__ == "__main__":
	main()