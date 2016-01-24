import nitime
#Import the time-series objects:
from nitime.timeseries import TimeSeries
#Import the analysis objects:
from nitime.analysis import CorrelationAnalyzer, CoherenceAnalyzer
#Import utility functions:
from nitime.utils import percent_change
from nitime.viz import drawmatrix_channels, drawgraph_channels, plot_xcorr

#This information (the sampling interval) has to be known in advance:
s_f = 40000
f_lb = 0.02 # low freq cut off.
f_ub = 100 #Hz high freq cut off


TR = float(h.dt)#already is float, just acknowledgement.



n_samples = int(tstop/h.dt)+2 #data_rec.shape[0]
#nseq=int(len(allrows))-1
nseq=numcell
data = np.zeros((numcell, numcell))
#Matrix = [[0 for x in xrange(0,int(num_cells+1))] for x in xrange(0,int(num_cells+1))]
data2= [[0 for x in xrange(0,int(numcell))] for x in xrange(0,int(numcell))]
numcell=int(numcell)
for i in xrange(1,numcell):
 for j in xrange(1,numcell):
  x=zip(np.array(time_courses[int(i)]),np.array(time_courses[int(j)]))
  data2[i][j]=statsmodels.tsa.stattools.  causalitytests(x, maxlag, addconst=True, verbose=True)
  print data2[i,j]
 
 #s=allrows[i]
 #if(int(len(s))==9): #//If there are enough columns in the file. Ie there should be nine 
  data[i,:] = np.array(time_courses[int(i)])


























def grangercausalitytests(x, maxlag, addconst=True, verbose=True):
    '''four tests for granger non causality of 2 timeseries

    all four tests give similar results
    `params_ftest` and `ssr_ftest` are equivalent based on F test which is
    identical to lmtest:grangertest in R

    Parameters
    ----------
    x : array, 2d, (nobs,2)
        data for test whether the time series in the second column Granger
        causes the time series in the first column
    maxlag : integer
        the Granger causality test results are calculated for all lags up to
        maxlag
    verbose : bool
        print results if true

    Returns
    -------
    results : dictionary
        all test results, dictionary keys are the number of lags. For each
        lag the values are a tuple, with the first element a dictionary with
        teststatistic, pvalues, degrees of freedom, the second element are
        the OLS estimation results for the restricted model, the unrestricted
        model and the restriction (contrast) matrix for the parameter f_test.

    Notes
    -----
    TODO: convert to class and attach results properly

    The Null hypothesis for grangercausalitytests is that the time series in
    the second column, x2, does NOT Granger cause the time series in the first
    column, x1. Grange causality means that past values of x2 have a
    statistically significant effect on the current value of x1, taking past
    values of x1 into account as regressors. We reject the null hypothesis
    that x2 does not Granger cause x1 if the pvalues are below a desired size
    of the test.

    The null hypothesis for all four test is that the coefficients
    corresponding to past values of the second time series are zero.

    'params_ftest', 'ssr_ftest' are based on F distribution

    'ssr_chi2test', 'lrtest' are based on chi-square distribution

    References
    ----------
    http://en.wikipedia.org/wiki/Granger_causality
    Greene: Econometric Analysis

    '''
    from scipy import stats # lazy import

    x = np.asarray(x)

    resli = {}

    for mlg in range(1, maxlag+1):
        result = {}
        if verbose:
            print '\nGranger Causality'
            print 'number of lags (no zero)', mlg
        mxlg = mlg #+ 1 # Note number of lags starting at zero in lagmat

        # create lagmat of both time series
        dta = lagmat2ds(x, mxlg, trim='both', dropex=1)

        #add constant
        if addconst:
            dtaown = add_constant(dta[:,1:mxlg+1], prepend=False)
            dtajoint = add_constant(dta[:,1:], prepend=False)
        else:
            raise ValueError('Not Implemented')
            dtaown = dta[:,1:mxlg]
            dtajoint = dta[:,1:]

        #run ols on both models without and with lags of second variable
        res2down = OLS(dta[:,0], dtaown).fit()
        res2djoint = OLS(dta[:,0], dtajoint).fit()

        #print results
        #for ssr based tests see: http://support.sas.com/rnd/app/examples/ets/granger/index.htm
        #the other tests are made-up

        # Granger Causality test using ssr (F statistic)
        fgc1 = (res2down.ssr-res2djoint.ssr)/res2djoint.ssr/(mxlg)*res2djoint.df_resid
        if verbose:
            print 'ssr based F test:         F=%-8.4f, p=%-8.4f, df_denom=%d, df_num=%d' % \
              (fgc1, stats.f.sf(fgc1, mxlg, res2djoint.df_resid), res2djoint.df_resid, mxlg)
        result['ssr_ftest'] = (fgc1, stats.f.sf(fgc1, mxlg, res2djoint.df_resid), res2djoint.df_resid, mxlg)

        # Granger Causality test using ssr (ch2 statistic)
        fgc2 = res2down.nobs*(res2down.ssr-res2djoint.ssr)/res2djoint.ssr
        if verbose:
            print 'ssr based chi2 test:   chi2=%-8.4f, p=%-8.4f, df=%d' %  \
              (fgc2, stats.chi2.sf(fgc2, mxlg), mxlg)
        result['ssr_chi2test'] = (fgc2, stats.chi2.sf(fgc2, mxlg), mxlg)

        #likelihood ratio test pvalue:
        lr = -2*(res2down.llf-res2djoint.llf)
        if verbose:
            print 'likelihood ratio test: chi2=%-8.4f, p=%-8.4f, df=%d' %  \
              (lr, stats.chi2.sf(lr, mxlg), mxlg)
        result['lrtest'] = (lr, stats.chi2.sf(lr, mxlg), mxlg)

        # F test that all lag coefficients of exog are zero
        rconstr = np.column_stack((np.zeros((mxlg-1,mxlg-1)), np.eye(mxlg-1, mxlg-1),\
                                   np.zeros((mxlg-1, 1))))
        rconstr = np.column_stack((np.zeros((mxlg,mxlg)), np.eye(mxlg, mxlg),\
                                   np.zeros((mxlg, 1))))
        ftres = res2djoint.f_test(rconstr)
        if verbose:
            print 'parameter F test:         F=%-8.4f, p=%-8.4f, df_denom=%d, df_num=%d' % \
              (ftres.fvalue, ftres.pvalue, ftres.df_denom, ftres.df_num)
        result['params_ftest'] = (np.squeeze(ftres.fvalue)[()],
                                  np.squeeze(ftres.pvalue)[()],
                                  ftres.df_denom, ftres.df_num)

        resli[mxlg] = (result, [res2down, res2djoint, rconstr])

    return resli
