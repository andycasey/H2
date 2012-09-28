from barak.stats import binomial_min_max_limits

nobs = 9
nsuccess = 2

if 1:
    # we have 6 sub-DLAs, and 2 with H2 what is the lowest mean rate that
    # could plausibly explain this?
      
    # mean number of events in nobs, this is tweaked via trial and error
    # to give p2_or_more = 5%
    mulo = 0.355 #    95% rate
    #mulo = 0.148   #  99% rate (p2_or_more 1%)
    #muhi = 6.5
     
    # probability of observing 2 events given mu
     
    #p2 = poisson.pmf(2, mu)
     
    # probability of observing 2 or more events is the same as 1 -
    # (probability of observing 1 or 0 events), so
     
    p2_or_more = 1 - poisson.cdf(1, mulo)
     
    # below is not really interesting, the quantity we want is the 95%
    # percentile (i.e. go 2.5% in from the smallest mu and highest mu
    # ends)
    #p2_or_less = poisson.cdf(2, muhi)
     
    # Therefore we can say with 95% confidence that low-z sub-DLAs show H2
    # at least 0.355/6 = 6% of the time.
     
    print 'mu %.3f, p(2+) %.5f, nobs %i' % (mulo, p2_or_more, nobs)
    # mulo is the mean number of events per nobs, so mean number for 1 obs
    # is mulo / nobs
    print '95%% confident minimum incidence rate is mu/nobs = %.4f' % (mulo/nobs)

