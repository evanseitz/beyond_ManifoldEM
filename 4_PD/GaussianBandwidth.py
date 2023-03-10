import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
import warnings
warnings.simplefilter(action='ignore', category=OptimizeWarning)

# =============================================================================
# D: Distance matrix
# logEps: Range of values to try
# =============================================================================
# History:   A. Dashti (UWM, 2016); adapted from Chuck (2011)
#            See ManifoldEM Matlab repository for syntax parallels. As well, a similar...
#            ...workflow will be publically released via the ManifoldEM Python suite...
#            ...(estimated 2021) with that code translated from Matlab to Python...
#            ...by H. Liao (CU, 2019) and modified therein by E. Seitz (CU, 2020).
# =============================================================================

def fun(xx, aa0, aa1, aa2, aa3):
    #aa3: y-value of tanh inflection point
    #aa2: y-value of apex (asymptote)
    #aa1: x-value of x-shift (inverse sign)
    #aa0: alters spread
    F = aa3 + aa2 * np.tanh(aa0 * xx + aa1)
    return F

def find_threshold(logEps, D2):
    eps = np.exp(logEps)
    d = 1. / (2. * np.max(eps)) * D2
    sg = np.sort(d)
    ss = np.sum(np.exp(-sg))
    thr = max(-np.log(0.01 * ss / len(D2)), 10)  # taking 1% of the average (10)
    return thr

def op(D,logEps,a0):
    # range of values to try:
    logSumWij = np.zeros(len(logEps))
    D2 = D * D
    thr = find_threshold(logEps, D2)
    for k in range(len(logEps)):
        eps = np.exp(logEps[k])
        d = 1. / (2. * eps) * D2
        d = -d[d < thr]
        Wij = np.exp(d) #see Coifman 2008
        logSumWij[k] = np.log(sum(Wij))

    # Curve fitting of a tanh():
    resnorm = np.inf
    cc = 0
    while (resnorm > 100):
        cc += 1
        popt, pcov = curve_fit(fun, logEps, logSumWij, p0=a0)#, maxfev=10000)
        resnorm = sum(np.sqrt(np.fabs(np.diag(pcov))))
        a0 = 1 * (np.random.rand(4, 1) - .5)
        residuals = logSumWij - fun(logEps, popt[0], popt[1], popt[2], popt[3])
        ss_res = np.sum(residuals**2) #residual sum of squares
        ss_tot = np.sum((logSumWij - np.mean(logSumWij))**2) #total sum of squares
        R_squared = 1 - (ss_res / ss_tot) #coefficient of determination

    return (popt, logSumWij, resnorm, R_squared)