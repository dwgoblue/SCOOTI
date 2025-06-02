import pandas as pd
import numpy as np



def hypergeom_test(M, N, n, k):
    """
    arr1: array with only 1, 0, -1
    arr2: array with only 1, 0, -1
    """
    from scipy.stats import hypergeom
    #M = len(arr1)
    #N = sum(arr1==label)
    #n = sum(arr2==label)
    #k = sum((arr1==label)*(arr2==label))
    # more overlaps higher pvalues
    hypergeom_p = hypergeom.sf(k-1, M, n, N)
    a, b, c, d = k, n-k, N-k, M-(n+N)+k
    table = np.array([[a, b], [c, d]])
    start, end = hypergeom.support(M, n, N)
    # print(sum(hypergeom.pmf(np.arange(k, end+1), M, n, N)))
    from scipy.stats import fisher_exact
    oddsr, fisher_p = fisher_exact(table, alternative='greater')
    # print(k-1, M, n, N, fisher_p, hypergeom_p)
    return hypergeom_p, fisher_p

#def hypergeom_test(arr1, arr2, label=1):
#    """
#    arr1: array with only 1, 0, -1
#    arr2: array with only 1, 0, -1
#    """
#    from scipy.stats import hypergeom
#    M = len(arr1)
#    N = sum(arr1==label)
#    n = sum(arr2==label)
#    k = sum((arr1==label)*(arr2==label))
#    # more overlaps higher pvalues
#    hypergeom_p = hypergeom.sf(k-1, M, n, N)
#    a, b, c, d = k, n-k, N-k, M-(n+N)+k
#    table = np.array([[a, b], [c, d]])
#    start, end = hypergeom.support(M, n, N)
#    # print(sum(hypergeom.pmf(np.arange(k, end+1), M, n, N)))
#    from scipy.stats import fisher_exact
#    oddsr, fisher_p = fisher_exact(table, alternative='greater')
#    # print(k-1, M, n, N, fisher_p, hypergeom_p)
#    return hypergeom_p, fisher_p, M, N, n, k


def ranksumtest(arr1, arr2, alternative='less'):
    from scipy.stats import ranksums
    w, p = ranksums(arr1, arr2, alternative=alternative)
    return w, p







def read_met_model(path):
    import scipy.io
    mat = scipy.io.loadmat(path)
    mat = mat[[k for k in mat.keys()][-1]]
    keys = mat.dtype.fields.keys()
    arrs = mat.item()
    model_dict = {k:arr for k, arr in zip(keys, arrs)}
    return model_dict




def RobustZScore(arr, lb=10**(-5), ub=30):
    
    ind, = np.nonzero(np.abs(arr)<lb)
    if len(ind)>0:
        arr[ind] = np.zeros(len(arr[ind]))
    
    STD = np.std(arr) if np.std(arr)!=0 else 1
    res = (arr-np.mean(arr))/STD
    ind2, = np.nonzero(np.abs(res)>ub)
    print('1')
    print(np.max(np.abs(res)))
    if len(ind2)>0:
        print('2')
        sign = arr[ind2]/np.abs(arr[ind2])
        print(sign)
        ind3, = np.nonzero(np.abs(res)<=ub)
        tmp = res[ind3]
        res[ind2] = sign*np.max(np.abs(tmp))
        print('3')
        print(np.max(np.abs(res)))
    if len(arr)==len(res):
        return res
