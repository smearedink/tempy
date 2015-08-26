import numpy as _np

def find_freq_clusters(freqs):
    # first make a histogram
    minf, maxf = freqs.min(), freqs.max()
    maxbins = 8  # related to the max colors defined...
    df = 4.0 # MHz
    if ((maxf - minf) < df):  # Only a single freq to our resolution
        return [[0.0, 'inf']]
    numbins = int((maxf - minf) / df) + 2
    lobound = minf - 0.5 * df
    hibound = lobound + numbins * df
    hist, edges = _np.histogram(freqs, numbins, [lobound, hibound])
    # Now choose the maxbins biggest bins where there are TOAs
    hibins = hist.argsort()[::-1]
    hibins = hibins[hist[hibins] > 0]
    if len(hibins) > maxbins:
        hibins = hibins[:maxbins]
    ctrs = edges[hibins] + 0.5 * df
    ctrs.sort()
    # and use these as starting points for kmeans
    kmeans, indices = kmeans2(freqs, ctrs)
    if len(kmeans)==1:
        return [[0.0, 'inf']]
    elif len(kmeans)==2:
        return [[0.0, kmeans.mean()], [kmeans.mean(), 'inf']]
    else:
        freqbands = [[0.0, kmeans[0:2].mean()]]
        for ii in range(len(kmeans)-2):
            freqbands.append([kmeans[ii:ii+2].mean(), kmeans[ii+1:ii+3].mean()])
        freqbands.append([kmeans[-2:].mean(), 'inf'])
        return freqbands

### This un2str code is taken near-verbatim from Lemming's reply at
### http://stackoverflow.com/questions/6671053/python-pretty-print-errorbars
def un2str(x, xe, precision=1):
    """pretty print nominal value and uncertainty

    x  - nominal value
    xe - uncertainty
    precision - number of significant digits in uncertainty

    returns shortest string representation of `x +- xe` either as
        x.xx(ee)e+xx
    or as
        xxx.xx(ee)"""
    # base 10 exponents
    x_exp = int(_np.floor(_np.log10(abs(x))))
    xe_exp = int(_np.floor(_np.log10(xe)))

    # uncertainty
    un_exp = xe_exp-precision+1
    un_int = round(xe*10**(-un_exp))

    # nominal value
    no_exp = un_exp
    no_int = round(x*10**(-no_exp))

    # format - nom(unc)exp
    fieldw = x_exp - no_exp
    fmt = '%%.%df' % fieldw
    result1 = (fmt + '(%.0f)e%d') % (no_int*10**(-fieldw), un_int, x_exp)

    # format - nom(unc)
    fieldw = max(0, -no_exp)
    fmt = '%%.%df' % fieldw
    result2 = (fmt + '(%.0f)') % (no_int*10**no_exp, un_int*10**max(0, un_exp))

    # return shortest representation
    if len(result2) <= len(result1):
        return result2
    else:
        return result1
