#=======================================
# DCT-1 : Discrete cosine transform
#---------------------------------------
# Coefficients are same as defined by 
# 6.14 in the book
#
#=======================================
from pylab import*

def dct1(x):

    M      = len(x)
    y      = vstack((x, flipud(x[1:M-1])))
    y_ft   = fft(y, axis = 0)
    y      = real(y_ft[0:M]/(M-1))
    y[0,:] = y[0,:]/2.
    y[M-1,:] = y[M-1,:]/2.

    return y

#=======================================
# DST-1 : Discrete sine transform
#---------------------------------------
# Coefficients are same as defined by 
# 6.17 in the book
#
#=======================================
def dst1(x):

    M, N = shape(x)
    x    = x[1:M-1,:]
    M, N = shape(x)
    y    = vstack((zeros([1,N]), x, zeros([1,N]), -flipud(x)))
    yy   = fft(y, axis = 0)
    y    = real(yy[1:M+1,:]/(-1j*(M+1)))

    y    = vstack((zeros([1,N]), y, zeros([1,N])))

    return y
