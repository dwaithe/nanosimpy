import numpy as np
import fib4
import time
import thread

"""FCS Bulk Correlation Software

    Copyright (C) 2015  Dominic Waithe

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""


def tttr2xfcs (y,num,NcascStart,NcascEnd, Nsub):
    """autocorr, autotime = tttr2xfcs(y,num,10,20)
     Translation into python of:
     Fast calculation of fluorescence correlation data with asynchronous time-correlated single-photon counting.
     Michael Wahl, Ingo Gregor, Matthias Patting, Jorg Enderlein
     """
 
    dt = np.max(y)-np.min(y)
    y = np.round(y[:],0)
    numshape = num.shape[0]
     
    autotime = np.zeros(((NcascEnd+1)*(Nsub+1),1));
    auto = np.zeros(((NcascEnd+1)*(Nsub+1), num.shape[1], num.shape[1])).astype(np.float64)
    shift = float(0)
    delta = float(1)
    
    
    
    for j in range(0,NcascEnd):
        
        #Finds the unique photon times and their indices. The division of 'y' by '2' each cycle makes this more likely.
        
        y,k1 = np.unique(y,1)
        k1shape = k1.shape[0]
        
        #Sums up the photon times in each bin.
        cs =np.cumsum(num,0).T
        
        #Prepares difference array so starts with zero.
        diffArr1 = np.zeros(( k1shape+1));
        diffArr2 = np.zeros(( k1shape+1));
        
        #Takes the cumulative sum of the unique photon arrivals
        diffArr1[1:] = cs[0,k1].reshape(-1)
        diffArr2[1:] = cs[1,k1].reshape(-1)
        
        #del k1
        #del cs
        num =np.zeros((y.shape[0],2))
        

        
        #Finds the total photons in each bin. and represents as count.
        #This is achieved because we have the indices of each unique time photon and cumulative total at each point.
        num[:,0] = np.diff(diffArr1)
        num[:,1] = np.diff(diffArr2)
        #diffArr1 = [];
        #diffArr2 = [];
        
        for k in range(0,Nsub):
            shift = shift + delta
            lag = np.round(shift/delta,0)
    
            
            #Allows the script to be sped up.
            if j >= NcascStart:
                

                #Old method
                #i1= np.in1d(y,y+lag,assume_unique=True)
                #i2= np.in1d(y+lag,y,assume_unique=True)
                
                #New method, cython
                i1,i2 = fib4.dividAndConquer(y, y+lag,y.shape[0]+1)

                #If the weights (num) are one as in the first Ncasc round, then the correlation is equal to np.sum(i1)
                i1 = i1.astype(np.bool);
                i2 = i2.astype(np.bool);

                #Now we want to weight each photon corectly.
                #Faster dot product method, faster than converting to matrix.
                jin = np.dot((num[i1,:]).T,num[i2,:])/delta
                   
                auto[(k+(j)*Nsub),:,:] = jin
            
            autotime[k+(j)*Nsub] =shift;
        
        #Equivalent to matlab round when numbers are %.5
        y = np.ceil(np.array(0.5*y))
        delta = 2*delta
    
    for j in range(0, auto.shape[0]):
        auto[j,:,:] = auto[j,:,:]*dt/(dt-autotime[j])
    autotime = autotime/1000000


    #Removes the trailing zeros.
    auto = auto[autotime.reshape(-1) != 0,:,:]
    autotime = autotime[autotime != 0]
    
    return auto, autotime


def delayTime2bin(dTimeArr, chanArr, chanNum, winInt):
    
    decayTime = np.array(dTimeArr)
    #This is the point and which each channel is identified.
    decayTimeCh =decayTime[chanArr == chanNum] 
    
    #Find the first and last entry
    firstDecayTime = 0;#np.min(decayTimeCh).astype(np.int32)
    tempLastDecayTime = np.max(decayTimeCh).astype(np.int32)
    
    #We floor this as the last bin is always incomplete and so we discard photons.
    numBins = np.floor((tempLastDecayTime-firstDecayTime)/winInt)
    lastDecayTime = numBins*winInt
    

    bins = np.linspace(firstDecayTime,lastDecayTime, int(numBins)+1)
    
    
    photonsInBin, jnk = np.histogram(decayTimeCh, bins)

    #bins are valued as half their span.
    decayScale = bins[:-1]+(winInt/2)

    #decayScale =  np.arange(0,decayTimeCh.shape[0])
   
    
    

    return list(photonsInBin), list(decayScale)
import numpy as np
import warnings

""" 
    A multiple-tau algorithm for Python 2.7 and 3.x.
    
    Copyright (c) 2014 Paul Muller

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

      1. Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.
       
      2. Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer in
         the documentation and/or other materials provided with the
         distribution.

      3. Neither the name of multipletau nor the names of its contributors
         may be used to endorse or promote products derived from this
         software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL INFRAE OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""



def autocorrelate(a, m=16, deltat=1, normalize=False,
                  copy=True, dtype=None):
    """ 
    Autocorrelation of a 1-dimensional sequence on a log2-scale.
    
    This computes the correlation according to 
    :py:func:`numpy.correlate` for positive :math:`k`  on a base 2
    logarithmic scale.

        numpy.correlate(a, a, mode="full")[len(a)-1:]  

        :math:`z_k = \Sigma_n a_n a_{n+k}`


    Note that only the correlation in the positive direction is
    computed.

    Parameters
    ----------
    a : array_like
        input sequence of real numbers
    m : even integer
        defines the number of points on one level, must be an
        even integer
    deltat : float
        distance between bins
    normalize : bool
        normalize the result to the square of the average input
        signal and the factor `M-k`.
    copy : bool
        copy input array, set to False to save memory
    dtype : dtype, optional
        The type of the returned array and of the accumulator in 
        which the elements are summed.  By default, the dtype of 
        `a` is used.


    Returns
    -------
    autocorrelation : ndarray
        Nx2 array containing lag time and autocorrelation


    Notes
    -----
    The algorithm computes the correlation with the convention of the
    curve decaying to zero.
    
    For experiments like e.g. fluorescence correlation spectroscopy,
    the signal can be normalized to `M-k` by invoking:
    
           normalize = True

    For emulating the numpy.correlate behavior on a logarithmic
    scale (default behavior) use:

           normalize = False


    Examples
    --------
    >>> from numpy import dtype
    >>> from multipletau import autocorrelate
    >>> autocorrelate(range(42), m=2, dtype=dtype(float))
    array([[  1.00000000e+00,   2.29600000e+04],
           [  2.00000000e+00,   2.21000000e+04],
           [  4.00000000e+00,   2.03775000e+04],
           [  8.00000000e+00,   1.50612000e+04]])

    """
    traceavg = np.average(a)
    if normalize and traceavg == 0:
        raise ZeroDivisionError("Normalization not possible. "+
                     "The average of the input *binned_array* is zero.")

    trace = np.array(a, dtype=dtype, copy=copy)
    dtype = trace.dtype

    if dtype.kind in ["b","i","u"]:
        warnings.warn("Converting input data type ({}) to float.".
                      format(dtype))
        dtype = np.dtype(float)
        trace = np.array(a, dtype=dtype, copy=copy)
    
    # Complex data
    if dtype.kind == "c":
        raise NotImplementedError(
              "Please use `multipletau.correlate` for complex data.")

    
    # Check parameters
    if np.around(m/2) != m/2:
        mold = 1*m
        m = int((np.around(m/2)+1) * 2)
        warnings.warn("Invalid value of m={}. Using m={} instead"
                      .format(mold,m))
    else:
        m = int(m)

    N = N0 = len(trace)
    
    # Find out the length of the correlation function.
    # The integer k defines how many times we can average over
    # two neighboring array elements in order to obtain an array of
    # length just larger than m.
    
    k = int(np.floor(np.log2(N/m)))


    # In the base2 multiple-tau scheme, the length of the correlation
    # array is (only taking into account values that are computed from
    # traces that are just larger than m):  
    lenG = np.int(np.floor(m + k*m/2))

    G = np.zeros((lenG, 2), dtype=dtype)

    normstat = np.zeros(lenG, dtype=dtype)
    normnump = np.zeros(lenG, dtype=dtype)
    
    # We use the fluctuation of the signal around the mean
    if normalize:
        trace -= traceavg
    if N < 2*m:
        # Otherwise the following for-loop will fail:
        raise ValueError("len(binned_array) must be larger than 2m.")
    ## Calculate autocorrelation function for first m bins
    # Discrete convolution of m elements
    for n in range(1,m+1):
        G[n-1,0] = deltat * n
        # This is the computationally intensive step
        G[n-1,1] = np.sum(trace[:N-n]*trace[n:], dtype=dtype)
        normstat[n-1] = N-n
        normnump[n-1] = N
    # Now that we calculated the first m elements of G, let us
    # go on with the next m/2 elements.
    # Check if len(trace) is even:
    if N%2 == 1:
        N -= 1
    # Add up every second element
    trace = (trace[:N:2]+trace[1:N+1:2])/2
    N /= 2
    ldx =10000
    ## Start iteration for each m/2 values
    for step in range(1,k+1):
        ## Get the next m/2 values via correlation of the trace
        for n in range(1,int(m/2)+1):
            idx = int(m + n - 1 + (step-1)*m/2)
            if len(trace[:N-(n+m/2)]) == 0:
                # This is a shortcut that stops the iteration once the
                # length of the trace is too small to compute a corre- 
                # lation. The actual length of the correlation function 
                # does not only depend on k - We also must be able to 
                # perform the sum with repect to k for all elements.
                # For small N, the sum over zero elements would be
                # computed here.
                #
                # One could make this for loop go up to maxval, where
                #   maxval1 = int(m/2)
                #   maxval2 = int(N-m/2-1)
                #   maxval = min(maxval1, maxval2)
                # However, we then would also need to find out which 
                # element in G is the last element...
                
                #I modified here because I want the sequence to be predictable in size before it is written.
                #G = G[:idx-1]
                #Finds the first value of idx in this part of loop.
                if ldx >= idx:
                    ldx = idx
                #Keeps indexing.
                G[idx,0] = deltat * (n+m/2) * 2**step
                #Freezes the normalisation.
                normstat = normstat[:ldx-1]
                normnump = normnump[:ldx-1]
                # Note that this break only breaks out of the current
                # for loop. However, we are already in the last loop
                # of the step-for-loop. That is because we calculated
                # k in advance.
                
            else:
                G[idx,0] = deltat * (n+m/2) * 2**step
                # This is the computationally intensive step
                G[idx,1] = np.sum(trace[:N-(n+m/2)]*trace[(n+m/2):],
                                  dtype=dtype)
                normstat[idx] = N-(n+m/2)
                normnump[idx] = N
        # Check if len(trace) is even:
        if N%2 == 1:
            N -= 1
        # Add up every second element
        trace = (trace[:N:2]+trace[1:N+1:2])/2
        N /= 2

    if normalize:
        #Added ldx-1 to insure only normalises valid sequences.
        G[:ldx-1,1] /= traceavg**2 * normstat[:ldx-1]
    else:
        G[:,1] *= N0/normnump 
    
    return G



def correlate(a, v, m=16, deltat=1, normalize=False,
              copy=True, dtype=None):
    """ 
    Cross-correlation of two 1-dimensional sequences
    on a log2-scale.
    
    This computes the cross-correlation according to
    :py:func:`numpy.correlate` for positive :math:`k`  on a base 2
    logarithmic scale.
    
        numpy.correlate(a, v, mode="full")[len(a)-1:]
        
        :math:`z_k = \Sigma_n a_n v_{n+k}`
    
    Note that only the correlation
    in the positive direction is computed.
    
    
    Parameters
    ----------
    a, v : array_like
        input sequences with equal length
    m : even integer
        defines the number of points on one level, must be an
        even integer
    deltat : float
        distance between bins
    normalize : bool
        normalize the result to the square of the average input
        signal and the factor `M-k`.
    copy : bool
        copy input array, set to False to save memory
    dtype : dtype, optional
        The type of the returned array and of the accumulator in 
        which the elements are summed.  By default, the dtype of 
        `a` is used.


    Returns
    -------
    crosscorrelation : ndarray
        Nx2 array containing lag time and cross-correlation
        
    
    Notes
    -----
    The algorithm computes the correlation with the convention of the
    curve decaying to zero.
    
    For experiments like e.g. fluorescence correlation spectroscopy,
    the signal can be normalized to `M-k` by invoking:
    
           normalize = True

    For emulating the numpy.correlate behavior on a logarithmic
    scale (default behavior) use:

           normalize = False


    Examples
    --------
    >>> from numpy import dtype
    >>> from multipletau import correlate
    >>> correlate(range(42), range(1,43), m=2, dtype=dtype(float))
    array([[  1.00000000e+00,   2.38210000e+04],
           [  2.00000000e+00,   2.29600000e+04],
           [  4.00000000e+00,   2.12325000e+04],
           [  8.00000000e+00,   1.58508000e+04]])

    """
    ## See `autocorrelation` for better documented code.
    traceavg1 = np.average(v)
    traceavg2 = np.average(a)
    if normalize and traceavg1*traceavg2 == 0:
        raise ZeroDivisionError("Normalization not possible. "+
                     "The average of the input *binned_array* is zero.")

    trace1 = np.array(v, dtype=dtype, copy=copy)
    dtype = trace1.dtype

    if dtype.kind in ["b","i","u"]:
        warnings.warn(
              "Converting input data type ({}) to float.".format(dtype))
        dtype = np.dtype(float)
        trace1 = np.array(v, dtype=dtype, copy=copy)
    
    # Prevent traces from overwriting each other
    if a is v:
        # Force copying trace 2
        copy = True
        
    trace2 = np.array(a, dtype=dtype, copy=copy)

    # Complex data
    if dtype.kind == "c":
        trace1.imag *= -1
   
    # Check parameters
    if np.around(m/2) != m/2:
        mold = 1*m
        m = int((np.around(m/2)+1) * 2)
        warnings.warn("Invalid value of m={}. Using m={} instead"
                      .format(mold,m))
    else:
        m = int(m)
        
    if len(a) != len(v):
        raise ValueError("Input arrays must be of equal length.")
        
    N = N0 = len(trace1)
    # Find out the length of the correlation function.
    # The integer k defines how many times we can average over
    # two neighboring array elements in order to obtain an array of
    # length just larger than m.
    k = int(np.floor(np.log2(N/m)))
    
    # In the base2 multiple-tau scheme, the length of the correlation
    # array is (only taking into account values that are computed from
    # traces that are just larger than m):   
    lenG = np.int(np.floor(m + k*m/2))
        
    G = np.zeros((lenG, 2), dtype=dtype)
    normstat = np.zeros(lenG, dtype=dtype)
    normnump = np.zeros(lenG, dtype=dtype)
    
    # We use the fluctuation of the signal around the mean
    if normalize:
        trace1 -= traceavg1
        trace2 -= traceavg2
    if N < 2*m:
        # Otherwise the following for-loop will fail:
        raise ValueError("len(binned_array) must be larger than 2m.")
    # Calculate autocorrelation function for first m bins
    for n in range(1,m+1):
        G[n-1,0] = deltat * n
        G[n-1,1] = np.sum(trace1[:N-n]*trace2[n:])
        normstat[n-1] = N-n
        normnump[n-1] = N
    # Check if len(trace) is even:
    if N%2 == 1:
        N -= 1
    # Add up every second element
    trace1 = (trace1[:N:2]+trace1[1:N+1:2])/2
    trace2 = (trace2[:N:2]+trace2[1:N+1:2])/2
    N /= 2

    for step in range(1,k+1):
        # Get the next m/2 values of the trace
        for n in range(1,int(m/2)+1):
            idx = int(m + n - 1 + (step-1)*m/2)
            if len(trace1[:N-(n+m/2)]) == 0:
                # Abort
                G = G[:idx-1]
                normstat = normstat[:idx-1]
                normnump = normnump[:idx-1]
                break
            else:
                G[idx,0] = deltat * (n+m/2) * 2**step
                G[idx,1] = np.sum(trace1[:N-(n+m/2)]*trace2[(n+m/2):])
                normstat[idx] = N-(n+m/2)
                normnump[idx] = N

        # Check if len(trace) is evn:
        if N%2 == 1:
            N -= 1
        # Add up every second element
        trace1 = (trace1[:N:2]+trace1[1:N+1:2])/2
        trace2 = (trace2[:N:2]+trace2[1:N+1:2])/2
        N /= 2

    if normalize:
        G[:,1] /= traceavg1*traceavg2 * normstat
    else:
        G[:,1] *= N0/normnump 
    
    return G
