from scipy import fftpack
import numpy as np
from skimage import filters

def GaussianMaskFit2(im,coo,s,optLoc=1,bgSub=2,winSize=13,convDelta=.01,nbIter=20):
    """Applies the algorithm from [Thompson et al. (2002) PNAS, 82:2775].
    Parameters:
    - im: a numpy array with the image
    - coo: approximate coordinates (in pixels) of the spot to localize and measure. Note, the coordinates are x,y!
    - s: width of the PSF in pixels
    - optLoc: If 1, applied the iterative localization refinement algorithm, starting with the coordinates provided in coo. If 0, only measures the spot intensity at the coordinates provided in coo.
    - bgSub: 0 -> no background subtraction. 1 -> constant background subtraction. 2 -> tilted plane background subtraction.
    - winSize: Size of the window (in pixels) around the position in coo, used for the iterative localization and for the background subtraction.
    - convDelta: cutoff to determine convergence, i.e. the distance (in pixels) between two iterations
    - nbIter: the maximal number of iterations.

    Returns
    - the intensity value of the spot.
    - the corrdinates of the spot.

    If convergence is not found after nbIter iterations, return 0 for both intensity value and coordinates.
    """
    coo=np.array(coo).astype('float')
    for i in range(nbIter):
        if not np.prod(coo-winSize/2.>=0)*np.prod(coo+winSize/2.<=im.shape[::-1]): return 0., np.r_[0., 0.], 0.
        winOrig=(coo-int(winSize)//2).astype(int)
        ix, iy=np.meshgrid(winOrig[0]+np.r_[:winSize], winOrig[1]+np.r_[:winSize])
        N=np.exp(-(ix-coo[0])**2/(2*s**2)-(iy-coo[1])**2/(2*s**2))/(2*np.pi*s**2)
        S=im[:, winOrig[0]:winOrig[0]+winSize][winOrig[1]:winOrig[1]+winSize]*1.
        if bgSub==2:
            xy=np.r_[:2*winSize]%winSize-(winSize-1)/2.
            bgx=np.polyfit(xy, np.r_[S[0], S[-1]], 1); S=(S-xy[:winSize]*bgx[0]).T
            bgy=np.polyfit(xy, np.r_[S[0], S[-1]], 1); S=(S-xy[:winSize]*bgy[0]).T
            bg=np.mean([S[0], S[-1], S[:, 0], S[:, -1],]); S-=bg
            bg=np.r_[bg, bgx[0], bgy[0]]
        if bgSub==1:
            bg=np.mean([S[0], S[-1], S[:, 0], S[:, -1],]); S-=bg
        S=S.clip(0) # Prevent negative values !!!!
        if optLoc:
            SN=S*N; ncoo=np.r_[np.sum(ix*SN), np.sum(iy*SN)]/np.sum(SN)
            #ncoo=ncoo+ncoo-coo # Extrapolation
            if abs(coo-ncoo).max()<convDelta: return np.sum(SN)/np.sum(N**2), coo, bg
            else: coo=ncoo
        else: return np.sum(S*N)/np.sum(N**2), coo, bg
    return 0., np.r_[0., 0.], 0.

def GaussianMaskFit3D(im,coo,s,sZ=None,optLoc=1,bgSub=2,winSize=None,winSizeZ=None,convDelta=.05,nbIter=20):
  """Applies the algorithm from [Thompson et al. (2002) PNAS, 82:2775] adapted to 3D images.
Parameters:
- im: a numpy array with the image
- coo: approximate z,y,x coordinates (in pixels) of the spot to localize and measure.
- s: width of the PSF in x,y in pixels
- sZ: width of the PSF in z in pixels. Defaults to the same value as s.
- optLoc: If 1, applied the iterative localization refinement algorithm, starting with the coordinates provided in coo. If 0, only measures the spot intensity at the coordinates provided in coo.
- bgSub: 0 -> no background subtraction. 1 -> constant background subtraction. 2 -> tilted plane background subtraction.
- winSize: Size of the x,y window (in pixels) around the position in coo, used for the iterative localization and for the background subtraction.
- winSizeZ: Same as winSize, in the z dimension.
- convDelta: cutoff to determine convergence, i.e. the distance (in pixels) between two iterations
- nbIter: the maximal number of iterations.

Returns
- the intensity value of the spot.
- the coordinates of the spot (z,y,x).
- the background level:
   - either a constant value if bgSub=1
   - or [offset, tilt in z, tilt in y, tilt in x] if bgSub=2

If convergence is not found after nbIter iterations, return 0 for both intensity value and coordinates.
"""
  coo=np.array(coo).astype('float')
  if sZ==None: sZ=s
  if winSize ==None: winSize =int(np.ceil(s*8./2))*2+1
  if winSizeZ==None: winSizeZ=int(np.ceil(sZ*4./2))*2+1
  for i in range(nbIter):
    if not (winSizeZ/2.<=coo[0]<=im.shape[0]-winSizeZ/2.)*np.prod([winSize/2.<=coo[j]<=im.shape[j]-winSize/2. for j in [1, 2]]):
      return 0., np.r_[0., 0., 0.], 0.
    winOrig=np.r_[coo[0]-int(winSizeZ//2), coo[1:]-int(winSize//2)].astype(int)
    iy, iz, ix=np.meshgrid(winOrig[1]+np.r_[:winSize], winOrig[0]+np.r_[:winSizeZ], winOrig[2]+np.r_[:winSize])
    N=np.exp(-(iz-coo[0])**2/(2*sZ**2)-(iy-coo[1])**2/(2*s**2)-(ix-coo[2])**2/(2*s**2))/((2*np.pi)**1.5*s*s*sZ)
    S=im[winOrig[0]:winOrig[0]+winSizeZ][:, winOrig[1]:winOrig[1]+winSize][:,:, winOrig[2]:winOrig[2]+winSize]*1.
    if bgSub==2:
      cxy=np.r_[:winSize]-(winSize-1)/2.
      cz=np.r_[:winSizeZ]-(winSizeZ-1)/2.
      bgx=np.polyfit(cxy, np.mean(np.r_[S[:, 0], S[:, -1]], 0), 1)[0]
      bgy=np.polyfit(cxy, np.mean(np.r_[S[:,:, 0], S[:,:, -1]], 0), 1)[0]
      bgz=np.polyfit(cz, np.mean(np.c_[S[:, 0], S[:, -1], S[:, 1:-1, 0], S[:, 1:-1, -1]], 1), 1)[0]
      S=np.rollaxis(np.rollaxis(np.rollaxis(S-cxy*bgx, 2)-cxy*bgy, 2)-cz*bgz, 2)
      bg=np.mean([S[:, 0], S[:, -1], S[:,:, 0], S[:,:, -1],]); S-=bg
      bg=np.r_[bg, bgz, bgy, bgx]
    elif bgSub==1:
      bg=np.mean([S[:, 0], S[:, -1], S[:,:, 0], S[:,:, -1],]); S-=bg
    else:
      bg=0
    #S=S.clip(0) # Prevent negative values !!!!
    if optLoc:
      SN=S*N; ncoo=np.r_[np.sum(iz*SN), np.sum(iy*SN), np.sum(ix*SN)]/np.sum(SN)
      #ncoo+=ncoo-coo # Extrapolation of localization step !!!!
      #ncoo+=(ncoo-coo)*.7 # Extrapolation of localization step !!!!
      #print(i,ncoo,abs(coo-ncoo).max())
      if abs(coo-ncoo).max()<convDelta: return np.sum(SN)/np.sum(N**2), coo, bg
      else: coo=ncoo
    else: return np.sum(S*N)/np.sum(N**2), coo, bg
  return 0., np.r_[0., 0., 0.], 0.


def gauss3(I, p, s, sz, bg, S):
    """ Create an image of a Gaussian peak
        I:  peak intensity
        p:  peak location (z, y, x)
        s:  peak width (x & y)
        sz: peak width (z)
        bg: tilt (offset, z, y, x) or just offset only
        S:  image size (z, y, x)

        returns: 3D array

        wp@tl20200713
    """
    yv, zv, xv = np.meshgrid(np.arange(S[1], dtype='float')-p[1],
                             np.arange(S[0], dtype='float')-p[0],
                             np.arange(S[2], dtype='float')-p[2])
    im = I * np.exp(-(xv**2/2/s**2 + yv**2/2/s**2 + zv**2/2/sz**2)) / ((2*np.pi)**1.5*s**2*sz)
    try:
        return im + bg[3]*xv + bg[2]*yv + bg[1]*zv + bg[0]
    except Exception:
        return im + bg


def GaussianMaskFit3D_err(im, I, coo, bg, s, sZ=None, winSize=None, winSizeZ=None):
    """ Calculate R2 and confidence intervals (SE) on the parameters from GaussianMaskFit3D
        im:       a numpy array with the image
        I:        fitted intensity of the localized spot.
        coo:      fitted z,y,x coordinates (in pixels) of the localized spot.
        bg:       fitted background of the localized spot.
        s:        width of the PSF in x,y in pixels
        sZ:       width of the PSF in z in pixels. Defaults to the same value as s.
        winSize:  Size of the x,y window (in pixels) around the position in coo, used for the iterative localization and for the background subtraction.
        winSizeZ: Same as winSize, in the z dimension.

        returns: R2, dI, dcoo, dbg

        wp@tl20200713
    """
    sZ = sZ or s
    winSize = winSize or int(np.ceil(s * 8. / 2)) * 2 + 1
    winSizeZ = winSizeZ or int(np.ceil(sZ * 4. / 2)) * 2 + 1
    winOrig = np.r_[coo[0] - int(winSizeZ / 2), coo[1:] - int(winSize / 2)].astype(int)
    coo = coo - winOrig
    S = im[winOrig[0]:winOrig[0] + winSizeZ][:, winOrig[1]:winOrig[1] + winSize][:,:, winOrig[2]:winOrig[2] + winSize] * 1.
    fun = lambda q: gauss3(q[0], q[1:4], s, sZ, q[4:], S.shape)
    q = [I]
    q.extend(coo)
    try:
        q.extend(bg)
    except Exception:
        q.append(bg)
    X2, dq, R2 = fminerr(fun, q, S)

    return R2, dq[0], dq[1:4], dq[4:]


def fminerr(fun, a, y, dy=None, diffstep=1e-6):
    """ Error estimation of a fit

        Inputs:
        fun: function which was fitted to data
        a:   function parameters
        y:   ydata
        dy:  errors on ydata

        Outputs:
        chisq: Chi^2
        da:    error estimates of the function parameters
        R2:    R^2

        Example:
        x = np.array((-3,-1,2,4,5))
        a = np.array((2,-3))
        y = (15,0,5,30,50)
        fun = lambda a: a[0]*x**2+a[1]
        chisq,dp,R2 = fminerr(fun,p,y)

        adjusted from Matlab version by Thomas Schmidt, Leiden University
        wp@tl2020
    """
    eps = np.spacing(1)
    a = np.array(a).flatten()
    y = np.array(y).flatten()
    if dy is None:
        dy = np.ones(np.shape(y))
    else:
        dy = np.array(dy).flatten()
    nData = np.size(y)
    nPar = np.size(a)
    dy = 1 / (dy + eps)
    f0 = np.array(fun(a)).flatten()
    chisq = np.sum(((f0 - y) * dy) ** 2) / (nData - nPar)

    # calculate R^2
    sstot = np.sum((y - np.nanmean(y)) ** 2)
    ssres = np.sum((y - f0) ** 2)
    R2 = 1 - ssres / sstot

    # calculate derivatives
    deriv = np.zeros((nData, nPar))
    for i in range(nPar):
        ah = a.copy()
        ah[i] = a[i] * (1 + diffstep) + eps
        f = np.array(fun(ah)).flatten()
        deriv[:, i] = (f - f0) / (ah[i] - a[i]) * dy

    hesse = np.matmul(deriv.T, deriv)

    if np.linalg.matrix_rank(hesse) == np.shape(hesse)[0]:
        da = np.sqrt(chisq * np.diag(np.linalg.inv(hesse)))
    else:
        try:
            da = np.sqrt(chisq * np.diag(np.linalg.pinv(hesse)))
        except Exception:
            da = np.zeros(a.shape)
        # da = np.full(np.shape(a),np.nan)
        # print('Hessian not invertible, size: {}, rank: {}'.format(np.shape(hesse)[0],np.linalg.matrix_rank(hesse)))
    return chisq, da, R2


sHS=fftpack.fftshift # Swap half-spaces. sHS(matrix[, axes]). axes=all by default
def hS(m,axes=None):
    if axes==None: axes=range(np.ndim(m))
    elif isinstance(axes, int): axes=[axes]
    elif axes==[]: return m
    return hS(m.swapaxes(0, axes[-1])[:m.shape[axes[-1]]/2].swapaxes(0, axes[-1]), axes[:-1])

def sHSM(m,axes=None):
    if axes==None: axes=range(np.ndim(m))
    elif isinstance(axes, int): axes=[axes]
    m=m.swapaxes(0, axes[0]); max=m[1]+m[-1]; m=(m+max/2)%max-max/2; m=m.swapaxes(0, axes[0])
    return sHS(m, axes)

def bpass(im,r1=1.,r2=1.7):
    ker1x=np.exp(-(sHS(sHSM(np.r_[:im.shape[1]]))/r1)**2/2); ker1x/=np.sum(ker1x); fker1x=fftpack.fft(ker1x)
    ker1y=np.exp(-(sHS(sHSM(np.r_[:im.shape[0]]))/r1)**2/2); ker1y/=np.sum(ker1y); fker1y=fftpack.fft(ker1y)
    ker2x=np.exp(-(sHS(sHSM(np.r_[:im.shape[1]]))/r2)**2/2); ker2x/=np.sum(ker2x); fker2x=fftpack.fft(ker2x)
    ker2y=np.exp(-(sHS(sHSM(np.r_[:im.shape[0]]))/r2)**2/2); ker2y/=np.sum(ker2y); fker2y=fftpack.fft(ker2y)
    fim=fftpack.fftn(im)
    return fftpack.ifftn((fim*fker1x).T*fker1y-(fim*fker2x).T*fker2y).real.T

def bpass3D(im,r1=1.,r2=1.7,rz1=1.,rz2=1.7,zMirror=False):
    psfPx = r2
    psfPxZ = rz2
    return filters.gaussian(im, 1., mode='mirror')-filters.gaussian(im, np.r_[psfPxZ, psfPx, psfPx], mode='mirror')


