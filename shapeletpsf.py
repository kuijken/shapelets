import numpy as np
import matplotlib.pyplot as plt

# could add fitting routines here too; so far only reading in and decoding

class shapeletmap:
    ''' Class representing a shapelet psf map, which allows the PSF at any
        position in an image to be expressed as a shapelet sum
    '''
    def __init__(self,mapfile=None):
        ''' Set up new shapeletmap; optionally load coefficients from mapfile'''
        if not mapfile is None:
            self.load(mapfile)
        return

    def load(self,mapfile):
        ''' Read a shapelet map from file and create shapeletmap object'''
        f=open(mapfile,'r')
        nsh,beta=np.loadtxt(f,max_rows=1)
        nsh=int(nsh)
        n,xmax,ymax=np.loadtxt(f,max_rows=1)
        n,xmax,ymax=int(n),int(xmax),int(ymax)
        a=np.loadtxt(f)
        if a.ndim<2: a=a.reshape(((nsh+1)*(nsh+2)//2,(n+1)*(n+2)//2))
        if len(a) != (nsh+1)*(nsh+2)//2:
            raise Exception('map file does not contain right number of rows')
        if len(a[0]) != (n+1)*(n+2)//2:
            raise Exception('map file does not contain right number of columns')
        self.nsh=nsh
        self.beta=beta
        self.n=n
        self.xmax=xmax
        self.ymax=ymax
        self.mapcoefs=a
        print('Loaded (nsh='+str(nsh)+', n='+str(n)+') shapelet map from',mapfile)
        return

    def atxy(self,x,y):
        ''' Evaluate the map at position (x,y), return shapeletseries object'''
        t=2*x/self.xmax-1
        u=2*y/self.ymax-1
        sh=shapeletseries(nsh=self.nsh,beta=self.beta)
        k=0
        for i in range(self.n+1):
            for m in range(i+1):
                l=i-m
                sh.coefs+=self.mapcoefs[:,k] * t**l * u**m
                k+=1
        sh.name='PSF at position %.1f,%.1f' % (x,y)
        return sh
    
class shapeletseries:
    ''' Class representing a sum of shapelets that describes a source.
        Initialise by specifying order nsh, scale radius beta, and coeffs array
    '''
    def __init__(self,nsh=10,beta=1.,coefs=None):
        ''' Create a shapelet series of given order, scale radius, coefs'''
        self.nsh=nsh
        self.beta=beta
        self.name=''
        if coefs is None:
            self.coefs=np.zeros( (nsh+1)*(nsh+2)//2 )
        else:
            if len(coefs) == (nsh+1)*(nsh+2)//2:
                self.coefs=coefs
            else:
                raise Exception('supplied coefficients array has wrong length')
        return

    def toimage(self,npix=25,xc=0.,yc=0.,display=False):
        ''' Convert shapeletseries to a npix x npix numpy array and return it.
            If npix is odd then the central pixel corresponds to position -xc,-yc.
            xc and yc are the pixel coordinates of the PSF peak'''
        z=np.zeros((npix,npix))
        nsh=self.nsh
        b=self.beta
        s=self.coefs
        mpix=(npix-1)/2     #  pixels count from -mpix to +mpix, npix in total
        x,y=np.meshgrid(np.linspace(-xc-mpix,-xc+mpix,npix),
                        np.linspace(-yc-mpix,-yc+mpix,npix))
        expfac=np.exp(-(x*x+y*y)/(2*b*b))/b
        hx=hermites(x/b,nsh)
        hy=hermites(y/b,nsh)
        i=0
        for n in range(nsh+1):
            for l in range(n+1):
                k=n-l
                z+=expfac*hx[k]*hy[l]*s[i]
                i+=1
        if display:
            plt.imshow(z,origin='lower',
                           extent=[-npix/2,npix/2,-npix/2,npix/2])
            plt.title(self.name)
            plt.xlabel('X [pix]')
            plt.ylabel('Y [pix]')
            plt.show()
        return z

def hermites(x,nsh):
        ''' construct Hermite polynomials Hn(x), for orders 0..nsh,
            and normalise so that Hn(x) Hm(x) exp(-x^2) integrate to delta_mn.
            return the Hn(x) as a dictionary indexed by order n'''
        h={}
        h[0]=0*x+1
        h[1]=2*x
        for n in range(1,nsh):
            h[n+1]=2*x*h[n] - 2*n*h[n-1]
        nfac=1
        for n in range(nsh+1):
            h[n]/=(np.pi**0.5*2**n*nfac)**0.5
            nfac*=n+1
        return h
    
