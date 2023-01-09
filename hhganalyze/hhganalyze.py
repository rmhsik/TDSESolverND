import numpy as np
import matplotlib.pyplot as plt

def LoadComplexData(file,**genfromtext_args):
    """
    Load complex data in the C++ format in numpy.
    """
    array_as_strings = np.genfromtxt(file,dtype=str,**genfromtext_args)
    complex_parser = np.vectorize(lambda x: complex(*eval(x)))
    return complex_parser(array_as_strings)


class HHGSpectrum:
    def __init__(self, parameters):
        self.path = parameters['path']
        self.w0_ = parameters['w0_']
        self.w0 = parameters['w0']
        self.period = 2*np.pi/self.w0
        self.tmax = parameters['tmax']*self.period
        self.dt = parameters['dt']
        self.Nt = int(self.tmax/self.dt)
        self.t = np.linspace(0,self.tmax,self.Nt)
        self.dw = 2*np.pi/(self.Nt*self.dt)
        self.wmax = self.dw*self.Nt/2
        self.w = np.arange(-self.wmax,self.wmax-self.dw,self.dw)
        
        self.acc = LoadComplexData(self.path)
        
        self.accF = np.fft.fft(self.acc)
        self.accF2 = np.conj(self.accF)*self.accF
        self.accF2 = self.accF2.real
        self.accF2 = np.fft.fftshift(self.accF2)
        
        
    def plot(self):
        fig = plt.figure(figsize=(9,4))
        self.ax = fig.add_subplot(111)
        print(np.ceil(len(self.w)/2))
        self.ax.plot(self.w[int(np.ceil(len(self.w)/2)):]/self.w0_,self.accF2[int(np.ceil(len(self.t)/2)):])
        self.ax.set_yscale('log')
        self.ax.set_xlim(-2,70)
        #ax.set_ylim(1E-11,1E1)
        self.ax.set_xlabel("Harmonic order",fontsize=15)
        self.ax.set_ylabel("Yield (arb. units)",fontsize=15)
        self.ax.tick_params(axis='both', which='major', labelsize=14)
    #     ax.set_title('Envolvente sin$^2$')

        fig.tight_layout()
        
class TimeFreq:
    def __init__(self,data,freq,wstep,w0=0.057):
        self.initdata = data
        self.w = freq
        self.wstep = wstep
        self.dw = np.abs(self.w[1]-self.w[0])
        self.dt = 2*np.pi/(self.w.shape[0]*self.dw*w0)
        self.t = np.arange(0,self.w.shape[0]*self.dt)
        self.period = 2*np.pi/w0
        self.freqs = []
        self.output = []
        self.calculate()

    def mask(self,w0):
        return np.exp(-4*np.log(2.0)*(self.w-w0)**2/self.wstep**2)
    def calculate(self):
        flast = 0.0
        for j in range(int(len(self.w))):
            f = self.w[j]#j*self.dw
            if (np.abs(f-flast)>self.wstep/1.5):
                flast = f
                self.freqs.append(f)
                temp = np.fft.ifft(self.mask(f)*self.initdata)
                #temp = np.fft.fftshift(temp)
                self.output.append(np.conjugate(temp)*(temp))
        self.output = np.fft.fftshift(np.array(self.output),axes=0)