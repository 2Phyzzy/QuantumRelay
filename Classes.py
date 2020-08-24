class Source:
	"""General PPLN source
	Values:
	
	pair_count				The rate at which entangled pairs are created by the source
	heralding_efficiency	The percentage of the Singles that are the Signal
	frequency1				The central frequency of output 1
	frequency2				The central frequency of output 2
	frequency_bandwidth	The frequency bandwidth of the two channels
	pulse_width				The time-like width of the pulse
	pulse_rate				The rate at which pulses occur
	
	Definitions:
	
	None
	"""
	
	def __init__(self, pair_count, frequency1_THz=193.4, frequency2_THz=193.4,
				 frequency_bandwidth_THz=0.1, pulse_width_ps=80, pulse_rate_MHz=2, singles_fluorescence = 3000,
				 dark_counts = 1000, afterpulsing = 0.15, channel_isolation = 500, DWDM_coupling_loss = 0.35):

		from scipy import interpolate
		import numpy as np
		"""These lists are for the JSA into correct channel entanglement"""
		
		#x=[100.0, 20.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0]
		#y=[0.5003342073340652, 0.500599828621913, 0.5014072421146843, 0.5016530028064845, 0.501991411119732, 0.5024745600881322, 0.503195889629716, 0.5043351171059877, 0.506268974364264, 0.5098804671466501, 0.5176579123564292, 0.5412619347598981, 0.5529911670071052, 0.5725002827235572, 0.6119168150821013, 0.7257876869657096, 0.8553024939862091, 0.95929858779096, 1]
		x = [0.1,35]
		y = [0.4091753018290429,0.49013461820231713]
		f = interpolate.interp1d(x,y, fill_value='extrapolate')
		
		self.frequency1 = frequency1_THz*(10**12)
		if frequency1_THz > 1000:
			print("Please check that frequency1 is in THz")
		self.frequency2 = frequency2_THz*(10**12)
		if frequency2_THz > 1000:
			print("Please check that frequency2 is in THz")
		self.frequency_bandwidth = frequency_bandwidth_THz
		if frequency_bandwidth_THz > 1:
			print("Please check that the Frequency Band Width is in THz")
		self.pulse_width = pulse_width_ps*(10**(-12))
		if pulse_width_ps < 0.01:
			print("Please check that the Pulse Width is in ps")
		self.pulse_rate = pulse_rate_MHz*(10**6)
		if pulse_rate_MHz > 1000:
			print("Please check that the Pulse Spacing in in MHz")
			
		fwhm = pulse_width_ps*(np.sqrt(2*np.log(2))/2)
		
		self.pair_count = pair_count*(1-DWDM_coupling_loss)
		#print(float(f(fwhm)))
		self.heralding_efficiency = float(f(fwhm))
		self.singles_fluorescence = singles_fluorescence
		self.dark_counts = dark_counts
		self.afterpulsing = afterpulsing
		self.channel_isolation = channel_isolation
		self.DWDM_coupling_loss = DWDM_coupling_loss
		noi_he = (pair_count/self.heralding_efficiency) - pair_count
		noi_ap = self.pair_count*self.afterpulsing
		noi_sf = self.singles_fluorescence
		noi_dc = self.dark_counts
		noi_ci = self.channel_isolation
		self.singles_count = (noi_he + noi_ap + noi_sf + noi_dc + noi_ci)*(1-DWDM_coupling_loss)
		





class Sync:
	"""The system that synchronises the pulses in the two sources
		Values:
		
		Source1						The 1st source
		Source2						The 2nd source
		jitter						The time synchronisation Jitter
		coincidence_window			The coincidence window required for photon interaction
		frequency_bandwidth1		The bandwidth of the 1st source
		frequency_bandwidth2		The bandwidth of the 2nd source
		coincidence_window_minimum1	The minimum possible coincidece window for Source 1
		coincidence_window_minimum2	The minimum possible coincidece window for Source 2
		coincidence_length			The minimum possible coincidece length for Source 1
		coincidence_length			The minimum possible coincidece length for Source 2				
		fwhm1						
		fwhm2
		sync_jitter						
		
		Definitions:
		
		CoincidenceOverlap(source_number, plot=False)
				This either calculates the probability of the output photons being within half of the coincidence window of the peak of the pulse 
				or produces 4 lists that are the pulse gaussian x and y and the coincidence window sectiond peak x and y.
		
		ArrivalProbabilityatBSMInput_NoLoss()
				Currently no definition
		"""
	
	def __init__(self,Source1, Link1, Source2, Link2, coincidence_window = 200E-12, sync_jitter = 40E-12):
		import numpy as np
		self.Source1 = Source1
		self.Source2 = Source2
		self.Link1 = Link1
		self.Link2 = Link2
		#self.jitter = jitter
		self.coincidence_window = coincidence_window #s
		self.frequency_bandwidth1 = self.Source1.frequency_bandwidth #1/s
		self.frequency_bandwidth2 = self.Source2.frequency_bandwidth #1/s
		self.coincidence_window_minimum1 = 1/self.frequency_bandwidth1 #s
		self.coincidence_window_minimum2 = 1/self.frequency_bandwidth2 #s
		self.coincidence_length1 = 299792458/self.frequency_bandwidth1 #m
		self.coincidence_length2 = 299792458/self.frequency_bandwidth2 #m
		self.sync_jitter = sync_jitter
		self.fwhm1 = self.Source1.pulse_width*(np.sqrt(2*np.log(2))/2)
		self.fwhm2 = self.Source2.pulse_width*(np.sqrt(2*np.log(2))/2)
		
	def FWHM_withdispersion(self, source_number):
		import numpy as np
		
		if source_number == 1:
			pulse_width = self.Source1.pulse_width
			pulse_width_new = np.sqrt(pulse_width**2 + self.Link1.Dispersion_Jitter()**2)
			
		elif source_number == 2:
			pulse_width = self.Source2.pulse_width
			pulse_width_new = np.sqrt(pulse_width**2 + self.Link2.Dispersion_Jitter()**2)
			
		else:
			pulse_width = self.Source1.pulse_width
			pulse_width_new = np.sqrt(pulse_width**2 + self.Link1.Dispersion_Jitter()**2)
			print("Error: Issue with source_number ",source_number,". Source set to source_number = 1")		
			
		fwhm = pulse_width_new*(np.sqrt(2*np.log(2))/2)
		
		return fwhm

	def CoincidenceOverlap(self, source_number, plot=False, Jitter_offset=0, preciseness=10000,dispersion=True,minmax =1E-9):
		import numpy as np
		
		preciseness=preciseness+1
		
		
		if source_number == 1:
			if dispersion is False:
				fwhm = self.fwhm1
				source = self.Source1
			elif dispersion is True:
				fwhm = self.FWHM_withdispersion(source_number)
				source = self.Source1
			else:
				print("Error: value 'dispersion' must be True or False. 'dispersion' set to True")
				fwhm = self.FWHM_withdispersion(source_number)
				source = self.Source1
				
			
		elif source_number == 2:
			if dispersion is False:
				fwhm = self.fwhm2
				source = self.Source2
			elif dispersion is True:
				fwhm = self.FWHM_withdispersion(source_number)
				source = self.Source2
			else:
				print("Error: value 'dispersion' must be True or False. 'dispersion' set to True")
				fwhm = self.FWHM_withdispersion(source_number)
				source = self.Source2


		else:
			if dispersion is False:
				fwhm = self.fwhm1
				source = self.Source1
			elif dispersion is True:
				fwhm = self.FWHM_withdispersion(source_number)
				source = self.Source1
			else:
				print("Error: value 'dispersion' must be True or False. 'dispersion' set to True")
				fwhm = self.FWHM_withdispersion(source_number)
				source = self.Source1
			print("Error: Issue with source_number ",source_number,". Source set to source_number = 1")		
		
		Gs = []
		ts = []
		for t in np.linspace(-minmax,minmax,preciseness):
			A = -(2*(t-Jitter_offset)/fwhm)**2
			G = (2**A)
			Gs.append(G)
			ts.append(t)
		Gs = np.array(Gs)
		ts = np.array(ts)
		step = 2*minmax/preciseness
		norm = np.sum(Gs*step)
		Gs = Gs/norm
		
		Gs_Co = []
		ts_Co = []
		
		for t in np.linspace(-self.coincidence_window/2.0,self.coincidence_window/2.0,preciseness):
			A = -(2*(t-Jitter_offset)/fwhm)**2
			G = (2**A)
			Gs_Co.append(G)
			ts_Co.append(t)
		Gs_Co = np.array(Gs_Co)
		ts_Co = np.array(ts_Co)
		step = self.coincidence_window/preciseness
		norm = np.sum(Gs_Co*step)
		Gs_Co = Gs_Co/norm
		
		Total         = np.sum(Gs*step)
		Total_Overlap = np.sum(Gs_Co*step)
		
		if plot is True:
			return ts, Gs, ts_Co, Gs_Co
			
		if plot is False:
			return Total, Total_Overlap

	def CrossOver(self,m1,m2,std1,std2):
		import numpy as np
		a = 1/(2*std1**2) - 1/(2*std2**2)
		b = m2/(std2**2) - m1/(std1**2)
		c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
		return np.roots([a,b,c])
	
	def find_nearest(self, array, value, location=False):
		import numpy as np
		array = np.asarray(array)
		idx = (np.abs(array - value)).argmin()
		if location is False:
			return array[idx]
		elif location is True:
			return idx
		else:
			print("location must be either True or False")
			
	def crossover_finder(self,y1,y2,preciseness,dif = 100):
		points=[]
		for i in range(len(y1)):
			if y1[i] - y2[i] <= dif:
				points.append(False)
			else:
				points.append(True)
		return points
				
	def ArrivalProbabilityatBSMInput(self,Plot=False,Rate=True,preciseness=10000,dif = 100,dispersion=True,minmax=1E-9):
		import numpy as np
		from itertools import compress 

		x1, y1, na1, na2 = self.CoincidenceOverlap(source_number=1,plot=True,
												   Jitter_offset=0, preciseness=10000,
												   dispersion=dispersion, minmax=minmax)
		x2, y2, na3, na4 = self.CoincidenceOverlap(source_number=2,plot=True,
												   Jitter_offset=self.sync_jitter, preciseness=10000,
												   dispersion=dispersion, minmax=minmax)

		
		sig1 = self.Link1.Transmission()*self.Source1.pair_count
		sig2 = self.Link2.Transmission()*self.Source2.pair_count
		
		if Rate is True:
			y1 = np.array(y1)*sig1
			y2 = np.array(y2)*sig2
		if Rate is False:
			y1 = np.array(y1)
			y2 = np.array(y2)
		
		bool_list = self.crossover_finder(y1,y2,preciseness,dif)
		
		y_new=[]
		for i in range(len(x1)):
			if bool_list[i] is True:
				y_new.append(y2[i])
			else:
				y_new.append(y1[i])
		
		y_new = np.array(y_new)

		cw = self.coincidence_window
		
		slice1 = self.find_nearest(x1,-cw/2.0,location=True)
		slice2 = self.find_nearest(x1,cw/2.0,location=True)
		
		x1    = x1[slice1:slice2]
		y_new = y_new[slice1:slice2]
		y1    = y1[slice1:slice2]
		y2    = y2[slice1:slice2]
		
		
		step = 2E-9/preciseness
		
		swapping_rate = sum(y_new*step)
		
		if Plot is True:
			return x1, y_new, y1, y2
		else:
			return swapping_rate 
		
		
		
		
			
		
class Link:
	"""The system that synchronises the pulses in the two sources
		Values:
		
		Length						The total length of fibre in this link
		Loss						The loss coefficient of this section of fibre.					
		
		Definitions:
		
		Transmission()
				The amount of signal that is seen at the output of the fibre
		
		Dispersion_Jitter()
				Currently no definition
		"""
	
	def __init__(self,Length,Loss=0.3,ITU_bandwidth_Hz = 100*10E9, peak_wavelength_nm=1550.12,
				 dispersion_coef=18, fibre_coupling_loss=0.08):
		"""Here Length is in km and Loss in in Db/km"""
		self.Length = Length
		self.Loss = Loss
		self.fibre_coupling_loss = fibre_coupling_loss
		self.dispersion_coef = dispersion_coef #ps/(nm*km)
		self.ITU_bandwidth = ITU_bandwidth_Hz
		self.peak_wavelength_nm = peak_wavelength_nm
		
	def Transmission(self):
		t = 10**(-(self.Length*self.Loss)/10)
		t_new = t-(2*t*self.fibre_coupling_loss)
		return t_new
		
	def Dispersion_Jitter(self):
		c = 299792458               #m/s
		d = self.dispersion_coef #ps/(nm*km)
		bandwidth_Hz = self.ITU_bandwidth
		peak_wavelength = self.peak_wavelength_nm*1E-9
		
		bandwidth_nm = (bandwidth_Hz*(peak_wavelength**2))/c
		length = self.Length
		
		return d*bandwidth_nm*length*1E-12 #s
		


	
	
class Detector:
	#This is taken from Alex Qui's project A
    """General Detector object.

    Public methods:
    paralyzable -- settable property, return True (False) if the
        detector is (not) paralyzable.
    base_efficiency -- settable property, return the base efficiency
        i.e. not accounting for dead time.
    jitter -- settable property, return detector jitter in picoseconds.
    dead_time -- settable property, return detector dead time in ps.
    dark rate -- settable property, return the rate of dark counts.
    efficiency(self, incidence_rate) --
        return the detector efficiency at a given incidence rate,
        accounting for dead time, dark rate and paralyzability.
    detection_rate(self, incidence_rate) --
        return the detection rate at a given incidence rate, accounting
        for dark counts and the detector efficiency at that incidence
        rate.
    """

    def __init__(self, paralyzable=True, base_efficiency=0.7,
                 jitter=100, dead_time=5000, dark_rate=1000):
        self.paralyzable = paralyzable
        self.base_efficiency = base_efficiency
        self.jitter = jitter
        self.dead_time = dead_time
        self.dark_rate = dark_rate

    @property
    def paralyzable(self):
        return self._paralyzable

    @paralyzable.setter
    def paralyzable(self, paralyzable):
        if paralyzable not in [True, False]:
            raise ValueError("paralyzable must be either True or False.")
        self._paralyzable = paralyzable
        return self._paralyzable

    @property
    def base_efficiency(self):
        return self._base_efficiency

    @base_efficiency.setter
    def base_efficiency(self, base_efficiency):
        if not 0 <= base_efficiency <= 1:
            raise ValueError("base_efficiency must be between 0 and 1.")
        self._base_efficiency = base_efficiency
        return self._base_efficiency

    @property
    def jitter(self):
        return self._jitter

    @jitter.setter
    def jitter(self, jitter):
        if jitter < 0:
            return ValueError("jitter cannot be less than 0.")
        self._jitter = jitter
        return self._jitter

    @property
    def dead_time(self):
        return self._dead_time

    @dead_time.setter
    def dead_time(self, dead_time):
        if dead_time < 0:
            return ValueError("dead_time cannot be less than 0.")
        self._dead_time = dead_time
        return self._dead_time

    @property
    def dark_rate(self):
        return self._dark_rate

    @dark_rate.setter
    def dark_rate(self, dark_rate):
        if dark_rate < 0:
            return ValueError("dark_rate cannot be less than 0.")
        self._dark_rate = dark_rate
        return self._dark_rate

    def _efficiency_p(self, incidence_rate):
        """Return the efficiency assuming a paralyzable detector.

        Calculate the effective efficiency from the incidence rate,
        taking into account the dead time and dark count rate.
        Assumes a paralyzable detector.
        """
        efficiency = (
            self.base_efficiency
            * math.exp(-self.dead_time*1E-12*incidence_rate))
        return efficiency

    def _efficiency_np(self, incidence_rate):
        """Return the efficiency assuming a non-paralyzable detector.

        Calculate the effective efficiency from the incidence rate,
        taking into account the dead time and dark count rate.
        Assumes a nonparalyzable detector.
        """
        efficiency = (
            self.base_efficiency
            * (1/(1 + incidence_rate*self.dead_time*1E-12)))
        return efficiency

    def efficiency(self, incidence_rate):
        """Return the detector efficiency given the incidence rate."""
        # Account for dark_rate here.
        incidence_rate += self.dark_rate
        if self.paralyzable is True:
            return self._efficiency_p(incidence_rate)
        elif self.paralyzable is False:
            return self._efficiency_np(incidence_rate)
        else:
            raise ValueError("self.paralyzable has been mis-set.")

    def detection_rate(self, incidence_rate):
        """Return the detection rate given the incidence rate.

        Accounts for dark rate and detector efficiency.
        """
        # Can't do an in-place addition because we need the raw
        # incidence_rate to plug into the efficiency method.
        total_rate = incidence_rate + self.dark_rate
        self._detection_rate = total_rate * self.efficiency(incidence_rate)
        return self._detection_rate





class User:
	"""The system that synchronises the pulses in the two sources
		Values:
		
		Length						The total length of fibre in this link
		Loss						The loss coefficient of this section of fibre.					
		
		Definitions:
		
		SignalCounts()
				The rate at which photons that are part of pairs
		
		SinglesCounts()
				The rate at which total photons are seen 
		"""
	
	def __init__(self,Link,Source):
		self.Link = Link
		self.Source = Source
		
	def SignalCounts(self):
		"""The signal seen be the User"""
		return self.Source.pair_count*self.Link.Transmission()
		
	def SinglesCounts(self):
		"""The Singles seen by the User"""
		return self.Source.singles_count*self.Link.Transmission()




		
class Repeater:
	"""The system that synchronises the pulses in the two sources
		Values:
		
		UserA
		UserB
		Source1
		Source2
		LinkA1
		LinkR1
		LinkB2
		LinkR2
		Length						The total length of fibre in this link
		Loss						The loss coefficient of this section of fibre.					
		
		Definitions:
		
		SignalA()
				The rate at which photon pairs are seen by Alice
				
		SignalB()
				The rate at which photon pairs are seen by Bob
		
		SinglesA()
				The rate at which total photons are seen by Alice
				
		SinglesB()
				The rate at which total photons are seen by Bob
				
		SignalAR()
				The Signal seen at the repeater from Source 1
		
		SignalBR()
				The Signal seen at the repeater from Source 2
		
		SinglesAR()
				The Singals seen at the repeater from Source 1
		
		SinglesBR()
				The Singals seen at the repeater from Source 2
				
		SwapSignal(minimum_coincidence)
				Returnes the rate at which entanglement is swap
				#Not curent fully correct
				
		Swapnoise()
				Calculates the noise rate in entanglement swapping
		
		find_nearest(array, value)
				Finds the array position of the nearet value to value
		
		JitterError(Tmin,TMax,TStep)
				This calculates the error associated with the swapping Jitter
		
		Jitter(Tmin,TMax,TStep,plot)
				This calculates the value of the swapping Jitter
		
		Accidentals()
				The Rate at which accidental photons appear to be swapped at the Reapeter
		"""
	
	def __init__(self,UserA,UserB,Source1,Source2,LinkA1,LinkR1,LinkB2,LinkR2,Sync):
		""""""
		self.UserA = UserA
		self.UserB = UserB
		self.Source1 = Source1
		self.Source2 = Source2
		self.LinkA1 = LinkA1
		self.LinkR1 = LinkR1
		self.LinkB2 = LinkB2
		self.LinkR2 = LinkR2
		self.Sync = Sync
		
	def SignalA(self):
		"""The Signal rate seen by Alice"""
		return self.UserA.SignalCounts()
	
	def SignalB(self):
		"""The Signal rate seen be Bob"""
		return self.UserB.SignalCounts()
	
	def SinglesA(self):
		"""The Singles rate seen by Alice"""
		return self.UserA.SinglesCounts()
	
	def SinglesB(self):
		"""The Singles rate seen by BoB"""
		return self.UserB.SinglesCounts()
		
	def SignalAR(self):
		"""The Signal rate seen at the Repeater that is entangled to Alice"""
		return self.Source1.pair_count*self.LinkR1.Transmission() 
	
	def SignalBR(self):
		"""The Signal rate seen at the Repeater that is entangled to Bob"""
		return self.Source2.pair_count*self.LinkR2.Transmission()
	
	def SinglesAR(self):
		"""The Singles rate seen at the Repeater from the shared source with Alice"""
		return self.Source1.singles_count*self.LinkR1.Transmission()
	
	def SinglesBR(self):
		"""The Singles rate seen at the Repeater from the shared source with Bob"""
		return self.Source2.singles_count*self.LinkR2.Transmission()

	def SwapSignal(self, minimum_coincidence = False, just_pc = False,dispersion=True):
		import numpy as np
		"""The Rate at which entanglement is swapped at the Reapeter"""
		
		c = 299792458               #m/s
		dz = abs(self.LinkR1.Length - self.LinkR2.Length)*0.0001              #m
		dOmega = abs(self.Source1.frequency2 - self.Source2.frequency2)    	#Hz
		if minimum_coincidence is True:
			coincidence_window = 2/(self.Source1.frequency_bandwidth + self.Source2.frequency_bandwidth)
		if minimum_coincidence is False:
			coincidence_window = self.Sync.coincidence_window

		Pc = 0.25*( 1 + np.exp(-0.5*(dz**2)*((1/(c*coincidence_window))**2))*np.exp(-0.5*((dOmega*coincidence_window)**2)))
		
		if just_pc is True:
			return Pc
		if just_pc is False:
			return Pc*self.Sync.ArrivalProbabilityatBSMInput(Plot=False,preciseness=100000,dif = 100,dispersion=dispersion)
	
	def SwapNoise(self,dispersion=True):
		import numpy as np
		"""The Possonian Noise on the BSM"""
		
		#Noise = np.sqrt((self.SinglesAR()*self.SinglesBR())*self.SwapSignal(minimum_coincidence = False, just_pc = True))
		Noise = np.sqrt(self.SwapSignal(minimum_coincidence = False, just_pc = False,dispersion=dispersion))
		
		return Noise
		
	def find_nearest(self, array, value, location=False):
		import numpy as np
		array = np.asarray(array)
		idx = (np.abs(array - value)).argmin()
		if location is False:
			return array[idx]
		elif location is True:
			return idx
		else:
			print("location must be either True or False")
		
		
	def JitterError(self,byX=False,Output_T=True,dispersion=True,minmax=1E-9):
		"""This calculates the Error associated with the synchronisation Jitter"""
		import numpy as np
		
		if self.Sync.FWHM_withdispersion(source_number = 1) >= self.Sync.FWHM_withdispersion(source_number = 2):
			x, y, na1, na2 = self.Sync.CoincidenceOverlap(source_number = 1, plot=True, Jitter_offset=0,
														  preciseness=10000,dispersion=dispersion,
														  minmax=minmax)			
		else:
			x, y, na1, na2 = self.Sync.CoincidenceOverlap(source_number = 2, plot=True, Jitter_offset=0,
														  preciseness=10000,dispersion=dispersion,
														  minmax=minmax)
		
		
		x=np.array(x)
		y=np.array(y)
		y=y/max(y)
		points=int(len(x))
		
		
		if byX is False:
			centre = np.where(y == max(y))[0]
		if byX is True:
			centre = np.where(x == min(abs(x)))[0]
		#print(len(centre))
		
		separations = []
		y_separations = []
		for i in np.arange(1,int((points-1)/2),1):
			separations.append(x[centre+i]-x[centre-i])
			y_separations.append(y[centre-i])
			
		loc = self.find_nearest(separations,self.Sync.sync_jitter,location=True)
		
		Error = (y[centre] - y_separations[loc][0])[0]
		
		RP = self.SwapSignal(minimum_coincidence = False, just_pc = False,dispersion=dispersion)
		Pn = self.SwapNoise(dispersion=dispersion)
		An = self.Accidentals(minimum_coincidence = False)
		
		T = (RP+Pn+An)/(1-Error)
		
		if Output_T is True:
			return Error*T
		else:
			return Error
		
				
	def Accidentals(self, minimum_coincidence = False):
		import numpy as np
		"""The Rate at which accidental photons appear to be swapped at the Reapeter"""
		
		c = 299792458               #m/s
		dz = abs(self.LinkR1.Length - self.LinkR2.Length)*0.0001              #m
		dOmega = abs(self.Source1.frequency2 - self.Source2.frequency2)    	#Hz
		
		if minimum_coincidence is True:
			coincidence_window = 2/(self.Source1.frequency_bandwidth + self.Source2.frequency_bandwidth)
		if minimum_coincidence is False:
			coincidence_window = self.Sync.coincidence_window
		
		Pc = 0.25*( 1 + np.exp(-0.5*(dz**2)*((1/(c*coincidence_window))**2))*np.exp(-0.5*((dOmega*coincidence_window)**2)))
		
		return Pc*(self.SinglesAR()*self.SinglesBR()*coincidence_window) #CHECK THIS WITH SID
		
	
	def FullBSMSignalAndNoise(self,dispersion=True,minmax=1E-9,SandN = True):
		
		RP = self.SwapSignal(minimum_coincidence = False, just_pc = False,dispersion=dispersion)
		Pn = self.SwapNoise(dispersion=dispersion)
		An = self.Accidentals(minimum_coincidence = False)
		Jn = self.JitterError(byX=False,Output_T=True,dispersion=dispersion,minmax=minmax)
		
		Signal = RP
		Noise = Pn+An+Jn
		
		if SandN is True:
			return Signal, Noise
		if SandN is "S":
			return Signal
		if SandN is "N":
			return Noise




			
class TwoPartyQKD:
	
	def __init__(self,Repeater):
		""""""
		self.UserA = Repeater.UserA
		self.UserB = Repeater.UserB
		self.Source1 = Repeater.Source1
		self.Source2 = Repeater.Source2
		self.LinkA1 = Repeater.LinkA1
		self.LinkR1 = Repeater.LinkR1
		self.LinkB2 = Repeater.LinkB2
		self.LinkR2 = Repeater.LinkR2
		self.Sync = Repeater.Sync
		self.Repeater = Repeater
		self.sync_jitter = self.Sync.sync_jitter
		
	def AandB_Signal(self):
		"""This is the rate at which the two users see an entangled pair."""
		
		#Unsure how to select Source to take the pair rate from.
		C = self.Source1.pair_count
		TA = self.LinkA1.Transmission()
		TB = self.LinkB2.Transmission()
		
		return  C*TA*TB
		
	def AandB_possonian(self):
		"""This returnes the False success rate."""
		import numpy as np
		
		return  np.sqrt(self.AandB_Signal())
		
	def AandB_Accidentals(self):
		""""""
		Sa = self.UserA.SinglesCounts()
		Sb = self.UserB.SinglesCounts()
		cw = self.Sync.coincidence_window
		Acc = Sa*Sb*cw
		
		return Acc
		
	def FWHM_withdispersion(self, Source, Link):
		import numpy as np
		
		pulse_width = Source.pulse_width
		pulse_width_new = np.sqrt(pulse_width**2 + Link.Dispersion_Jitter()**2)
		
		fwhm = pulse_width_new*(np.sqrt(2*np.log(2))/2)
		
		return fwhm
		
	def PulseShape(self, Plot=False, AandB=False, preciseness=10000, dispersion=True,minmax=1E-9):
		import numpy as np
		
		preciseness=preciseness+1
		Jitter_offset=self.sync_jitter#*1E-12
		

		if dispersion is False:
			A_fwhm = self.Sync.fwhm1
			S1 = self.Source1
			B_fwhm = self.Sync.fwhm2
			S2 = self.Source2
		elif dispersion is True:
			A_fwhm = self.FWHM_withdispersion(self.Source1,self.LinkA1)
			S1 = self.Source1
			B_fwhm = self.FWHM_withdispersion(self.Source2,self.LinkB2)
			S2 = self.Source2
		else:
			print("Error: value 'dispersion' must be True or False. 'dispersion' now set to True")
			A_fwhm = self.FWHM_withdispersion(self.Source1,self.LinkA1)
			S1 = self.Source1
			B_fwhm = self.FWHM_withdispersion(self.Source2,self.LinkB2)
			S2 = self.Source2

		#Alice
		Gs_Alice = []
		ts_Alice = []
		for t in np.linspace(-minmax,minmax,preciseness):
			A = -(2*(t)/A_fwhm)**2
			G = (2**A)
			Gs_Alice.append(G)
			ts_Alice.append(t)
			 
		Gs_Alice = np.array(Gs_Alice)
		ts_Alice = np.array(ts_Alice)
		step = 2*minmax/preciseness
		norm = np.sum(Gs_Alice*step)
		Gs_Alice = Gs_Alice/norm
		
		#Bob
		Gs_Bob = []
		ts_Bob = []
		Gs_Bob_norm = []
		for t in np.linspace(-minmax,minmax,preciseness):
			A = -(2*(t-Jitter_offset)/B_fwhm)**2
			A_norm = -(2*(t)/B_fwhm)**2
			G = (2**A)
			G_norm = (2**A_norm)
			Gs_Bob.append(G)
			Gs_Bob_norm.append(G_norm)
			ts_Bob.append(t)
		Gs_Bob = np.array(Gs_Bob)
		Gs_Bob_norm = np.array(Gs_Bob_norm)
		ts_Bob = np.array(ts_Bob)
		step = 2*minmax/preciseness
		norm = np.sum(Gs_Bob_norm*step)
		Gs_Bob = Gs_Bob/norm
		
		yA = Gs_Alice
		xA = ts_Alice
		yB = Gs_Bob
		xB = ts_Bob
		
		
		TraA = 1#self.LinkA1.Transmission()#*self.Source1.pair_count
		TraB = 1#self.LinkB2.Transmission()#*self.Source2.pair_count
		
		yA_atuser = np.array(yA)*TraA
		yB_atuser = np.array(yB)*TraB
		
		bool_list = self.Sync.crossover_finder(yA_atuser,yB_atuser,preciseness,10)
		
		y_new=[]
		for i in range(len(xA)):
			if bool_list[i] is True:
				y_new.append(yB_atuser[i])
			else:
				y_new.append(yA_atuser[i])
		
		y_new = np.array(y_new)

		cw = self.Sync.coincidence_window
		
		slice1 = self.Sync.find_nearest(xA,-cw/2.0,location=True)
		slice2 = self.Sync.find_nearest(xA,cw/2.0,location=True)
		
		x_new    = xA[slice1:slice2]
		y_new = y_new[slice1:slice2]
		yA_atuser    = yA_atuser[slice1:slice2]
		yB_atuser    = yB_atuser[slice1:slice2]
		
		
		step = 2*minmax/preciseness
		
		swapping_trans = sum(y_new*step)
		
		if Plot is True:
			if AandB is False:
				return x_new, y_new, yA_atuser, yB_atuser
			if AandB is True:
				return xA[slice1:slice2], yA[slice1:slice2], xB[slice1:slice2], yB[slice1:slice2]
		else:
			return swapping_trans
		
	def FourFold_withswap(self,minmax=1E-9):
		"""This returns The rate at which 4-fold coincidences are seen
		between Alice, Bob, Repeater_alpha, and Repeater_beta"""
		
		TA = self.LinkA1.Transmission()
		TB = self.LinkB2.Transmission()
		
		#Trans_overlap = self.PulseShape(Plot=False, preciseness=10000, dispersion=True,Jitter_offset_ps=Jitter_offset_ps):
		pulse_overlap = self.PulseShape(Plot=False, preciseness=10000, dispersion=True, minmax=minmax)
		
		QR_sig = self.Repeater.SwapSignal(minimum_coincidence = False, just_pc = False)
		
		return QR_sig*TA*TB*pulse_overlap
		
	def FourFold_trueAB_falseQR(self,minmax=1E-9, Possonian=True, preciseness=10000):
		"""This returns The rate at which 4-fold coincidences are seen
		between Alice, Bob, Repeater_alpha, and Repeater_beta where there
		was a flase swap at the repeater"""
		
		TA = self.LinkA1.Transmission()
		TB = self.LinkB2.Transmission()
		
		QR_sig, QR_noi = self.Repeater.FullBSMSignalAndNoise(minmax=minmax)
		if Possonian is False:
			QR_noi = QR_noi - self.Repeater.SwapNoise()
		
		pulse_overlap = self.PulseShape(Plot=False, preciseness=preciseness, dispersion=True, minmax=minmax)
		
		return QR_noi*TA*TB*pulse_overlap
		
	def FourFold_accidentals_trueswap(self):
		"""This returns the rate at which 4-fold accidentals are seen
		between Alice, Bob, Repeater_alpha, and Repeater_beta where
		a true swap happens at the Repeater and a false swap happens
		between Alica and Bob."""
		
		QR_sig = self.Repeater.SwapSignal(minimum_coincidence = False, just_pc = False)
		Sa = self.UserA.SinglesCounts()
		Sb = self.UserB.SinglesCounts()
		cw = self.Sync.coincidence_window
		
		Acc = QR_sig*Sa*Sb*cw*cw
		
		return Acc
			
	def FourFold_accidentals_falseswap(self):
		"""This returns the rate at which 4-fold accidentals are seen
		between Alice, Bob, Repeater_alpha, and Repeater_beta"""
		
		Sa = self.UserA.SinglesCounts()
		Sb = self.UserB.SinglesCounts()
		
		QR_noi = self.Repeater.Accidentals()
		
		cw = self.Sync.coincidence_window
		
		Acc = QR_noi*Sa*Sb*cw*cw
		
		return Acc
		
	def FourFold_possonian(self):
		"""This is the 4-fold possonian noise that is seen
		between Alice, Bob, Repeater_alpha, and Repeater_beta"""
		import numpy as np
		return np.sqrt(self.FourFold_withswap())
		
	def FourFold_full(self,minmax=1E-9, Possonian=True, preciseness=10000, SandN=True):
		
		Signal  = self.FourFold_withswap(minmax=minmax)
		NnoiatQR = self.FourFold_trueAB_falseQR(minmax=minmax, Possonian=Possonian, preciseness=preciseness)
		Nacc_true = self.FourFold_accidentals_falseswap()
		Nacc_false = self.FourFold_accidentals_trueswap()
		Npos = self.FourFold_possonian()
		#Nsye = self.FourFold_syncerror()
		Noise  = NnoiatQR+Nacc_true+Nacc_false+Npos#+Nsye
		
		if SandN is True:
			return Signal, Noise
		if SandN is "S":
			return Signal
		if SandN is "N":
			return Noise
	
	def QBER_estimate(self,minmax=1E-9, Possonian=True, preciseness=10000):
		Sig, Err_with_pos = self.FourFold_full(minmax=minmax,Possonian=Possonian, preciseness=preciseness)
		Npos = self.FourFold_possonian()
		NnoiatQR = self.FourFold_trueAB_falseQR(minmax=minmax, Possonian=Possonian, preciseness=preciseness)
		Nacc_true = self.FourFold_accidentals_falseswap()
		Nacc_false = self.FourFold_accidentals_trueswap()
		
		Err = NnoiatQR+Nacc_true+Nacc_false
		
		QBER = Err/(Sig-Npos)
		
		return QBER
		
		
		
	


if __name__ == "__main__":
	S1  = Source(10**6,10**7,193.4,193.4)
	S2  = Source(10**6,10**7,193.4,193.4)
	LA1 = Link(25)
	LR1 = Link(25)
	LB2 = Link(25)
	LR2 = Link(25)
	UA  = User(LA1,S1)
	UB  = User(LB2,S2)
	Sy  = Sync(S1,S2,5)
	R   = Repeater(UA,UB,S1,S2,LA1,LR1,LB2,LR2)
	TP  = TwoPartyQKD(R)
