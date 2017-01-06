# -*- coding: utf-8 -*-
"""
Created on Sat June 12 15:26:38 2016
@author: AF
"""
import numpy as np
from scipy.fftpack import rfft, irfft, fftfreq
import matplotlib.pyplot as plt


class two_body:
    def __init__(self,T,dt,xe_0,vey_0):
        self.xe = [xe_0]
        self.ye = [0]
        self.vxe = [0]
        self.vye = [vey_0]
        self.t = [0]
        self.dt = dt
        self.N = int(T/(self.dt)) + 1
        
    def calculate(self):
        M_s = 2.0e30
        M_e = 6.0e24
        self.r = [np.sqrt(self.xe[0]**2+self.ye[0]**2)]
        for i in range(1,self.N):
            r_e = np.sqrt(self.xe[i - 1]**2+self.ye[i - 1]**2)
            ### velocity ----
            self.vxe.append(self.vxe[i - 1] - 4*np.pi**2*self.xe[i - 1]*self.dt/r_e**3)
            self.vye.append(self.vye[i - 1] - 4*np.pi**2*self.ye[i - 1]*self.dt/r_e**3)
            ## new position ---
            self.xe.append(self.xe[i - 1] + self.vxe[i]*self.dt)
            self.ye.append(self.ye[i - 1] + self.vye[i]*self.dt)
            ## time ----
            self.t.append(self.t[i - 1] + self.dt)
            self.r.append(np.sqrt(self.xe[i]**2+self.ye[i]**2))
        return 0 
        
    def plot(self):
        plt.figure(figsize = (8,8))
        plt.title('3-body simulation')
        plt.xlabel('x(AU)')
        plt.ylabel('y(AU)')
        plt.plot(self.xe,self.ye,label = 'Earth')
        plt.plot(self.xj,self.yj,label = 'Jupiter')
        r = 0.1
        theta = np.linspace(0,2*np.pi,36)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        plt.fill(x,y,'y',label = 'sun')
        plt.legend()
        plt.savefig('chapter4_4.18.png', dpi = 144)
        plt.show()
        
    def plot_rdis(self):
        #plt.figure(figsize = (8,8))
        plt.title('Distance distribution',fontsize = 20)
        plt.xlabel('Time(yr)',fontsize = 20)
        plt.ylabel('Distance(AU)',fontsize = 20)
        plt.plot(self.t,self.r, lw = 1)
        plt.legend()

class three_body:
    def __init__(self,T,dt,xe_0,vey_0,xj_0,vyj_0):
        self.xe = [xe_0]
        self.ye = [0]
        self.xj = [xj_0]
        self.yj = [0]
        self.vxe = [0]
        self.vye = [vey_0]
        self.vxj = [0]
        self.vyj = [vyj_0]
        self.t = [0]
        self.dt = dt
        self.N = int(T/(self.dt)) + 1
        
    def calculate(self):
        M_s = 2.0e30
        M_e = 6.0e24
        M_j = 1.9e27
        for i in range(1,self.N):
            r_e = np.sqrt(self.xe[i - 1]**2+self.ye[i - 1]**2)
            r_j = np.sqrt(self.xj[i - 1]**2+self.yj[i - 1]**2)
            r_ej = np.sqrt((self.xe[i - 1] - self.xj[i - 1])**2+(self.ye[i - 1] - self.yj[i - 1])**2)
            ### velocity ----
            self.vxe.append(self.vxe[i - 1] - 4*np.pi**2*self.xe[i - 1]*self.dt/r_e**3 - 4*np.pi**2*(M_j/M_s)*(self.xe[i - 1]-self.xj[i - 1])*self.dt/r_ej**3)
            self.vye.append(self.vye[i - 1] - 4*np.pi**2*self.ye[i - 1]*self.dt/r_e**3 - 4*np.pi**2*(M_j/M_s)*(self.ye[i - 1]-self.yj[i - 1])*self.dt/r_ej**3)
            self.vxj.append(self.vxj[i - 1] - 4*np.pi**2*self.xj[i - 1]*self.dt/r_j**3 - 4*np.pi**2*(M_e/M_s)*(self.xj[i - 1]-self.xe[i - 1])*self.dt/r_ej**3)
            self.vyj.append(self.vyj[i - 1] - 4*np.pi**2*self.yj[i - 1]*self.dt/r_j**3 - 4*np.pi**2*(M_e/M_s)*(self.yj[i - 1]-self.ye[i - 1])*self.dt/r_ej**3)
            ## new position ---
            self.xe.append(self.xe[i - 1] + self.vxe[i]*self.dt)
            self.ye.append(self.ye[i - 1] + self.vye[i]*self.dt)
            self.xj.append(self.xj[i - 1] + self.vxj[i]*self.dt)
            self.yj.append(self.yj[i - 1] + self.vyj[i]*self.dt)
            ## time ----
            self.t.append(self.t[i - 1] + self.dt)
        return 0 
        
    def plot(self):
        #plt.figure(figsize = (8,8))
        plt.title('Effect of Jupiter on asteroid', fontsize = 20)
        plt.xlabel('x(AU)', fontsize = 20)
        plt.ylabel('y(AU)', fontsize = 20)
        plt.plot(self.xe,self.ye,label = 'Asteroid')
        #plt.plot(self.xj,self.yj,label = 'Jupiter')
        r = 0.1
        theta = np.linspace(0,2*np.pi,36)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        plt.fill(x,y,'y',label = 'sun',fontsize = 20)
        plt.legend()
        #plt.savefig('chapter4_4.18.png', dpi = 144)
        #plt.show()
        
class asteriod(three_body):
     def calculate(self):
        M_s = 2.0e30
        M_e = 1.0
        M_j = 1.9e27
        r_ej = np.sqrt((self.xe[0] - self.xj[0])**2+(self.ye[0] - self.yj[0])**2)
        r_e = np.sqrt(self.xe[0]**2+self.ye[0]**2)
        self.r = [np.sqrt(self.xe[0]**2+self.ye[0]**2)]
        self.e_o = [0.5*M_e*(self.vxe[0]**2+self.vye[0]**2)-4*np.pi**2*M_e/r_e - 4*np.pi**2*(M_j/M_s)*M_e/r_ej]
        self.e_p = [-4*np.pi**2*M_e/r_e - 4*np.pi**2*(M_j/M_s)*M_e/r_ej]
        self.e_k = [0.5*M_e*(self.vxe[0]**2+self.vye[0]**2)]
        for i in range(1,self.N):
            r_e = np.sqrt(self.xe[i - 1]**2+self.ye[i - 1]**2)
            r_j = np.sqrt(self.xj[i - 1]**2+self.yj[i - 1]**2)
            r_ej = np.sqrt((self.xe[i - 1] - self.xj[i - 1])**2+(self.ye[i - 1] - self.yj[i - 1])**2)
            ### velocity ----
            self.vxe.append(self.vxe[i - 1] - 4*np.pi**2*self.xe[i - 1]*self.dt/r_e**3 - 4*np.pi**2*(M_j/M_s)*(self.xe[i - 1]-self.xj[i - 1])*self.dt/r_ej**3)
            self.vye.append(self.vye[i - 1] - 4*np.pi**2*self.ye[i - 1]*self.dt/r_e**3 - 4*np.pi**2*(M_j/M_s)*(self.ye[i - 1]-self.yj[i - 1])*self.dt/r_ej**3)
            self.vxj.append(self.vxj[i - 1] - 4*np.pi**2*self.xj[i - 1]*self.dt/r_j**3 - 4*np.pi**2*(M_e/M_s)*(self.xj[i - 1]-self.xe[i - 1])*self.dt/r_ej**3)
            self.vyj.append(self.vyj[i - 1] - 4*np.pi**2*self.yj[i - 1]*self.dt/r_j**3 - 4*np.pi**2*(M_e/M_s)*(self.yj[i - 1]-self.ye[i - 1])*self.dt/r_ej**3)
            ## new position ---
            self.xe.append(self.xe[i - 1] + self.vxe[i]*self.dt)
            self.ye.append(self.ye[i - 1] + self.vye[i]*self.dt)
            self.xj.append(self.xj[i - 1] + self.vxj[i]*self.dt)
            self.yj.append(self.yj[i - 1] + self.vyj[i]*self.dt)
            ## time ----
            self.t.append(self.t[i - 1] + self.dt)
            ## radius ---
            self.r.append(np.sqrt(self.xe[i]**2+self.ye[i]**2)) 
            ## energy ---
            self.e_o.append(0.5*M_e*(self.vxe[i]**2+self.vye[i]**2)-4*np.pi**2*M_e/(np.sqrt(self.xe[i]**2+self.ye[i]**2)))
            self.e_p.append(-4*np.pi**2*M_e/r_e - 4*np.pi**2*(M_j/M_s)*M_e/r_ej)
            self.e_k.append(0.5*M_e*(self.vxe[i]**2+self.vye[i]**2))
        print np.mean(self.e_p)
        return 0 
        
     def plot(self):
         #plt.figure(figsize = (8,8))
         plt.title('Effect of Jupiter on asteroid', fontsize = 20)
         plt.xlabel('x(AU)', fontsize = 20)
         plt.ylabel('y(AU)', fontsize = 20)
         plt.plot(self.xe,self.ye,'o',markersize=0.1)
         #plt.plot(self.xj,self.yj,label = 'Jupiter')
         r = 0.1
         theta = np.linspace(0,2*np.pi,36)
         x = r*np.cos(theta)
         y = r*np.sin(theta)
         plt.fill(x,y,'y',label = 'sun')
        # plt.legend()
        # plt.savefig('chapter4_4.18a.png', dpi = 144)
        # plt.show() 
         
     def plot_rdis(self):
         #plt.figure(figsize = (8,8))
         plt.title('Distance distribution',fontsize = 20)
         plt.xlabel('Time(yr)',fontsize = 20)
         plt.ylabel('Distance(AU)',fontsize = 20)
         plt.plot(self.t,self.r, lw = 1)
         plt.legend()
         #plt.show()
        
     def radius_fft(self):
         radius_fft = np.fft.fft(self.r)
         plt.figure(figsize = (10,8))
         plt.plot(radius_fft.real)
         plt.xlim(0,200)
         plt.ylim(-5000,18000)
         '''t = np.arange(256)
         sp = np.fft.fft(np.sin(t))
         plt.plot(sp.real)
         plt.plot(sp.imag)'''
         plt.savefig('chapter4_4.18_fft.png', dpi = 256)
         plt.show()
         return 0
    
     def energy_plot(self):
         #plt.figure(figsize = (8,8))
         plt.title('Orbital Energy',fontsize = 20)
         plt.xlabel('Time(yr)',fontsize = 20)
         plt.ylabel('Energy', fontsize = 20)
         plt.plot(self.t,self.e_o, lw = 1)
         #plt.plot(self.t,self.e_p, lw = 1,label = 'Potential energy')
         #plt.plot(self.t,self.e_k, lw = 1,label = 'Kinetic energy')
         #plt.legend()
         #plt.savefig('chapter4_4.18_energy_1.png', dpi = 256)
         #plt.show()
         
     def energy_dispose(self):
         energy_record = []
         time_record = []
         #print self.e_o
         for i in range (1,len(self.e_o) - 1):
             if (self.e_o[i] > self.e_o[i - 1] and self.e_o[i] > self.e_o[i + 1]):
                 energy_record.append(self.e_o[i])
                 time_record.append(self.t[i])
         #print np.mean(energy_record)
         #print np.std(energy_record) 
         #plt.figure(figsize = (8,8))
         plt.title('Energy',fontsize = 20)
         plt.xlabel('Time(yr)',fontsize = 20)
         plt.ylabel('Energy', fontsize = 20)
         plt.plot(time_record,energy_record)
         #plt.savefig('chapter4_4.18_energy_d276.png', dpi = 256)
         #plt.show()
         '''
         x=np.array(self.e_o)
         n = len(self.e_o)
         freq = fftfreq(len(self.e_o), d=self.dt)
         freq = np.array(abs(freq))
         #y=np.fft.fft(x)
         y = rfft(x)
         y=np.array(abs(y))
         '''
         '''
         plt.xlim(1,n)
         plt.ylim(0,0.3)
         plt.bar(range(n),y,width = 0.2)
         plt.show()'''
         return 0
#### -----------------------    
'''        
plt.figure(figsize = (8,8))         
distance = 3.2    
for i in range(7): 
    v = np.sqrt(4*np.pi**2/distance)        
    A = asteriod(400,0.001,distance,v,5.45,0.42*2*np.pi)
    A.calculate()
    #A.plot_rdis()
    #A.energy_plot()
    data = A.energy_dispose() 
    n = len(data)
    plt.plot(data,lw = 1, label = str(distance))
    distance+=0.02
    #A.radius_fft()
for i in range(7):   
    v = np.sqrt(4*np.pi**2/distance)       
    A = asteriod(400,0.001,distance,v,5.45,0.42*2*np.pi)
    A.calculate()
    #A.plot_rdis()
    #A.energy_plot()
    data = A.energy_dispose() 
    n = len(data)
    plt.plot(data,':',linewidth = 1, label = str(distance))
    distance+=0.02
plt.xlim(1,n)
plt.xlabel('Frequency')
plt.ylabel('Intensity')
plt.ylim(0,0.5) 
plt.legend()
plt.savefig('chapter4_4.18_e_fft_all.png',dpi = 256)
plt.show()
'''
'''
radius = 3.276
v = np.sqrt(4*np.pi**2/radius)
A = asteriod(400,0.001,radius,v,5.45,0.42*2*np.pi)
#radius = 3.42
#v = np.sqrt(4*np.pi**2/radius)
#B = asteriod(400,0.001,radius,v,5.45,0.42*2*np.pi)
A.calculate()
#B.calculate()
#A.plot_rdis()
plt.figure(figsize = (16,8))
plt.subplot(121)
A.plot_rdis()
#B.plot_rdis()
plt.legend(['Initial radius = 3.276'])
plt.subplot(122)
A.energy_plot()
#B.energy_plot()
plt.legend(['Initial radius = 3.276'])
plt.savefig('chapter4_4.18_re_3.276.png',dpi = 256)
plt.show()
'''
## -------- orbital energy curve ------------
'''
M_s = 2.0e30
M_e = 1.0
M_j = 1.9e27
e_o = []
xe = 3.2
xe_record = []
ye = 0
xj = 5.45
yj = 0
vye = 3.471
vxe = 0
vxj = 0
vye = 0.42*2*np.pi 
for i in range(20):
    r_ej = np.sqrt((xe - xj)**2+(ye - yj)**2)
    r_e = np.sqrt(xe**2+ye**2)
    xe_record.append(xe)
    e_o.append(0.5*M_e*(vxe**2+vye**2)-4*np.pi**2*M_e/r_e - 4*np.pi**2*(M_j/M_s)*M_e/r_ej)
    xe+=0.01
plt.figure(figsize = (8,8))
plt.ylabel('Orbital energy')
plt.xlabel('Initial orbital radius')
plt.plot(xe_record,e_o)
plt.savefig('chapter4_4.18_initial_orbital_e.png',dpi = 256)
plt.show()'''
### ------------------------------------------
'''
plt.figure(figsize = (24,8))
plt.subplot(131)
radius = 3.276
v = np.sqrt(4*np.pi**2/radius)
A = two_body(400,0.001,radius,v)
A.calculate()
plt.ylim(3.0,3.6)
A.plot_rdis()
plt.subplot(132)
radius = 3.1
v = np.sqrt(4*np.pi**2/radius)
A = asteriod(400,0.001,radius,v,5.45,0.42*2*np.pi)
A.calculate()
plt.ylim(3.0,3.6)
A.plot_rdis()
plt.subplot(133)
radius = 3.276
v = np.sqrt(4*np.pi**2/radius)
A = asteriod(400,0.001,radius,v,5.45,0.42*2*np.pi)
A.calculate()
plt.ylim(3.0,3.6)
A.plot_rdis()
plt.savefig('final_rd3c.png',dpi = 256)
plt.show()
'''
#### ----------o_c -------------------
plt.figure(figsize = (16,8))
plt.subplot(121)
radius = 2.85
v = np.sqrt(4*np.pi**2/radius)
A = asteriod(200,0.001,radius,v,5.45,0.42*2*np.pi)
A.calculate()
#A.energy_plot()
A.plot()
'''
radius = 3.29
v = np.sqrt(4*np.pi**2/radius)
A = asteriod(200,0.001,radius,v,5.45,0.42*2*np.pi)
A.calculate()
A.plot()
radius = 3.8
v = np.sqrt(4*np.pi**2/radius)
A = asteriod(200,0.001,radius,v,5.45,0.42*2*np.pi)
A.calculate()
A.plot()
plt.legend(['Initial radius = 3.000 AU','Initial radius = 3.29 AU','Initial radius = 3.8 AU'],loc = 10, fontsize = 18)
'''

plt.subplot(122)
radius = 2.85
v = np.sqrt(4*np.pi**2/radius)
A = asteriod(400,0.001,radius,v,5.45,0.42*2*np.pi)
A.calculate()
#A.plot()
plt.ylim(2.8,2.9)
A.plot_rdis() 
#plt.plot(data[0],data[1],linewidth = 1)
#plt.legend(['Initial radius = 3.29 AU'],fontsize = 18)

plt.savefig('final_oe_c.png',dpi = 256)
plt.show()