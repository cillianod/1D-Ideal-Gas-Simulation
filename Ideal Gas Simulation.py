print("Cillian O'Donnell")
print('20333623')
print()
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.stats import maxwell

#Useful constants, initial conditions
R = 1 # for simplicity here, R=1
T0 = 1 #Temperature
F = 10 #Constant Force
m = 1 #particle mass
M = 100 #piston mass
N = 1000 #particle number
#t = 0
Leq = (R*N*T0)/F #equilibrium length
kb = 1
forces = [0.1,0.3,1,3,10,53,71,100] #for varying forces

X0 = 2*Leq
dist = np.linspace(-X0,X0,1000)
xi0 = np.random.choice(dist,N)
equipart = np.sqrt(1*kb*T0)
#Using this random gaussian dist. of velocities

vi0 = np.random.normal(loc=0,scale = equipart,size=N)
V0 = np.random.normal(loc=0,scale = equipart,size=1)
t=0


v_c = np.random.normal(loc=0.1,scale=equipart,size=1)
v_cs = [float(-v_c),float(v_c)]
#vi0 = np.random.choice(v_cs,size=N)
#Uncomment this to use an intiial bipolar distribution

#Equations of motion - lowercase for particles, uppercase for Piston
def move_xi_t(xi0,vi0,t):
    return xi0 + vi0*t

def move_Xt (X,V,t,F):
    return X + V*t - (F/(2*M))*t**2

def move_Vt(V,t,F):
    return V - (F/M)*t

def wait(x,v,X,V,F):
    if v > 0:
        tau = (M/F)*((V-v)+np.sqrt(((V-v)**2)-2*(F/M)*(x-X)))
        return tau
    if v < 0:
        tau = M/F*((V+v)+np.sqrt(((V+v)**2)+2*(F/M)*(x+X)))
        return tau

#collision of particle and piston function
def collision(v,V):
    if v > 0:
        w = (2*M*V+(m-M)*v)/(m+M)
        W = (2*m*v+(M-m)*V)/(m+M)
        return w,W
    elif v < 0:
        w = (-(2*M*V)+(m-M)*v)/(m+M)
        W = (-(2*m*v)+(M-m)*V)/(m+M)
        return w, W


avgX_list = []
X_ideal_list = []


#In this form, this will range through many forces, looping over numerous
#oscillations of the piston. The motion of the piston is plotted. Also plotted
#is the average piston position after a number of collisions for varying 
#forces.
for F in forces:
    microstate=[0]*N
    for i in range(N):
        microstate[i] = [xi0[i],vi0[i]]
    
    Leq = (R*N*T0)/F

    Vt = V0
    Xt = X0
    t=0
    v_list = []
    series_X = []
    series_t = []
    series_V = []
    series_X.append(X0)
    series_V.append(V0)
    series_t.append(t)
    v_list_list = []
    dist_times = []
    #loop over time intervals for part(c)
    for l in range(1):
        #create list of all velocities after time interval for part (c)
        for i in range(N):
            particle = microstate[i]
            velocity = particle[1]
            v_list.append(float(velocity))
        #plot velocity distribution after this time interval     
    
        v_list_list.append(v_list)
        v_list = []
        dist_times.append(t)
    
    #loop over a large number of collisions
        for c in range(1000):
            time0 = 1*10**10
            #find the time to next collision
            for j in range(N):
                particle = microstate[j]
                x_i = particle[0]
                v_i = particle[1]
                tau1 = wait(x_i,v_i,Xt,Vt,F)
                if tau1 <= time0:
                    time0 = tau1
                    index = j
                
            shortwait = time0
            bouncer = microstate[index]
        #update positions and velocities by this wait time
            for i in range(N):
                particle = microstate[i]
                x_i = particle[0]
                v_i = particle[1]
                x_it = move_xi_t(x_i,v_i,shortwait)
                microstate[i] = [x_it,v_i]
            Xt = move_Xt(Xt,Vt,shortwait,F)
            Vt = move_Vt(Vt,shortwait,F)
        
            #find bounces for piston and particle
            v_before = bouncer[1]
            V_before = Vt
            collide =  collision(v_before,V_before)
            v_after = collide[0]
            Vt = collide[1]
            t = t+shortwait
            series_V.append(Vt)
            series_X.append(Xt)
            series_t.append(t)
            X_final = Xt
            microstate[index] = bouncer[0],v_after 
        
#plots are broken into cells so the above can be run once and plots can be
#changed / run multiple times without long processing time
#%%
#plotting for piston position over time
avgX = sum(series_X)/len(series_X)
n = 500
line = [avgX]*len(series_t)
line2 = [Leq]*len(series_t)
plt.title('Piston position over time')
plt.plot(series_t,series_X,label='Piston position')
plt.plot(series_t[::n],line[::n],'.',label='Average position',ms = 5)
plt.plot(series_t[::n],line2[::n],'.',label='Expected equilibrium position',ms = 5)
plt.xlabel('Time')
plt.ylabel('Piston position')
plt.legend()
plt.show()
#%%
#plotting for average piston position for a range of forces
plt.loglog(avgX_list,forces,'.',label='Calculated position')
plt.xlabel('Average X position')
plt.ylabel('Force')
plt.title('Average Piston position across a range of Forces')
plt.loglog(X_ideal_list,forces,'--',label='Ideal Gas Law Estimation')
plt.legend()
plt.show()

#%%
#plotting for the equilibration of an intial bipolar distribution.
x = np.linspace(4,-4,1000)
fig = plt.figure()
ax =  plt.axes(projection='3d')
nbins = 30
data = v_list_list[6]
params = maxwell.fit(data,floc=-4)
stdev = maxwell.std(params[0],params[1])
mean = maxwell.mean(params[0],params[1])

print('The mean of the maxwell distribution is = %.3g' %mean)
print('The standard deviation of the maxwell distribution is = %.3g' %stdev)

for b in range(7):
    vs = v_list_list[b]
    ts = dist_times[b]
    hist,bins = np.histogram(vs,bins=nbins)
    xline = (bins[:-1] + bins[1:])/2
    ax.bar(xline,hist,ts,zdir='y',alpha = 0.8,label = 'T = %.1g'%ts)
tconst = [0]*len(x)

ax.set_title('Velocity distribution over time for initial bipolar distribution',loc='right')
ax.set_xlabel('Velocity')
ax.set_ylabel('Time')
ax.set_zlabel('Velocity Distribution')
plt.legend(loc='upper left',prop={'size': 7})
plt.show()
fig = plt.figure()
ax =  plt.axes(projection='3d')
ax.plot(x,tconst,maxwell.pdf(x, params[0],params[1]), label="Maxwell Distribution",color='purple')
ax.set_title('Theoretical Maxwell Distribution of Velocities')
ax.set_xlabel('Velocity')
ax.set_ylabel('Time')
ax.set_zlabel('Velocity Distribution')
plt.show()