# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 15:39:18 2014

@author: HP
"""
from math import sin, cos, log, ceil,pi
from scipy import optimize
import numpy
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def f(u):
    v = u[0]
    theta = u[1]
    x = u[2]
    y = u[3]
    return numpy.array([-g*sin(theta) - C_D/C_L*g/v_t**2*v**2,
                      -g*cos(theta)/v + g/v_t**2*v,
                      v*cos(theta),
                      v*sin(theta)])

def euler_step(u, f, dt):
    return u + dt * f(u)    
def rk2_step(u,f,dt):
    u_star=u+0.5*dt*f(u)
    return u+dt*f(u_star)
  
#Параметры бумажного самолетика:
g = 9.8     
v_t = 4.9     
C_D = 1./5
C_L = 1.0   

#Начальные условия
v0 = 10
theta0 = -0.5
x0 = 0.0  
y0 = 2.
T = 10.0
dt = 0.001                           
N = int(T/dt) +1  #шаг с которым мы разбиваем +1                
t = numpy.linspace(0.0, T, N) 

#Нарисуем путь самолетика
#plt.figure(figsize=(8,6))
#plt.grid(True)
#plt.xlabel(r'x', fontsize=18)
#plt.ylabel(r'y', fontsize=18)
#plt.title('Glider trajectory, flight time = %.2f , Distance = %.2f' % (T, x[-1]), fontsize=18) # форматирование строк 
#plt.plot(x,y, 'k-', lw=2)

def find_dist(v0, t0):
    u = numpy.zeros((1, 4)) #столбец по 4 элемента СОЗДАНИЕ

    u[0] = numpy.array([v0, theta0, x0, y0]) #u[0] -ЕГО ЗАПОЛНЕНИЕ
    
    #Метод Эйлера
    y=y0
    n=0#колчисе0ство столбцов
    while y >0:
        u = numpy.append(u,[rk2_step(u[n], f, dt)],axis =0) # 
        y = u[n+1,3] 
        n=n+1
        T= n*dt
    
    #Найдем расположение самолетика с течением времени
    x = u[:,2]
    return x[-1] 

def find_dist_for_circle(v0, t0):
   # u = numpy.zeros((1, 4)) #столбец по 4 элемента СОЗДАНИЕ

    u = numpy.array([v0, theta0, x0, y0]) #u[0] -ЕГО ЗАПОЛНЕНИЕ
    
    #Метод Эйлера
    y=y0
    n=0#колчисе0ство столбцов
    while y >0:
        u=rk2_step(u, f, dt)
        #u = numpy.append(u,[rk2_step(u[n], f, dt)],axis =0) # 
        y = u[3] 
        n=n+1
       
    
    #Найдем расположение самолетика с течением времени
    
    return u[2] 

print find_dist(v0, theta0)            
#dists = []
#v0s =numpy.linspace(0.1, 15, 1000)
#t0s = numpy.linspace(-1 , pi/3, 10)
#for v0 in v0s:
#    for t0 in t0s:
#        
#        dists.append(find_dist(v0, t0))
#        
#print max(dists)

def inverse_find_dist(x):
    v0=x[0]
    theta0=x[1]
    
    return 1./find_dist_for_circle(v0, theta0)

y = [v0,theta0]

res = optimize.brute(inverse_find_dist,((0.1,15),(-1,pi/3)),Ns=30,finish=optimize.fmin)
print res
   
#dists = []
#v0s =numpy.linspace(0.1, 15, 1000)
#t0s = numpy.linspace(-1 , pi/3, 10)
#for v0 in v0s:
#    for t0 in t0s:
#        
#        dists.append(find_dist(v0, t0))
#        
#print max(dists)
#
#v0=x[0]
#theta0=x[1]
