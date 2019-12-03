# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:45:18 2019

@author: pcast
"""
import numpy as np
import MiniLFCollector as mlfc



#%%   Creacion del objeto Colector
    
#Diseño concentrador
#dimensiones en metros
W = 5.4
N_m = 11
w_m = 0.40
alt_col = 3.9
L = 12


#Diseño receptor
#Tomando como referencia el 0,0 en la esquina inf izq del receptor
p1 = np.array([0, 0])
p2 = np.array([40, 120])
#p5 = np.array([5, 60])
p3 = np.array([400, 120])
p4 = np.array([440, 0])

coord_recep = np.array([p1,p2,p3,p4])

#Descripción de Absorbedor
#tomando como referencia el 0,0 en el punto medio del lado inferior del receptor
dim_abs = 360           #dimension: ancho en caso de placa, radio del tubo. dimension en mm.
origen_abs = np.array([[0,100]])    #arreglo con origen de los absorbedores 
#Thermal geometry
w_port  = 2.5/1000        #Ancho de minicanales. Metro
h_port  = 2/1000        #Altura de minicanales. Metro
e_mc    = 0.3/1000        #espesor de minicanales. Metro


## Todo lo de arriba es fijo

#Condiciones de operación
DNI = 900              #Radiacion Solar[W/m2]
v_wind = 7              #Velocidad del viento [m/s]
T_amb = 273.15 + 30     #Temp. ambiente

#.......Condiciones de Operación
#Condiciones de Entrada
P_in = 0.2        #Presion entrada [MPa] 
m_in = 0.1    # equivalente a Re=2300   0.07 = 250
T_in = 273.15 + 90


#Angulo del sol
theta_sol = 10           #Angulo del sol en grados, donde 0° es el mediodia solar

#%%
#Colector
colector = mlfc.MiniLFCollector()
#Separar el angulo del sol

colector.construccion(W, w_m, N_m, alt_col, L, coord_recep, dim_abs, origen_abs, w_port, h_port, e_mc)

#colector.CondInicial(DNI, v_wind, T_amb, T_in, P_in, m_in*0.8)

eff, T_f, x_f, h_transf, h = colector.simulacion(theta_sol, DNI, v_wind, T_amb, T_in, P_in, m_in, plot = "y", corr="gungar")

#colector.plotColector()
