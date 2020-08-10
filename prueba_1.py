# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:45:18 2019

@author: pcast
"""
import numpy as np
import MiniLFCollector as mLFC
import matplotlib.pyplot as plt

#%%   Definicion de las variables y parametros del colector.
    
#Diseño concentrador
#dimensiones en metros
W = 6.4                 #Ancho W
N_m = 9                 #Numero de Espejos, n
w_m = 0.50              #Ancho de cada espejo, w_m
alt_col = 4.5           #Altura del colector, a (hasta base del receptor)
L = 100                 #Largo del colector, L


#Diseño receptor

#El receptor esta formado por cuatro puntos, los cuales cierran un cuadrilatero.
#La referencia es el 0,0 en la esquina inferior izquierda del receptor. 


p1 = np.array([0, 0])
p2 = np.array([0.125*w_m*1000, 250])
#p5 = np.array([5, 60])
p3 = np.array([1.075*w_m*1000, 250])
p4 = np.array([1.2*w_m*1000, 0])

coord_recep = np.array([p1,p2,p3,p4])

#Descripción de Absorbedor
#El absorbedor corresponde al tubo de minicanales
#Referencia: el 0,0 está en el punto medio del lado inferior del receptor. 

dim_abs = 470                       #dimension: ancho en caso de placa, radio del tubo. dimension en mm.
origen_abs = np.array([[0,240]])    #arreglo con origen de los absorbedores 

#Thermal geometry
w_port  = 2.5/1000                  #Ancho de minicanales. Metros
h_port  = 2/1000                    #Altura de minicanales. Metros
e_mc    = 0.3/1000                  #Espesor de minicanales. Metros


## Todo lo de arriba es fijo. Diseño/Construcción del Colector. 


#.......Condiciones de ambientales

DNI = 750                       #Radiacion Solar[W/m2]
v_wind = 1                      #Velocidad del viento [m/s]
T_amb = 273.15 + 25             #Temp. ambiente (kelvin)
theta_sol = 10                  #Angulo del sol en grados, donde 0° es el mediodia solar.

#.......Condiciones de Operación

#Condiciones de Entrada
P_in = 4.0                          #Presion entrada [MPa] 
m_in = 0.7                          #Flujo masico en (kg/s)equivalente a Re=2300   0.07 = 250
T_in = 273.15 + 175                 #Temperatura de Entrada en Kelvin


#%%
#Ejemplo Colector

#Se crea un Objeto tipo MiniLFCollector. Esta vacio, no tiene parametros. 
clt = mLFC.MiniLFCollector()

#Se utiliza el metodo Contruccion para asignarle valores a los parametros de diseño del colector.
clt.construccion(W, w_m, N_m, alt_col, L, coord_recep, dim_abs, origen_abs, w_port, h_port, e_mc)

#Se utiliza el metodo Simulacion para simular optica y termicamente el colector con las variables entregadas.
clt.simulacion(theta_sol, DNI, v_wind, T_amb, T_in, P_in, m_in, plot = "n", corr="gungar")

#Ejemplo de la temperatura final del liquido y la perdidas termicas del colector.
t_f = clt.T_fl
q_loss = clt.Q_loss


