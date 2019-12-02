# -*- coding: utf-8 -*-
"""
Created on Fri May 10 11:03:11 2019

@author: pablo
"""


import numpy as np

#Funcion que genera una linea a partir de dos puntos
def line(p1, p2): 
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    
    return A, B, -C

#Funcion que genera un poligono a partir de:
    #type0: 0 para circunferencia, 1 para rectangulo
    #X1: 0, origen de circunferencia; 1, para una esquina
    #X2: 0, origen de circunferencia; 1, para la sig. esquina
    #dim: 0, diametro, 1, espesor
    #N: 0, no considerar; 1, direccion del poligono, UP-DOWN
    # 'up' = theta - 90, 'down' = theta + 90
    #por defecto la circunferencia tiene 14 lados
    
def polygon(type0, X1, X2, dim, up_down):
     
    if type0 == 0:
        #circunferencia
        theta = np.linspace(0, 2*np.pi, 21)
        X = dim*np.cos(theta) + X1[0]
        Y = dim*np.sin(theta) + X1[1]
        return X, Y
    
    elif type0 == 1:
        #rectangulo
        delta_y = X2[1] - X1[1]
        delta_x = X2[0] - X1[0]
        theta = np.arctan(delta_y/delta_x)%(2*np.pi)
        #theta: angulo de inclinacion entre los puntos X1, X2
        
        if up_down == 'up':
            theta = theta  - np.pi*0.5
        elif up_down == 'down':
            theta = theta  + np.pi*0.5
        
        #'up0  y 'down' definen si la cara generado mira hacia arriba o abajo
        
        X3 = [X2[0] + dim*np.cos(theta), X2[1] + dim*np.sin(theta)]
        X4 = [X1[0] + dim*np.cos(theta), X1[1] + dim*np.sin(theta)]
        
        X = [X1[0], X2[0], X3[0], X4[0], X1[0]] #Se agrega el punto 1 para cerrar la figura
        Y = [X1[1], X2[1], X3[1], X4[1], X1[1]]
        
        return X, Y
    
    else: 
        print('Valor fuera de rango')    
        
        return False
        

class Superficie:
    
    def __init__(self, absortividad, reflectividad, transmisividad):
        self.alfa = absortividad
        self.rho = reflectividad
        self.tau = transmisividad
        self.ray_abs = 0
        self.ray_ref = 0
        self.ray_trans = 0
        self.intercept = []
        
    def posicion(self, x1, y1, x2, y2):
        self.xplot = [x1, x2]
        self.yplot = [y1, y2]
        delta_y = y2-y1
        delta_x = x2-x1 
        self.origen = [x1 + delta_x*.5, y1 + delta_y*.5]
        self.largo = np.round(np.sqrt(delta_y**2 + delta_x**2),6)
        
        if delta_x == 0:
            self.angulo = 0.5     
        else:
            self.angulo = np.arctan(delta_y/delta_x)%(2*np.pi)
        
        self.line = line([x1,y1], [x2,y2]) 
        
    def interseccion(self, L2, theta, x0, y0):
        #L2 la linea del rayo desde el sol
        #theta es el angulo, x0, y0 origen 
        L1 = self.line
        D  = L1[0]*L2[1] - L1[1]*L2[0]
        Dx = L1[2]*L2[1] - L1[1]*L2[2]
        Dy = L1[0]*L2[2] - L1[2]*L2[0]
        
        x1 = self.origen[0]
        y1 = self.origen[1]
        
        if D != 0:
        #Si se intersecta. Siempre se va a intersectar, a menos que sean paralelas 
            xi = Dx/D
            yi = Dy/D
            
            if np.sqrt((x1 - xi)**2 + (y1 - yi)**2) <= self.largo/2:
            #Esta en el rango del objeto    
                
                if np.round(np.cos(theta),4) == 0: 
                #si cos(ang) es cero, significa que el rayo es vertical
                
                    lambda_y = np.round((yi - y0)/np.sin(theta), 5)
                    
                    if lambda_y > 0:
                        #si lambda es positivo, significa que el rayo avanza
                        #luego puedo retornar el valor de lambda
                            return [lambda_y, xi, yi]
                    
                    else: return [100,0,10]
                    
                else:
                    
                    lambda_x = np.round((xi - x0)/np.cos(theta), 3) 
                    lambda_y = np.round((yi - y0)/np.sin(theta), 3)
                    #calcular el lambda del vector
                
                    if lambda_x == lambda_y:
                    #si son iguales, estan bien calculados
                    
                        if lambda_x > 0:
                        #si lambda es positivo, significa que el rayo avanza
                        #luego puedo retornar el valor de lambda
                            return [lambda_x, xi, yi]
            
                        else:
                            return [100,10,10]
                    else:
                        return [100,20,10]
            else:
                return [100,30,10]
        else:
            return [100,40,10]
        
    def reflection(self, L2, theta, xi, yi):
    #En este objeto se intersecta, ahora hay que ver si se refleja, absorbe o transmite
        p = np.random.uniform(0,1)
        #distribución uniforme(0,1)
        if p < self.rho:
            theta1 = (2*self.angulo - theta)%(2*np.pi)
            tupla = (xi,yi)
            self.intercept.append(tupla)
            return xi, yi, theta1
        
        elif p < self.rho + self.alfa:
            self.ray_abs += 1
            tupla = (xi,yi)
            self.intercept.append(tupla)
            return xi, yi, 10
        else:
            self.ray_trans += 1
            theta1 = theta%(2*np.pi)
            tupla = (xi,yi)
            self.intercept.append(tupla)
            return xi, yi, theta1
        


class Colector:
    
    def __init__(self):
        pass
        
    def posicionEspejos(self, W, w_m, N_m, a):  
        #Ancho total W, ancho de cada espejo w_m, numero de espejo n_m, altura a 
        W_util = W - w_m
        delta_w = W_util/(N_m-1)
        x = np.zeros(N_m)
        theta = np.zeros(N_m)
        for i in range(N_m):
            x[i]= (w_m/2 - W/2) + i*delta_w 
            if x[i] == 0:
                theta[i] = (np.pi/2) - 0.00001
            else:
                theta[i]= np.arctan(a/x[i])
                if theta[i] < 0:
                    theta[i] = -theta[i]
                elif theta[i] < np.pi:
                    theta[i] = np.pi - theta[i] 
       
        self.centros = x
        self.theta = theta
        self.height = a
        self.w_m = w_m
        self.ancho = W
    
    def GeometriaReceptor(self, coord):
        #coord_receptor tiene que ser np.array
        coord = coord/1000
        
        a1 = (coord[0][0] - coord[-1][0]*0.5)
        b1 = coord[0][1] + self.height
        coord[0] = np.array([a1,b1])

        coord[1:] = coord[1:] + coord[0]
        
        self.receptor = coord

    def GeometriaAbsorbedor(self, id0, cantidad, dimension, origen, lados = 14):
        
        #Mueve la placa al sistema cartesiano maestro del colector
        if id0 == 0:
            origen = origen[0]
            x1 = origen[0]/1000 - dimension*0.5/1000
            y1 = origen[1]/1000 + self.height
            x2 = origen[0]/1000 + dimension*0.5/1000
            y2 = origen[1]/1000 + self.height
        
            self.absorbedor = [x1, y1, x2, y2]
            
        elif id0 == 1:
            XY_cil = []
            i = 0
            #Resolver despues para mas de una tuberia
            if len(origen) == cantidad:
                for i in range(cantidad):
                    theta = np.linspace(0,2*np.pi, (lados+1))
                    dim = dimension/1000
                    X1 = origen[i]/1000
                    X = dim*np.cos(theta) + X1[0]
                    Y = dim*np.sin(theta) + X1[1] + self.height
                    XY_cil.append([X, Y])
            
                self.absorbedor = XY_cil
            
            else: 
                print ("Valor no valido")
        else: 
            print ('Valor indicado no valido')
        

    def anguloRefleccion(self, theta_sol):
        N = np.shape(self.theta)
        j = theta_sol
        gamma = np.zeros(N)
        for i in range(N[0]):
            O = self.theta[i]
            beta = O - np.pi/2 + j
            n = O - beta/2
            gamma[i] = n + np.pi/2
        self.inclinacion = gamma


    



