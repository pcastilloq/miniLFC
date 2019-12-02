# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:53:37 2019

@author: pablo
"""
import numpy as np
import matplotlib.pyplot as plt
import Raytrace as rt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
        
##Colector de la memoria##
def simulacion(theta_sol, rayos, lados_poly):        
    #Diseño concentrador
    #dimensiones en metros
    W = 7.5
    N_m = 11
    w_m = 0.50
    alt_col = 3.9
    
    
    #El cielo xd
    #Altura del "cielo" es 1.1 veces la altura del colector
    alt_cielo = alt_col*1.1
    
    
    #Diseño receptor
    #Tomando como referencia el 0,0 en la esquina inf izq del receptor
    p1 = np.array([0, 0])
    p2 = np.array([50, 150])
    #p5 = np.array([5, 60])
    p3 = np.array([550, 150])
    p4 = np.array([600, 0])
    
    coord_recep = np.array([p1,p2,p3,p4])
    
    #Descripción de Absorbedor
    #tomando como referencia el 0,0 en el punto medio del lado inferior del receptor
    id0 = 0                 #0 si minichannel plate, 1 si copper tube
    amount = 1              #1 si es una placa. Especificar numeros de tubos
    dim_abs = 500           #dimension: ancho en caso de placa, radio del tubo. dimension en mm.
    origen_abs = np.array([[0,100]])    #arreglo con origen de los absorbedores 
    lados = lados_poly

    #Angulo del Sol
    
    #Maxima altura emision de rayos
    max_alt_cielo = (alt_cielo + W*(np.cos(np.radians(theta_sol))*np.sin(np.radians(theta_sol))))*1.1
    
    ##Crear objeto Colector
    colector = rt.Colector()
    colector.posicionEspejos(W, w_m, N_m, alt_col)
    colector.GeometriaReceptor(coord_recep)
    colector.anguloRefleccion(np.radians(theta_sol))
    colector.GeometriaAbsorbedor(id0,amount,dim_abs, origen_abs, lados)
    
    #Lista con la lista de las 4+1 esquinas de cada espejo
    XY_mir = []
    XY_rec = []
    for i in range(N_m):
        x =  colector.centros[i]
        y = 0 
        x1 = x + colector.w_m*np.cos(colector.inclinacion[i])*0.5
        x2 = x - colector.w_m*np.cos(colector.inclinacion[i])*0.5
        y1 = y + colector.w_m*np.sin(colector.inclinacion[i])*0.5
        y2 = y - colector.w_m*np.sin(colector.inclinacion[i])*0.5
        
        #up y down definen hacia donde se "extruye" el cuerpo
        x_mir, y_mir = rt.polygon(1, [x1,y1], [x2,y2], 0.005, 'up')
        
        XY_mir.append([x_mir, y_mir])
    
    for i in range(len(colector.receptor)-1):
        X1 = colector.receptor[i]
        X2 = colector.receptor[i+1]
        x_rec, y_rec = rt.polygon(1, X1, X2, 0.01, 'down')
        
        XY_rec.append([x_rec, y_rec])
    
        
    
    #Diccionario creador de objetos del colector
    surface = {}
    #Cada espejo tiene 4 lados, cada pared del receptor tiene 4 lados, el vidrio se considera 1-dim, 
    N_surf = (N_m + (len(colector.receptor)-1))*4 + 1 + amount*(lados**id0) + 2
    
    for x in range(N_surf):
        if x < N_m*4:
        #espejos
    
            #Cuenta el espejo en la que va el index
            x_mir = x%4
            i = int(np.floor(x/4))
            if x_mir == 0:
                alpha = 0.0
                rho = 1.0
                tau = 0   
                          
            else:
                alpha = 1
                rho = 0
                tau = 0
                
            x1 = XY_mir[i][0][x_mir]
            y1 = XY_mir[i][1][x_mir]
            x2 = XY_mir[i][0][x_mir+1]
            y2 = XY_mir[i][1][x_mir+1]
    
            
        elif x < (N_m + len(colector.receptor) - 1)*4:
            #paredes
            alpha = 1
            rho = 0
            tau = 0
            #Cuenta la pared del receptor en la que va el index
            x_rec = (x - N_m*4)%4
            i = int(np.floor((x- N_m*4)/4))
            
            #Receptor y espejos secundarios
            if x_rec == 0:
                alpha = 0.0
                rho = 1.0
                tau = 0
            
            x1 = XY_rec[i][0][x_rec]
            y1 = XY_rec[i][1][x_rec]
            x2 = XY_rec[i][0][x_rec+1]
            y2 = XY_rec[i][1][x_rec+1]
            
        elif x == (N_m + len(colector.receptor) - 1)*4:
            #Modificar esto si se agrega opcion de remover cubierta de vidrio
            #ultimo cuerpo del receptor: cubierta de vidrio
            alpha = 0.0
            rho = 0.0
            tau = 1.0
            
            x1 = colector.receptor[-1][0]
            y1 = colector.receptor[-1][1]
            x2 = colector.receptor[0][0]
            y2 = colector.receptor[0][1]
            
        elif x < (N_m + (len(colector.receptor)-1))*4 + 1 + amount*(lados**id0):
        #ante-penultimo cuerpo a generar: absorbedor
            alpha = 1
            rho = 0.0
            tau = 0
            
            i = x - (N_m + (len(colector.receptor)-1))*4 + 1
            if id0 == 0:
                #minichannel plate
                x1, y1, x2, y2 = colector.absorbedor
            else: 
                #tuberias
                i = i%lados
                j = int(np.floor(i/lados))
                
                #Las coord. X,Y del lado del dodecagano
                x1 = colector.absorbedor[j][0][i]
                y1 = colector.absorbedor[j][1][i]
                x2 = colector.absorbedor[j][0][i+1]
                y2 = colector.absorbedor[j][1][i+1]
                
                
        elif x < (N_m + (len(colector.receptor)-1))*4 + 1 + amount*(lados**id0) + 1:  
        #penultimo cuerpo a generar: el suelo
            alpha = 1
            rho = 0
            tau = 0
            x1, y1, x2, y2 = [-(W*6), -0.3, (W*6), -0.3]
        
        else:
        #ultimo cuerpo a generar: el cielo
            alpha = 1
            rho = 0
            tau = 0
            x1, y1, x2, y2 = [-(W*1.5), max_alt_cielo, (W*1.5), max_alt_cielo]
            
        surface["surf_{0}".format(x)] = rt.Superficie(alpha, rho, tau)
        surface["surf_" + str(x)].posicion(x1, y1, x2, y2)
    
    
            
    fig = plt.figure()
    for key in surface:
        plt.plot(surface[key].xplot, surface[key].yplot)
#    
        
    ### Ahora hay que hacer el raytracing
    
    #Se calcula el offset en X, cuan corrido esta el plano de proyeccion
    if theta_sol == 0:
        offset = 0
    else: 
        offset = alt_cielo/np.tan(np.pi*0.5 - np.radians(theta_sol))
    
    #Se define el numero de rayos
    N_ray = rayos
    
    #Para cada rayo
    for i in range(N_ray):
        x0 = offset + W*0.5 - i*(W/(N_ray-1))*np.cos(np.radians(theta_sol))**2
        y0 = alt_cielo + i*(W/(N_ray-1))*np.cos(np.radians(theta_sol))*np.sin(np.radians(theta_sol))
        theta0 = 1.5*np.pi - np.radians(theta_sol)
        
        #contador de reflexiones de cada rayo
        n_ref = 0
        
        while n_ref < 10: 
            #La punta del vector del rayo, de esta forma se determina la ecuacion de la recta
            x1 = x0 + 0.01*np.cos(theta0)
            y1 = y0 + 0.01*np.sin(theta0)
            
            L2 = rt.line([x0, y0], [x1, y1])
            
            #Lista con la informacion de los punto de intersección de los cuerpos
            #Se reinicia para cada rayo
            lambda_int = []
            punto_int  = []

            for key in surface:
                #Se calcula la informacion.
                #Retorna (lambda, xi, yi) si se intersecta, Falso si no. 
                info = surface[key].interseccion(L2, theta0, x0, y0)
                lambda_int.append(info[0])
                punto_int.append(info[1:])
            
            #busca el lambda con el minimo valor, es decir, el primer cuerpo que intersecta
            minimo = min(lambda_int)
            
            if minimo == 100:
                n_ref = 10
            else:
                body_int = lambda_int.index(minimo)
    #        print ('El rayo impacto en el cuerpo ' + str(body_int))
            #el cuerpo que intersecta 
            
            xi = punto_int[body_int][0]
            yi = punto_int[body_int][1]
            
            
            if xi == 10:
                print (i, theta0, minimo)
            
            #limites para el grafico, up y don van de 0.0 a 10.0
            up = 10
            down = 0
            
            #plotear trozos de los rayos
            if i < N_ray*up/10 and i > N_ray*down/10:
                plt.plot([x0, xi], [y0, yi], color = (1.0 - (i/N_ray)*0.8, 0.2, 0.3 + (i/N_ray)*0.6))
                plt.axis('equal')
                plt.xlim(-4, 4)
                plt.ylim(-0.3, max_alt_cielo)
            
            #Luego se calcula la posicion inicial y angulo del rayo
            x0, y0, theta0 = surface['surf_'+ str(body_int)].reflection(L2, theta0, xi, yi)
            
            if theta0 == 10:
                n_ref = 10
    #            print ('El rayo fue absorbido por el cuerpo' + str(body_int))
            n_ref +=1
    
    
    #Calcular el factor de interceptacion
    index_abs = (N_m + (len(colector.receptor)-1))*4 + 1
    x_y_pos=[]
    
    if id0 == 0:
        ray_interc = (surface['surf_'+str(index_abs)].ray_abs/N_ray)*np.cos(np.radians(theta_sol))*W/(N_m*w_m)
        #Numero de rayos absorbidos / Numero total de rayos. 
        #Luego se divide por el coseno del angulo del sol, ya que es la proyección paralela del sol.
        #Por ultimo se pondera por la razon entre el area total y el area de espejos. 
        
        x_y_pos.append(surface['surf_'+str(index_abs)].intercept)
        
    elif id0 == 1:
        ray_int = 0
        for i in range(amount*lados):
            index = (N_m + (len(colector.receptor)-1))*4 + 1
            ray_int += surface['surf_' + str(index+i)].ray_abs
            x_y_pos.append(surface['surf_'+str(index+i)].intercept)
            
        ray_interc = ray_int*np.cos(np.radians(theta_sol))*W/(N_m*w_m)/N_ray
#        print (ray_int)
    #En caso de tubo, hacer que cuente las 14 caras
    
    
    print ('El factor de interceptación es de ' + str(np.round(ray_interc, 4)))

    return ray_interc, surface, x_y_pos

###################################################
#IAM_t     
#angulos = np.linspace(0, 80, 9)
#pruebas = [angulos]
#for i in range(1):
#    gamma = []
#    for angulo in angulos:
#        valor, surface, x_y_flux = simulacion(angulo, 600, 20)
#        gamma.append(valor)
#        
#    pruebas.append(gamma)
#    plt.plot(angulos, gamma)
    
factor_int, surface, x_y_flux = simulacion(0, 600, 20)


##########################################################################
#Analisis de lugar de absorcion de los rayos en el absorbedor
#x_pos=np.array([])
#y_pos=np.array([])
#for i in x_y_flux:
#    for j in i:
#        x_pos = np.append(x_pos, j[0])
#        y_pos = np.append(y_pos, j[1])
#
#y_pos-=4
#
##plt.plot(x_pos, y_pos, '*')
#
#theta = np.array([])
#for i in range(len(x_pos)):
#    theta_i = np.arctan(y_pos[i]/x_pos[i])
#    
#    if theta_i < 0:
#        if y_pos[i] < 0:
#            theta_i += (2*np.pi)  
#        else: 
#            theta_i = np.pi + theta_i
#    else:
#        if x_pos[i] < 0:
#            theta_i = np.pi + theta_i
#   
#    theta = np.append(theta, theta_i)
#    
#x_energy = 75*theta
#
#n, bins, patches = plt.hist(x_energy, bins = 20)
#
#bins = np.delete(bins, -1)
#
#plt.plot(bins, n)
#########################################################################


#########################################################################
# ANALISIS DE SENSIBILIDAD DE RAYOS
#rayos = [100, 200, 300, 380, 500, 600]
#lados = [6, 10, 16, 20]
#x, y = np.meshgrid(lados, rayos)
#f_int = np.zeros([len(rayos), len(lados)])
#
#for i in range(len(rayos)):
#    for j in range(len(lados)):
#        ray_interc, surface, x_y_pos = simulacion(0, rayos[i], lados[j])
#        f_int[i,j] = ray_interc
#        print(i,j)
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#f_int.reshape(y.shape)
#surf = ax.plot_surface(x, y, f_int, cmap = cm.coolwarm)
# 
#plt.show()
#    MEJOR RESULTADO, 20 LADOS Y 600 RAYOS
#######################################################################