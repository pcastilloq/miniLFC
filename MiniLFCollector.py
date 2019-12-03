# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:29:41 2019

@author: pablo
"""
from iapws import IAPWS97 as IAPWS
import numpy as np
import matplotlib.pyplot as plt
import Raytrace as rt


#%% Objeto Colector
class MiniLFCollector:
    
    #Se crea un objeto colector
    def __init__(self):
        pass
    
    #Input: Parametros de geometria, Output: Diseño del colector
    
    #Primero los espejos
    def posicionEspejos(self, W, w_m, N_m, a, L):  
        #Ancho total W, ancho de cada espejo w_m, numero de espejo n_m, altura a 
        W_util = W - w_m
        delta_w = W_util/(N_m-1)
        x = np.zeros(N_m)
        theta = np.zeros(N_m)
        for i in range(N_m):
            x[i]= (w_m/2 - W/2) + i*delta_w 
            if x[i] == 0:
                theta[i] = (np.pi/2) - 0.001
            else:
                theta[i]= np.arctan(a/x[i])
                if theta[i] < 0:
                    theta[i] = -theta[i]
                elif theta[i] < np.pi:
                    theta[i] = np.pi - theta[i] 
       
        self.centros = x
        self.theta = theta  #radianes
        self.height = a
        self.w_m = w_m
        self.ancho = W
        self.largo = L
        self.N_m = N_m
    
    #Luego el receptor trapezoidal
    def GeometriaReceptor(self, coord):
        #coord_receptor tiene que ser np.array 
        coord = coord/1000
        
        a1 = (coord[0][0] - coord[-1][0]*0.5)
        b1 = coord[0][1] + self.height
        coord[0] = np.array([a1,b1])

        coord[1:] = coord[1:] + coord[0]
        
        self.receptor = coord
    
    #Por ultimo el absorbedor de minicanales
    def GeometriaAbsorbedor(self, dimension, origen, id0 = 0, cantidad = 1, lados_poly = 20):
        
        #Mueve la placa al sistema cartesiano maestro del colector
        #0 si minichannel plate, 1 si copper tube
        #dimension: ancho en caso de placa, radio del tubo. dimension en mm.
        #arreglo con origen de los absorbedores
        self.id0 = id0 
        self.amount = cantidad
        self.lados = lados_poly
        self.dim_abs = dimension
        
        if id0 == 0:
            origen = origen[0]/1000
            x1 = origen[0] - dimension*0.5/1000
            y1 = origen[1] + self.height
            x2 = origen[0] + dimension*0.5/1000
            y2 = origen[1] + self.height
        
            self.absorbedor = [x1, y1, x2, y2]
            
        elif id0 == 1:
            XY_cil = []
            i = 0
            #Resolver despues para mas de una tuberia
            if len(origen) == cantidad:
                for i in range(cantidad):
                    theta = np.linspace(0,2*np.pi, (self.lados+1))
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
            

    
    ##Angulo de inclinacion de los espejos para un theta del sol        
    def anguloReflexion(self, theta_sol):
        self.theta_sol = np.radians(theta_sol)
        N = self.N_m
        j = self.theta_sol
        gamma = np.zeros(N)
        for i in range(N):
            O = self.theta[i]
            beta = O - np.pi/2 + j
            n = O - beta/2
            gamma[i] = n + np.pi/2
            
        self.inclinacion = gamma



    def GeometriaMinicanal(self, w_port, h_port, e):
        w_tot        = self.dim_abs                      #Ancho total de la placa
        self.w_port  = w_port
        self.h_port  = h_port
        self.e_mc    = e
        self.N_port  = np.round((w_tot - e)/(w_port + e))#Numero de puertos 
        self.P_port  = 2*(w_port + h_port)               #Perimetro minicanal
        self.A_trans = self.P_port*e                     #Area transversal de conduccion
        self.A_port  = w_port*h_port                     #Area transversal de conveccion
        self.D_h     = 4*self.A_port/self.P_port         #Diametro hidraulico minicanal



    def construccion(self, W, w_m, N_m, alt_col, L, coord_recep, dim_abs, origen_abs, w_port, h_port, e_mc, theta_sol = 0):
        self.posicionEspejos(W, w_m, N_m, alt_col, L)
        self.GeometriaReceptor(coord_recep)
        self.anguloReflexion(theta_sol)
        self.GeometriaAbsorbedor(dim_abs, origen_abs)
        self.GeometriaMinicanal(w_port, h_port, e_mc)
        self.dict_LFC_bodies()


    #Creacion de diccionario de objetos del colector
    ##Simulacion de raytraicing
    def dict_LFC_bodies(self, rayos=300):        
        
        W = self.ancho
        N_m = self.N_m
        theta_sol = self.theta_sol
        alt_cielo = self.height*1.1
        
        
        amount = self.amount
        id0 = self.id0
        lados = self.lados
        
        #Maxima altura emision de rayos
        max_alt_cielo = (alt_cielo + W*(np.cos(theta_sol)*np.sin(theta_sol)))*1.5
        
        self.max_alt_cielo = max_alt_cielo
        
        #Lista con la lista de las 4+1 esquinas de cada espejo
        XY_mir = []
        XY_rec = []
        for i in range(N_m):
            x =  self.centros[i]
            y = 0 
            x1 = x + self.w_m*np.cos(self.inclinacion[i])*0.5
            x2 = x - self.w_m*np.cos(self.inclinacion[i])*0.5
            y1 = y + self.w_m*np.sin(self.inclinacion[i])*0.5
            y2 = y - self.w_m*np.sin(self.inclinacion[i])*0.5
            
            #up y down definen hacia donde se "extruye" el cuerpo
            x_mir, y_mir = rt.polygon(1, [x1,y1], [x2,y2], 0.005, 'up')
            
            XY_mir.append([x_mir, y_mir])
        
        for i in range(len(self.receptor)-1):
            X1 = self.receptor[i]
            X2 = self.receptor[i+1]
            x_rec, y_rec = rt.polygon(1, X1, X2, 0.01, 'down')
            
            XY_rec.append([x_rec, y_rec])
        
        #Diccionario creador de objetos del colector
        surface = {}
        #Cada espejo tiene 4 lados, cada pared del receptor tiene 4 lados, el vidrio se considera 1-dim, 
        N_surf = (N_m + (len(self.receptor)-1))*4 + 1 + amount*(lados**id0) + 2
        
        for x in range(N_surf):
            if x < N_m*4:
            #espejos
        
                #Cuenta el espejo en la que va el index
                x_mir = x%4
                i = int(np.floor(x/4))
                
                #primera superficie es el espejo, el resto las caras.
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
        
                
            elif x < (N_m + len(self.receptor) - 1)*4:
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
                
            elif x == (N_m + len(self.receptor) - 1)*4:
                #Modificar esto si se agrega opcion de remover cubierta de vidrio
                #ultimo cuerpo del receptor: cubierta de vidrio
                alpha = 0.0
                rho = 0.0
                tau = 1.0
                
                x1 = self.receptor[-1][0]
                y1 = self.receptor[-1][1]
                x2 = self.receptor[0][0]
                y2 = self.receptor[0][1]
                
            elif x < (N_m + (len(self.receptor)-1))*4 + 1 + amount*(lados**id0):
            #ante-penultimo cuerpo a generar: absorbedor
                alpha = 1
                rho = 0.0
                tau = 0
                
                i = x - (N_m + (len(self.receptor)-1))*4 + 1
                if id0 == 0:
                    #minichannel plate
                    x1, y1, x2, y2 = self.absorbedor
                else: 
                    #tuberias
                    i = i%lados
                    j = int(np.floor(i/lados))
                    
                    #Las coord. X,Y del lado del dodecagano
                    x1 = self.absorbedor[j][0][i]
                    y1 = self.absorbedor[j][1][i]
                    x2 = self.absorbedor[j][0][i+1]
                    y2 = self.absorbedor[j][1][i+1]
                    
                    
            elif x < (N_m + (len(self.receptor)-1))*4 + 1 + amount*(lados**id0) + 1:  
            #penultimo cuerpo a generar: el suelo
                alpha = 1
                rho = 0
                tau = 0
                x1, y1, x2, y2 = [-(W), -0.3, (W), -0.3]
            
            else:
            #ultimo cuerpo a generar: el cielo
                alpha = 1
                rho = 0
                tau = 0
                x1, y1, x2, y2 = [-W*0.8, max_alt_cielo*1.2, W*0.8, max_alt_cielo*1.2]
                
            surface["surf_{0}".format(x)] = rt.Superficie(alpha, rho, tau)
            surface["surf_" + str(x)].posicion(x1, y1, x2, y2)
        
            self.surface = surface
        

    ####Comentar/Descomentar para plotear rayos####

    def plotColector(self, rayos = 300):

        fig = plt.figure()
        for key in self.surface:
            plt.plot(self.surface[key].xplot, self.surface[key].yplot)

        
        plt.show()    
        ### Ahora hay que hacer el raytracing

    def rotacionEspejos(self, theta_sol):
        
        surface = self.surface
        mir_center = self.centros
        gamma_0 = self.inclinacion        
        
        self.anguloReflexion(theta_sol)
        
        gamma_1 = self.inclinacion
        
        theta = gamma_1 - gamma_0
        
        for x in range(self.N_m*4):
            i = int(np.floor(x/4))
            theta_i = theta[i]
            
            #Coordenadas X1, Y1 y X2, Y2 de la linea de la superficie
            x01 = surface["surf_" + str(x)].xplot[0] - mir_center[i]
            y01 = surface["surf_" + str(x)].yplot[0]
            
            x02 = surface["surf_" + str(x)].xplot[1] - mir_center[i]
            y02 = surface["surf_" + str(x)].yplot[1]
            
            #nuevas coordenadas 
            x1 = x01*np.cos(theta_i) - y01*np.sin(theta_i) + mir_center[i]
            y1 = x01*np.sin(theta_i) + y01*np.cos(theta_i)
            
            x2 = x02*np.cos(theta_i) - y02*np.sin(theta_i) + mir_center[i]
            y2 = x02*np.sin(theta_i) + y02*np.cos(theta_i)
            
            surface["surf_" + str(x)].posicion(x1, y1, x2, y2)
 


    def simulacionRaytraicing(self, plot = "y", rayos=300):
        
        surface = self.surface
        
        W = self.ancho
        N_m = len(self.centros)
        w_m = self.w_m
        theta_sol = self.theta_sol
        alt_cielo = self.height*1.1
        max_alt_cielo = self.max_alt_cielo
        
        amount = self.amount
        id0 = self.id0
        lados = self.lados
        
        
        #Se calcula el offset en X, cuan corrido esta el plano de proyeccion
        if theta_sol == 0:
            offset = 0
        else: 
            offset = alt_cielo/np.tan(np.pi*0.5 - theta_sol)
        
        #Se define el numero de rayos
        N_ray = rayos
        
        #Para cada rayo
        for i in range(N_ray):
            #punto de inicio desde cada rayo
            x0 = offset + W*0.5 - i*(W/(N_ray-1))*np.cos(theta_sol)**2
            y0 = alt_cielo + i*(W/(N_ray-1))*np.cos(theta_sol)*np.sin(theta_sol)
            theta0 = 1.5*np.pi - theta_sol
            
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
                
                ## Graficar rayos ##
                #limites para el grafico, up y don van de 0.0 a 10.0
                if plot == "y":
                
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
        index_abs = (N_m + (len(self.receptor)-1))*4 + 1
        x_y_pos=[]
        
        if id0 == 0:
            ray_interc = (surface['surf_'+str(index_abs)].ray_abs/N_ray)*np.cos(theta_sol)*W/(N_m*w_m)
            #Numero de rayos absorbidos / Numero total de rayos. 
            #Luego se divide por el coseno del angulo del sol, ya que es la proyección paralela del sol.
            #Por ultimo se pondera por la razon entre el area total y el area de espejos. 
            
            x_y_pos.append(surface['surf_'+str(index_abs)].intercept)
            
        elif id0 == 1:
            ray_int = 0
            for i in range(amount*lados):
                index = (N_m + (len(self.receptor)-1))*4 + 1
                ray_int += surface['surf_' + str(index+i)].ray_abs
                x_y_pos.append(surface['surf_'+str(index+i)].intercept)
                
            ray_interc = (ray_int/N_ray)*W/(N_m*w_m)*np.cos(theta_sol)
    #        print (ray_int)
        #En caso de tubo, hacer que cuente las 14 caras
        
        print ('El factor de interceptación es de ' + str(np.round(ray_interc, 4)))
        
        self.intercept_factor = ray_interc
        self.surfaces = surface
        self.x_y_pos = x_y_pos      

        if plot == "y":
            for key in self.surface:
                plt.plot(self.surface[key].xplot, self.surface[key].yplot)
                plt.axis('equal')
                plt.xlim(-4, 4)
                plt.ylim(-0.3, max_alt_cielo*1.4)
                



####____Funciones Termicas_____####

## pasar asignacion de theta sol y rotacion de espejos a funcion
        
    def CondInicial(self, DNI, v_wind, T_amb, T_in, P_in, m_in):
        #Propiedades termodinamicas del fluido = Agua. Prueba
        est1 = IAPWS(T=T_in, P=P_in)  #Estado inicial. h = enthalpy, rho = density, mu = viscosity, k= conductivity
        
        self.estado_inicial = est1
        
        #Variables de entrada
        self.T_in = T_in
        self.T_amb = T_amb
        self.DNI = DNI
        self.P_in = P_in
        self.m_in = m_in
        
        #Variables calculadas
        self.Re = m_in*self.D_h/(self.A_port*self.N_port*est1.mu)      #Flujo masico de entrada o Numero de Reynolds
        self.G = self.Re*est1.mu/self.D_h
        self.x_0 = est1.x 
        self.h_0 = est1.h
        self.h_wind = 2.8 + 3*v_wind                         #Coef Transf Viento
    
         #Velocidad Inicial, Flujo Masico, Gasto masico, titulo inicial, entalpia inicial
     
    
    #Se define la funcion que calcula el coef de transferencia de calor        
    def CoefTrans(self, Pr, k, x, rho_v, rho_l, q_pto, G, h_lv, P_in, corr = "gungar"):
        
        w = self.w_port
        h = self.h_port
        D = self.D_h
        
        a=np.minimum(w/h, h/w)
        
        #Correlacion para laminar, flujo uniforme por 4 paredes, single-phase       
        b0 = 1
        b1 = -2.0421*a
        b2 = 3.0853*np.float_power(a, 2)
        b3 = -2.4765*np.float_power(a, 3)
        b4 = 1.0578*np.float_power(a, 4)
        b5 = -0.1861*np.float_power(a, 5)            
                
        Nu=8.235*(b0+b1+b2+b3+b4+b5)           
        h_lam = Nu*k/D
    
        self.h_trans = h_lam
         
        Re = self.Re
        
        
        ##Correlacion para turbulento, flujo uniforme por 4 paredes, single-phase
        if Re > 3000:
            f = np.float_power(0.79*np.log(Re)-1.64, -2)
            
            #Si es menor a 1e4 pero mayor a 3000, se resta 1000
            if Re < 10000:    
                Re = Re-1000
                
            Nu = (f/8)*(Re)*Pr/(1+12.7*np.float_power((f/8), 0.5)*(np.float_power(Pr, 2/3)-1))
            h_tur = Nu*k/D
            self.h_trans = h_tur
                      
        ## si es mayor a 1600 pero menor a 3000, es de transicion
        elif Re > 1600:
            
            Re_tr = 3000
            
            f = np.float_power(1.82*np.log10(Re_tr)-1.64, -2)
            Nu = (f/8)*(Re_tr - 1000)*Pr/(1+12.7*np.float_power((f/8), 0.5)*(np.float_power(Pr, 2/3)-1))
            h_tur = Nu*k/D
            
            self.h_trans = ((h_tur - h_lam)/1400)*(Re-1600) + h_lam
        
        
    #two-phase flow    
        if x < 0.001:
            pass
        elif x < 0.99:
                h_lo = self.h_trans
                #Razon entre viscosidad liquido y gas
                mu_l_v = IAPWS(P=P_in, x=0).mu/IAPWS(P=P_in, x=1).mu
                
                #Numero de martinelli y boiling
                Xtt = (((1-x)/x)**0.9)*((rho_v/rho_l)**0.5)*(mu_l_v**0.1)
                Bo = q_pto/(G*h_lv) 
                
                
                #Correlacion de Gungar y Winterton
                E = 1 + 24000*Bo**1.16 + 1.37*(1/Xtt)**0.86
                
                S = 1/(1 + 1.15*1e-6*(E**2)*Re**1.17)
                
                Pr = P_in/22.064            #Presion reducida = Presion / Presion critica = 22.06 MPa
                M = 18                      #peso molecular del agua wtf xd cooper qlo
                
                h_pool = 0.55*(q_pto**0.67)*(Pr**0.12)*((-1*np.log10(Pr))**-0.55)*(M**-0.55)
                
                h_gungar = E*h_lo + S*h_pool
                
                #Correlaciones de Kandlikar
                Co = np.float_power(((1-x)/x), 0.8)*np.float_power((rho_v/rho_l), 0.58)
                
                h_kand = (0.6683*np.float_power(Co, -0.2) + 1058*np.float_power(Bo, 0.7))*(np.float_power((1-x), 0.8)*h_lo)

                if corr == "gungar":
                    self.h_trans = h_gungar
                    
                elif corr == "kandlikar":
                    self.h_trans = h_kand

       
        return self.h_trans     

    #Numeero de Rayleigh
    def Ra(self, T_h, T_c, L): #Numero de Rayleigh
        T_m = (T_h + T_c)/2.0
        g = 9.8
        nu = 0.00001796         #viscosidad de cinematica. A 50ºC
        a = 0.00002546         #difusividad termica. A 50ºC 
        Ra = g*(T_h - T_c)*np.float_power(L,3)/(nu*a*T_m)
        return Ra

#Coef de transf de calor conveccion natural dentro de la cavidad
    def CoefTransAir(self, Ra, k, L):
        k = 0.028
        if Ra < 10000:
            Nu = 1.95
        elif Ra < 400000:
            Nu = 0.195*np.float_power(Ra,0.25)
        elif Ra < 10^7:
            Nu = 0.068*np.float_power(Ra,0.33)
        else:
            Nu=3.82 
        h=Nu*k/L
        
        return h


    def simulacion_thermal(self, corr ,N = 50):
        
        #Geometria y C.Inicial heredadas
        w_tot = self.dim_abs
        L = self.largo
        T_in = self.T_in
        T_amb = self.T_amb
        P_port = self.P_port
        N_port = self.N_port        
        T_sky = 0.5*np.float_power((T_amb - 273.15),1.5) + 273.15
        D_h = self.D_h
        m_in = self.m_in
        e_mc = self.e_mc
        A_trans = self.A_trans
        x_0 = self.x_0
        h_0 = self.h_0
        G = self.G
        Q_in_o = self.DNI*self.w_m*self.N_m*self.intercept_factor
        
        #Resolucion de discretizado
        #Cantidad de subdivisiones en el largo
        
        A_dif = w_tot*(L/N)/1000             #Area diferencial
        
        #Vector de temperaturas y variables
        T_fl = np.zeros(N+1)
        T_p_ext = np.zeros(N+1)
        T_p_int = np.zeros(N+1)
        T_cov = np.zeros(N+1)
        x = np.zeros(N+1)
        h = np.zeros(N+1)
        Q_u = np.zeros(N)
        coef_trans = np.zeros(N)
        
        #Temporal. Propiedades del Aire
        k_air =  0.03299              #conductividad
    
    
        #Propiedades del absorbedor
        k_cu = 400                  # conductividad cobre. 400 [W/mK]
        a_cu = 0.94                 #Absortividad de la tuberia. Coating
        eps_abs = 0.12              #emisividad del absorbedor. Coating.      
        
        e_cov = 0.10                #espacio entre vidrio y minicanal 
        eps_cov = 0.85              #Emisividad del vidrio
        tau_cov = 0.95             #transmisividad del vidrio
        rho_cov = 0.08              #reflectividad del vidrio
        
        trans = tau_cov*a_cu/(1-(1-a_cu)*rho_cov)           #transmision del vidrio y absorcion
        sigma = 5.6*np.float_power(10,(-8))     # Constante de Stefan-Boltzmann. [W/m2K4]
        
        #Integradores de Calor
        Q_loss_amb_t= 0
        Q_loss_air_t= 0
        Q_util_t = 0
    
        #Condicion de borde
        T_fl[0]= T_in                   #temperatura del fluido al inicio
        T_p_ext[0] = T_in               #temperatura de pared exterior
        T_p_int[0] = T_in               #temperatura de pared interior
        T_cov[0]= T_amb + (T_in-273)    #temperatura del cover (vidrio)
        x[0] = x_0
        h[0] = h_0
        
        
        #Iteracion
        for z in range(N):
            Q_in = Q_in_o*(L/N)   #Q_ingresa =  G_t * Area ext * delta Z
        #    if x[z] > 0:
        #        est_ini = IAPWS(x=x[z], P=P_in)
        #    else:
            est_ini = IAPWS(T=T_fl[z], P= self.P_in)          #Estado liquido del tramo
            h_in = h[z]
            k = est_ini.k
            Pr = est_ini.Prandt
            if x[z] == 1:
                self.Re = G*D_h/est_ini.mu
                
            h_f = h_in + (Q_in/self.m_in)/1000               #Entalpia estimada final del tramo
            est_out_o = IAPWS(P=self.P_in, h=h_f)            #Estado estimado final del tramo 1 
        
            T_a = T_amb
            h_a = IAPWS(T=T_a, x = 0).h
            T_b = est_out_o.T
            h_b = h_f
            
        #Aunque se puede resolver de forma analitica, mejor aplicar directamente la biseccion. 
            T_c = (T_b + T_a)/2.0
            h_c = (h_b + h_a)/2.0
            st_l = IAPWS(P = self.P_in, x = 0)
            st_v = IAPWS(T=T_fl[z], x = 1)
            h_lv = (st_v.h - st_l.h)*1000                                               #Temperatura media
            h_trans = self.CoefTrans(Pr, k, x[z], st_v.rho, st_l.rho, Q_in/A_dif, G, h_lv, self.P_in, corr)        
            #Coef Transf Calor Liquid Only
                            
            R_cu = (e_mc)/(k_cu*(P_port*N_port*(L/N)))                                      #Resistencia cobre
            R_fl = 1/(h_trans*P_port*(L/N)*N_port)                          #Resistencia flujo
            R_t = R_cu + R_fl                                               #Resistencia total
            
            T_p_ext_o = Q_in*trans*R_t  + T_c                         #Temperatura tuberia exterior inicial
            Ra_o = 100000
            h_air_o = self.CoefTransAir(Ra_o, k_air, e_cov)
            
            T_cov_o = (Q_in*(1-trans) + (self.h_wind*T_amb + h_air_o*T_p_ext_o)*A_dif)/((self.h_wind + h_air_o)*A_dif)
            
            while (h_b - h_a) > (0.01/N):
                #Perdida de calor del cover al ambiente
                Q_conv_amb = self.h_wind*(T_cov_o - T_amb)*A_dif
                Q_rad_amb = eps_cov*sigma*(np.float_power(T_cov_o,4) - np.float_power(T_sky,4))*A_dif
            
                #Perdida de calor del absorbedor al cover
                Q_conv_air = h_air_o*(T_p_ext_o - T_cov_o)*A_dif
                Q_rad_air = sigma*(np.float_power(T_p_ext_o,4) - np.float_power(T_cov_o,4))*A_dif/((1/eps_cov)+(1/eps_abs)-1)
                
                 #Calor que cede hacía adelante y atras
                Q_cond_z1 = (T_p_ext_o - T_p_ext[z])*k_cu*A_trans*N_port/(L/N)
                
                if z == 0:
                    Q_cond_z0 = 0
                else:
                    Q_cond_z0 = (T_p_ext[z] - T_p_ext[z-1])*k_cu*A_trans*N_port/(L/N)
                    
                #Calor que absorbe el absorbedor
                Q_abs = Q_in*trans
                Q_cu_0 = Q_abs - Q_conv_air - Q_rad_air + Q_cond_z0 - Q_cond_z1 
                
                #Q_loss = Q_conv_amb + Q_rad_amb - (Q_conv_air + Q_rad_air)           
            
                h_out_1 = h[z] + (Q_cu_0/m_in)/1000
                
                if h_out_1 - h_c == 0: 
                    T_fl[z+1] = T_c
                    h[z+1] = h_c
                    x[z+1] = IAPWS(h=h_c, P= self.P_in).x
                    
                elif  h_out_1 - h_c > 0:     #Cover y Absorbedor estan "helados". Luego solucion debe estar a mas Tº  
                    h_a = h_c                   #modificar esto 
                    
                else: 
                    h_b = h_c                       #Cover y Absorbedor estan "caliente". Luego solucion debe estar a menos Tº
    
                h_c = (h_a + h_b)/2.0
                est_out_2= IAPWS(h = h_c, P = self.P_in)
                T_c = est_out_2.T
                x_f_2 = est_out_2.x
        #        print('La calidad del vapor es', x_f_2)
                st_l = IAPWS(P=self.P_in, x = 0)
                st_v = IAPWS(T=T_c, x = 1)
                h_lv = (st_v.h - st_l.h)*1000
                h_trans = self.CoefTrans(Pr, k, x[z], st_v.rho, st_l.rho, Q_in/A_dif, G, h_lv, st_l.sigma, corr)
            
        #            print ('El coef transfer es ',h_kan)
                R_fl = 1/(h_trans*P_port*(L/N)*N_port)              #Resistencia flujo
                R_t = R_cu + R_fl                                   #Resistencia total
                    
                T_p_ext_o = Q_cu_0*R_t  + T_c
        #        print ('La temperatura externa es ', T_p_ext_o)
                T_p_int_o = T_p_ext_o - Q_cu_0*R_cu
                Ra_1 = abs(self.Ra(T_p_ext_o, T_cov_o, e_cov))
                h_air_o = self.CoefTransAir(Ra_1, k_air, e_cov)
        
                h_rad_amb = eps_cov*sigma*(np.float_power(T_cov_o,2) + np.float_power(T_sky,2))*(T_cov_o + T_sky)
                h_rad_air = sigma*(np.float_power(T_p_ext_o,2) + np.float_power(T_cov_o,2))*(T_p_ext_o + T_cov_o)/((1/eps_cov)+(1/eps_abs)-1)
                T_cov_o = (Q_in*(1-trans) + ((self.h_wind+h_rad_amb)*T_amb + (h_air_o+h_rad_air)*T_p_ext_o)*A_dif)/((self.h_wind + h_air_o + h_rad_amb + h_rad_air)*A_dif)
        #        print ('La temperatura del vidrio es ', T_cov_o)
            Q_u[z] = Q_cu_0    
            T_fl[z+1] = T_c
            T_p_ext[z+1]=T_p_ext_o
            T_p_int[z+1]=T_p_int_o
            T_cov[z+1] = T_cov_o
            x[z+1] = x_f_2
            h[z+1] = h_c
            coef_trans[z] = h_trans
            
            Q_loss_amb_t= Q_loss_amb_t + Q_conv_amb + Q_rad_amb
            Q_loss_air_t= Q_loss_air_t + (Q_conv_air + Q_rad_air)
            Q_loss_tot = Q_loss_amb_t + Q_loss_air_t
            Q_util_t = Q_util_t + Q_cu_0
                   
            
        T_fl=np.add(T_fl, -273.15)
        T_p_ext=np.add(T_p_ext, -273.15)
        T_p_int=np.add(T_p_int, -273.15)
        T_cov=np.add(T_cov, -273.15)
        
        Q_util = m_in*(h[N] - h[0])*1000
        eficiencia = Q_util/(Q_in_o*L)
        print ('La eficiencia del concentrador es de', np.round((eficiencia*100),2), '%')
        
        return eficiencia, T_fl, x, coef_trans, h


    def simulacion(self, theta_sol, DNI, v_wind, T_amb, T_in, P_in, m_in, plot = "y", corr="gungar"):
        #Se ajustan los espejos
        self.rotacionEspejos(theta_sol)
        #Luego se simulan los rayos
        self.simulacionRaytraicing(plot)
        #Se setean las condiciones iniciales
        self.CondInicial(DNI, v_wind, T_amb, T_in, P_in, m_in)
        #por ultimo se simula la parte termica
        eficiencia, T_fl, x, coef_trans, h = self.simulacion_thermal(corr)

        return eficiencia, T_fl, x, coef_trans, h
