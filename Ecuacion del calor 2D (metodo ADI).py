# -*- coding: utf-8 -*-
"""
Fecha: 22/05/2023
Autor: cheng
Contacto: bigdata.cheng@gmail.com

Objetivo: resolver y graficar la difusión del calor en 2 dimensiones
"""


import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

dx = 0.1 # diferencial en x (coordenada x)
dy = 0.1 # diferencial en y (coordenada y)
dt = 0.005 # diferencia en t (tiempo)

Nx = 20; Ny = 20; dimt = 250 # son el numero de particiones para x,y,t


alphax = 0.4 
alphay = 0.4

sx = alphax*dt/(dx**2)
sy = alphay*dt/(dy**2)

#tanto alphax, alphay, sx, sy son coeficientes que aparecen cuando se trata de resolver la ecuacion del calor



#-------ESTABILIDAD:--------------------------------------------------

print('sx + sy = %f' %(sx+sy))

if (sx + sy <= 0.5):
    print('Se cumple la condicion de estabilidad')
else:
    print('No se cumple la condicion de estabilidad')


#condiciones iniciales:--------------------------------------------------
  
    
T = np.zeros(shape=(Nx+1,Ny+1))


T[8:12,8:12] = 100 # inicialmente hay una temperatura de 100ºC en el centro

Tnew = np.zeros(shape=(Nx+1,Ny+1))
    
    
   
    
   
# Por el metodo ADI hago primero implicita en direccion OX y explicito en OY, luego hago explicito en OX e implicito en OY. 
#---------------------------------------- Defino las matrices A1 y A2 que viene de (A1)T* = (A2)T^n
A1 = np.zeros(shape=(Nx+1,Nx+1)) 
A1[0][0] = 1 ; A1[Nx][Nx]=1
for i in range(1,Nx):
    A1[i][i] = 1 + sx
    A1[i][i-1] = -0.5*sx
    A1[i][i+1] = -0.5*sx
        
A1_inv = np.linalg.inv(A1) # hago la inversa de la matriz A1



A2 = np.zeros(shape=(Ny+1,Ny+1)) 
A2[0][0] = 1 ; A2[Ny][Ny]=1
for i in range(1,Ny):
    A2[i][i] = 1 - sy
    A2[i][i-1] = 0.5*sy
    A2[i][i+1] = 0.5*sy
    
#---------------------------------------- Defino las matrices A3 y A4 que viene de (A3)T^{n+1} = (A4)T*
A3 = np.zeros(shape=(Nx+1,Nx+1)) 
A3[0][0] = 1 ; A3[Nx][Nx]=1
for i in range(1,Nx):
    A3[i][i] = 1 + sy
    A3[i][i-1] = -0.5*sy
    A3[i][i+1] = -0.5*sy
        
A3_inv = np.linalg.inv(A3) # hago la inversa de la matriz A3



A4 = np.zeros(shape=(Ny+1,Ny+1)) 
A4[0][0] = 1 ; A4[Ny][Ny]=1
for i in range(1,Ny):
    A4[i][i] = 1 - sx
    A4[i][i-1] = 0.5*sx
    A4[i][i+1] = 0.5*sx
    
#-----------------------------------------------------------------------------------------------


 #Taux es una matriz para el estado intermedio entre T y Tnew


for n in range(dimt):
    
        
            
    Taux = (np.dot( np.dot(A1_inv,A2) , T.T)) # Taux es producto de A1_inv, A2 y la transpuesta de T
    
    Tnew = (np.dot( np.dot(A3_inv,A4) , Taux.T)) # Tnew es producto de A3_inv, A4 y la transpuesta de Taux
    

    
    #condiciones de frontera de Dirichlet---
        
    # Tnew[0,0:Ny] = 10
    # Tnew[Nx,0:Ny] = 0
    
    # Tnew[0:Nx,0] = 0
    # Tnew[0:Nx,Ny] = 0
    
    
    
    #condiciones de frontera: flujo nulo ---
        
    Tnew[0,0:Ny+1] = Tnew[1,0:Ny+1]
    Tnew[Nx,0:Ny+1] = Tnew[Nx-1,0:Ny+1]
    
    Tnew[0:Nx+1,0] = Tnew[0:Nx+1,1]
    Tnew[0:Nx+1,Ny] = Tnew[0:Nx+1,Ny-1]
    
    
    
    #actualizamos variables----------------------------
    
    for m in range(0,Nx+1):
        for l in range(0,Ny+1):
            T[m,l]= Tnew[m,l]
        
        
    
    
    
    
    #grafica con pcolor-----------------------------------

    if n%10 == 0:
        # plt.clf()
        

        fig = plt.figure()  #creo una figura vacia, le puedo introducir argumentos (resolucion, nombre, tamaño, ...)
        ax = fig.add_subplot(111) #1 fila, 1 columna, 1 elemento. si pongo 122 tendremos dos columnas


        # ax.minorticks_on()
        # ax.tick_params(axis='x', lengh=6, width = 1, labelsize=12)
        # ax.tick_params(which='both',direction='in',top='on',right='on')
        ax.ticklabel_format(axis='y', style='plain')
        ax.set_title('t = %.2f' %(n*dt))
        ax.set_xlabel('X'); ax.set_ylabel('Y')
    
        plt.pcolor(T) 

        plt.colorbar(label='Intensidad')

        plt.clim(0,10) #los limites fijos de la barra de color. 
        # plt.plot(T)
        plt.pause(0.01)
        plt.show()
