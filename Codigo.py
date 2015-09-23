

#Aqui resolveremos los puntos de la tarea


import time
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from astropy import units as un
from scipy import integrate






#PRIMERA PARTE

#Cargamos los datos y definimos arreglos para la Longitud y el Flujo
Datos = np.loadtxt("sun_AM0.dat")
Longitud = Datos[:,0]
Flujo = Datos[:,1]


#declaramos las unidades y las convertimos

UniLongitud = Longitud*un.nm   #longitud en nanometros
UniFlujo = Flujo*un.W*(un.m**-2)*(un.nm**-1)   #Flujo en watts, dividido metros cuadrados, dividido nanometros

Longitud_um= UniLongitud.to('um')    #convertimos a [um]
Flujo_cgs= UniFlujo.to('erg/(s cm2 um)')    #convertimos a [cgs]


plt.clf()
plt.plot(Longitud_um, Flujo_cgs)
plt.xlim(0,8)
plt.xlabel('Longitud de onda [$ \mu m $]')
plt.ylabel('Flujo de Energia [$ erg / s * cm^2 * \mu m$]')
plt.title('Grafico de Flujo de Energia en relacion con la Longitud de onda incidente')
plt.savefig('Grafico1.png', bbox_inches='tight')
#plt.show()















#Segunda Parte,  Tenemos que integrar la funcion anterior


#lo haremos usando el metodo del trapecio visto enclases

m = len(Longitud_um)
n = len(Flujo_cgs)
CSolar=0


Ttrapecio=time.time()       #Enserramos la funcion, con contadores para ver cuanto se demora
for i in range(n-1):
    paso = (Longitud_um[i+1] - Longitud_um[i])
    trapecio = ((Flujo_cgs[i+1] + Flujo_cgs[i]) * paso /2)
    CSolar += trapecio  # El area calculada corresponde a la constante Solar

Ttrapecio = time.time()-Ttrapecio   #asignamos a una variable el tiempo que demora nuestro metodo


#En paralelo usamos el metodo de python para calcular la misma integral

TcompTrapecio = time.time()   #Iniciamos el contador para le metodo de python
ConstanteComparacionT1 = np.trapz(Flujo_cgs , Longitud_um)
TcompTrapecio = time.time()-TcompTrapecio   #Cerramos el contador


print 'constantes solares con el metodo del trapecio propio y el de python, respectivamente'
print (CSolar)
print(ConstanteComparacionT1)

















#TERCERA PARTE

#Buscamos calcular el flujo energetico del sol atravez de una unidad de superficie en la atmosfera solar en una unidad de tiempo


CantiInter = input('Indique la cantidad de intervalos para la integracion (maximo 100)')
Salto = (np.pi/2-0.01)/CantiInter     #Me indica la distancia entre los cuadros a integrar
Intervalo = np.arange(0.01, np.pi/2, Salto)   #Intervalo discreto a integrar,no se puede partir de 0 asi que elejimos 0,01
Paso = Intervalo[1] - Intervalo[0]  #la distancia entre los elementos del intervalo es la misma asi que usamos un salto cualquiera.
AreaS = 0
T = 5778*un.K
Constantes = ((2*np.pi*const.h)/((const.c)**2)) * ((const.k_B*T)/(const.h))**4  #constantes que acompañan la integral
Tamano = len(Intervalo)



TSimpson=time.time()   #Iniciamos el contador para el metodo de simpson

def Integral(y):  #Definimos el argumento de la integral como una funcion para simplificar los calculos
    funcion = (np.tan(y)**3 + np.tan(y)**5) / ((np.exp(np.tan(y)))-1)

    return funcion


#Ahora iteramos para integrar a los largo de los elementos [k] del intervalo evaluados en la funcion

for k in range(0, (Tamano-2)):

    simpson = (Paso/6.0)*((Integral(Intervalo[k])) + 4*Integral(Intervalo[k+1]) + Integral(Intervalo[k+2]))
    AreaS += simpson


FlujoSolar = Constantes*AreaS    # las constantes por el area calculada (integral)
TSimpson= time.time() - TSimpson   #Cerramos el contador


#Y ahora usamos el metodo de comparacion Quad de python


TCompSimpson=time.time()    #iniciamos el contador
FlujoCom = integrate.quad(Integral, 0, np.pi/2)
FlujoCom = FlujoCom * Constantes
TCompSimpson=time.time() - TCompSimpson     #cerramos el contador

print 'flujos solar'
print FlujoSolar
print FlujoCom


#Ahora calculamos el radio del sol en base al flujo de energia en una seccion de la atmosfera terrestre y la constante a0.

a0= const.au
CSolar= CSolar.to('J /(m2 s)')  # cambio de unidades de la constante Solar, para que calzen
#El radio esta dado por la raiz cuadrada de la relacion entre la constante solar y el flujo, multiplicada por la constante a0
Radio = (np.sqrt((CSolar / FlujoSolar)))*a0
print 'Radio'
print Radio


print ' Tiempo que demoran las integracion que realizamos, con sus respectivas comparaciones'
print 'Constante Solar(metodo del trapecio)'
print Ttrapecio
print TcompTrapecio
print 'Flujo Solar (Metodo de Simpson y Quad)'
print TSimpson
print TCompSimpson
