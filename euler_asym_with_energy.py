import math
import matplotlib.pyplot as plt
import numpy as np

#Création des dictionnaires
planete={
    'la Terre':[(5.9722)*pow(10,24),(147098291)*pow(10,3),0.0167,149598023*pow(10,3)],
    'Mercure':[(0.330)*pow(10,24),46001009*pow(10,3),0.21,57909050*pow(10,3)],
    'Venus':[(4.87)*pow(10,24),107476170*pow(10,3),0.00678,108208475*pow(10,3)],
    'Mars' : [(6.39)*pow(10,23) , (206655000)*pow(10,3) , 0.09341233 , (227944000)*pow(10,3)],
  }

vit_perihelie={'la Terre':[],'Mercure':[],'Venus':[]} #en m/s

#Définition des constantes
G=6.67*pow(10,-11)               #constante gravitationnelle en m^3/kg/s^2
M=1.98892*pow(10,30)         #masse du soleil en kg

#Calcul de la vitesse des planètes à la périhélie:
for i in planete.keys():
    v_perihelie=math.sqrt(G*M*(1+planete[i][2])/planete[i][3]*(1-planete[i][2]))
    vit_perihelie[i]=[v_perihelie]
    print('La vitesse à la périhélie de', i ,'est', v_perihelie*pow(10,-3),'km/s')

def euler_asym_ener(j):
    #Définition du nombre de valeurs calculées
    dt = 3600                                       # Pas de temps en secondes
    t_ini = 0                                          # Temps de début en secondes
    t_fin= 365 * 24 * 3600*2         # Temps de fin en secondes (1 an)
    n_pas = int((t_fin - t_ini) / dt)  #nombre de valeurs a faire
    
    #Création des listes nécessaires
    t = np.linspace(t_ini, t_fin, n_pas)    #liste contenant les n_pas valeurs de t
    x = np.zeros(n_pas)     #liste des positions en x contenant initialemnt n_pas 0
    y = np.zeros(n_pas)     #liste des positions en y contenant initialemnt n_pas 0
    z = np.zeros(n_pas)
    vx=np.zeros(n_pas)                          #liste des vitesses en x contenant initialemnt n_pas 0
    vy=np.zeros(n_pas)                          #liste des vitesses en y contenant initialemnt n_pas 0
    vz=np.zeros(n_pas)

    # Initialisation des conditions initiales dans les listes
    x[0] = planete[j][1]         # Position initiale en m (périhélie)
    y[0] = 0                             # Position initiale en m
    z[0]=0
    vx[0]=0                           # Vitesse initiale en x en m/s
    vy[0]=vit_perihelie[j][0]   # Vitesse initiale en y en m/s
    vz[0]=0
    ep=[- G*planete[j][0]*M / abs(x[0])]
    ec=[0.5 *planete[j][0]*vy[0]**2]
    em=[ep[0]+ec[0]]
    
    # Boucle de calcul des positions
    for i in range(1, n_pas):
        # Calcul des positions à tn+1 à partir des vitesses à tn
        x[i]=x[i-1]+vx[i-1]*dt
        y[i]=y[i-1]+vy[i-1]*dt
        z[i]=z[i-1]+vz[i-1]*dt
       
       #Calcul de r
        r=(np.sqrt(x[i] ** 2 + y[i] ** 2+ z[i]**2))
       # Calcul des accélérations à partir des nouvelles positions
        ax = -G * M * x[i] /  r** 3
        ay = -G * M * y[i] / r ** 3
        az = -G * M * z[i] / r** 3
        
         # Calcul des vitesses au temps tn+1 à partir des accélérations
        vx[i] = vx[i-1] + ax * dt
        vy[i] = vy[i-1] + ay * dt
        vz[i] = vz[i-1] + az * dt
        
        #Calcul de l'énergie
        ec.append(0.5*planete[j][0]* (vx[i]**2 + vy[i]**2))
        ep.append(-G*planete[j][0]*M/r)
        em.append(ec[i]+ep[i])
    
    #Affichage de la trajectoire
    plt.subplot(1,2,1)
    plt.plot(x, y,label=j)
    plt.title("Trajectoire des planètes autour du Soleil")
    plt.legend()
    plt.xlabel("Position x (m)")
    plt.ylabel("Position y (m)")
    plt.axis('square')
    #Création du soleil 
    plt.scatter(149598023*pow(10,3)   -(147098291)*pow(10,3),0,s=300,c="yellow")
    
    #Affichage des courbes de Em
    plt.subplot(1,2,2)
    plt.plot(t, em, label=f"Em {j}")
    plt.title('Conservation énergie mécanique')
    plt.ylabel("Em (en J)")
    plt.xlabel("Temps (en s)")
    plt.legend()


print(euler_asym_ener('la Terre'))
print(euler_asym_ener('Venus'))
print(euler_asym_ener('Mercure'))
print(euler_asym_ener('Mars'))
plt.show()


        
        