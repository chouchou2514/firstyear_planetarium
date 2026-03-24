import math
import matplotlib.pyplot as plt
import numpy as np

#Création des dictionnaires
planete={
    'la Terre':[(5.9722)*pow(10,24),(147098291)*pow(10,3),0.0167,149598023*pow(10,3)],
    'Mercure':[(0.330)*pow(10,24),46001009*pow(10,3),0.21,57909050*pow(10,3)],
    'Venus':[(4.87)*pow(10,24),107476170*pow(10,3),0.00678,108208475*pow(10,3)],
    'Mars' : [(6.39)*pow(10,23) , (206655000)*pow(10,3) , 0.09341233 , (227944000)*pow(10,3)]}

vit_perihelie={'la Terre':[],'Mercure':[],'Venus':[]} #en m/s

#Définition des constantes
G=6.67*pow(10,-11)               #constante gravitationnelle en m^3/kg/s^2
M=1.98892*pow(10,30)         #masse du soleil en kg

#Calcul de la vitesse des planètes à la périhélie:
for i in planete.keys():
    v_perihelie=math.sqrt(G*M*(1+planete[i][2])/planete[i][3]*(1-planete[i][2]))
    vit_perihelie[i]=[v_perihelie]
    print('La vitesse à la périhélie de', i ,'est', v_perihelie*pow(10,-3),'km/s')

def euler_positions(j):
    
    #Définition du nombre de valeurs calculées
    dt = 3600                                       # Pas de temps en secondes
    t_ini = 0                                          # Temps de début en secondes
    t_fin= 365 * 24 * 3600*2          # Temps de fin en secondes (1 an)
    n_pas = int((t_fin - t_ini) / dt)  #nombre de valeurs a faire
    
    #Création des listes nécessaires
    t = np.linspace(t_ini, t_fin, n_pas)    #liste contenant les n_pas valeurs de t
    x = np.zeros(n_pas)     #liste des positions en x contenant initialemnt n_pas 0
    y = np.zeros(n_pas)     #liste des positions en y contenant initialemnt n_pas 0
    vx=np.zeros(n_pas)                                                        #liste  des vitesses en x 
    vy=np.zeros(n_pas)                            #liste des vitesses en y    



    # Initialisation des conditions initiales dans les listes
    t[0] = t_ini
    x[0] = planete[j][1]         # Position initiale en m (périhélie)
    y[0] = 0
    vx[0]=0                                                           #liste vide des vitesses en x 
    vy[0]=vit_perihelie[j][0]
    ep=[- G*planete[j][0]*M / abs(x[0])]
    ec=[0.5 *planete[j][0]*vy[0]**2]
    em=[ep[0]+ec[0]]

    # Boucle de calcul des positions
    for i in range(1, n_pas):
        r=(np.sqrt(x[i-1] ** 2 + y[i-1] ** 2))
        # Calcul de l'accélération
        ax = -G * M * x[i-1] / r**3 
        ay = -G * M * y[i-1] /r**3
       #on retrouve cela en projetant les équations sur les axes x et y
        
        # Mise à jour des positions
        x[i] = x[i-1] + vx[i-1] * dt
        y[i] = y[i-1] + vy[i-1] * dt
        
        # Mise à jour des vitesses
        vx[i] = vx[i-1] + ax * dt
        vy[i] = vy[i-1] + ay * dt
        
        #Calcul de l'énergie
        ec.append(0.5*planete[j][0]* (vx[i]**2 + vy[i]**2))
        ep.append(-G*planete[j][0]*M/r)
        em.append(ec[i]+ep[i])
    
    #Affichage de la trajectoire
    plt.subplot(1,2,1)
    plt.plot(x, y)
    plt.title("Trajectoire des planètes autour du Soleil")
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


#Création des trajectoires des 3 planètes
print(euler_positions('la Terre'))
print(euler_positions('Venus'))
print(euler_positions('Mercure'))
print(euler_positions('Mars'))
plt.show()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



    
