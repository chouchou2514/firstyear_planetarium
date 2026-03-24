import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

#Création des dictionnaires
planete={
    'la Terre':[(5.9722)*pow(10,24),(147098291)*pow(10,3),0.0167,149598023*pow(10,3)], #[masse, périhélie, excentricité (e), demi grand-axe (a)]
    'Mercure':[(0.330)*pow(10,24),46001009*pow(10,3),0.21,57909050*pow(10,3)],
    'Venus':[(4.87)*pow(10,24),107476170*pow(10,3),0.00678,108208475*pow(10,3)],
     'Mars' : [(6.39)*pow(10,23) , (206655000)*pow(10,3) , 0.09341233 , (227944000)*pow(10,3)]}

vit_perihelie={'la Terre':[],'Mercure':[],'Venus':[],'Mars':[],'Jupiter':[]} # dictionnaire des vitesses à la périhélie en m/s

#Définition des constantes
G=6.67*pow(10,-11)               #constante gravitationnelle en m^3/kg/s^2
M=1.98892*pow(10,30)         #masse du soleil en kg

#Calcul de la vitesse des planètes à la périhélie:
for i in planete.keys():
    v_perihelie=math.sqrt(G*M*(1+planete[i][2])/planete[i][3]*(1-planete[i][2]))    #racine [GM*(1+e)/a*(1-e)]
    vit_perihelie[i]=[v_perihelie]                                                                                               #on rentre les valeurs de chaque vitesse dans le dictionnaire

#Calcul des positions
def euler_asym_positions(j): #avec j la planète, elle retourne les positions x et y 

    #Définition du nombre de valeurs calculées
    dt = 3600*12                                      # Pas de temps en secondes
    t_ini = 0                                                 # Temps de début en secondes
    t_fin= 365 * 24 * 3600*2                      # Temps de fin en secondes (1 an)
    n_pas = int((t_fin - t_ini) / dt)       # nombre de valeurs à faire
    
    #Création des listes nécessaires
    t = np.linspace(t_ini, t_fin, n_pas)     #liste contenant les n_pas valeurs de t
    x = np.zeros(n_pas)                                     #liste des positions en x contenant initialemnt n_pas 0
    y = np.zeros(n_pas)                                     #liste des positions en y contenant initialemnt n_pas 0
    vx = np.zeros(n_pas)                                   #liste des vitesses en x contenant initialemnt n_pas 0
    vy = np.zeros(n_pas)                                   #liste des vitesses en y contenant initialemnt n_pas 0

    # Initialisation des conditions initiales dans les listes
    x[0] = planete[j][1]               # Position initiale en x (périhélie) en m
    y[0] = 0                                   # Position initiale en y en m
    vx[0]=0                                   # Vitesse initiale en x en m/s
    vy[0]=vit_perihelie[j][0]     # Vitesse initiale en y en m/s
    
    # Boucle de calcul des positions avec euler asymétrique
    for i in range(1, n_pas):
        # Calcul des positions à tn+1 à partir des vitesses à tn
        x[i]=x[i-1]+vx[i-1]*dt
        y[i]=y[i-1]+vy[i-1]*dt
       
       #Calcul de r en réalité norme de r
        r= np.sqrt(x[i] ** 2 + y[i] ** 2)
       # Calcul des accélérations à partir des nouvelles positions
        ax = -G * M * x[i] / r** 3
        ay = -G * M * y[i] / r ** 3
        
         # Calcul des vitesses au temps tn+1 à partir des accélérations
        vx[i] = vx[i-1] + ax * dt
        vy[i] = vy[i-1] + ay * dt
        
    return x, y

#Création de la figure et de l'axe
fig, ax= plt.subplots()

#Echelle du graphique
ax.set_xlim(-2.5e11, 2.5e11)           #choisi
ax.set_ylim(-2.5e11, 2.5e11)           #choisi
ax.set_aspect('equal')                #on souhaite obtenir à graphique carré

#Création des trajectoires mises à jour au fur et à mesure et des points
ligne_terre, = ax.plot([], [], 'c-', label='Terre')                       #Initialisation du tracé de la trajectoire
ligne_venus, = ax.plot([], [], 'w-', label='Venus')                    #Initialisation du tracé de la trajectoire
ligne_mercure, = ax.plot([], [], 'm-', label='Mercure')           #Initialisation du tracé de la trajectoire
ligne_mars, =ax.plot([],[],'r-', label ='Mars')
scatter_terre, = ax.plot([], [], 'co')                                             #Initialisation du point de la planète
scatter_venus, = ax.plot([], [], 'wo')                                           #Initialisation du point de la planète
scatter_mercure, = ax.plot([], [], 'mo')                                      #Initialisation du point de la planète
scatter_mars, = ax.plot([],[],'ro')
ax.legend()                                                                                         #Affichage des légendes

#Liste des trajectoires des planètes
lignes = [ligne_terre, ligne_mercure, ligne_venus,ligne_mars]

#Initialisation des lignes de trajectoire et des points avant l'animation avec coordonnées vides
def init():  # initialisation des lignes de trajectoire et des points en position vide
    for ligne in lignes:
        ligne.set_data([], [])
    scatter_terre.set_data([], [])
    scatter_mercure.set_data([], [])
    scatter_venus.set_data([], [])
    scatter_mars.set_data([],[])
    
    return lignes +[scatter_terre, scatter_mercure, scatter_venus,scatter_mars]    #retourne la liste des trajectoires et la liste des points vides

#Création de la fonction des nouvelles images
def update(frame):     #met à jour les positions des trajectoires et des points avec frame le nombre d'images
    x_terre, y_terre = euler_asym_positions('la Terre')          #appel de la fonction euler_asym_positions pour définir les positions
    x_mercure, y_mercure = euler_asym_positions('Mercure')
    x_venus, y_venus = euler_asym_positions('Venus')
    x_mars, y_mars = euler_asym_positions('Mars')

    lignes[0].set_data(x_terre[:frame], y_terre[:frame])                        #:frame va de l'indice 0 a frame pour mettre à jour, à chaque étape il y a de plus en plus d'indices un par un
    #appel de l'indice 0 de la liste 'lignes' qui est la trajectoire de la terre (ligne_terre), met à jour la trajectoire
    lignes[1].set_data(x_mercure[:frame], y_mercure[:frame])
    lignes[2].set_data(x_venus[:frame], y_venus[:frame])
    lignes[3].set_data(x_mars[:frame], y_mars[:frame])
    
    #Mise à jour des points représentants la position des planètes
    scatter_terre.set_data([x_terre[frame]], [y_terre[frame]])                 #même principe pour les points, rempli la liste des positions successives de la planète
    scatter_mercure.set_data([x_mercure[frame]], [y_mercure[frame]])
    scatter_venus.set_data([x_venus[frame]], [y_venus[frame]])
    scatter_mars.set_data([x_mars[frame]], [y_mars[frame]])

    return lignes +[scatter_terre, scatter_mercure, scatter_venus, scatter_mars]     #retourne la liste des trajectoires et la liste des points


ani = FuncAnimation(fig, update, frames=9000, init_func=init, interval = 20,blit=True)
#(figure à animer, fonction appelée à chaque étape, le nombre total d'étapes, la fonction initiale, l'intervalle entre chaque étape en ms)
#le blit=True permet à l'animation de ne pas changer ce qui a déjà été déssiné

#Création du soleil 
plt.scatter(149598023*pow(10,3) -(147098291)*pow(10,3),0,s=300,c="yellow")    #position du soleil en x = a - périhélie

#Création d'étoiles aléatoirement
plt.scatter(0.5*10**11,1.5*10**11,s=1,marker ='d',c="lemonchiffon")
plt.scatter(-0.5*10**11,-1.5*10**11,s=1,marker ='d',c="lemonchiffon")
plt.scatter(-1*10**11,1*10**11,s=1,marker ='d',c="lemonchiffon")
plt.scatter(1*10**11,-1*10**11,s=1,marker ='d',c="lemonchiffon")
plt.scatter(-1*10**11,-0.5*10**11,s=1,marker ='d',c="lemonchiffon")
plt.scatter(-1.5*10**11,-1.5*10**11,s=1,marker ='d',c="lemonchiffon")
plt.scatter(-1.5*10**11,1.5*10**11,s=1,marker ='d',c="lemonchiffon")
plt.scatter(1.5*10**11,0.5*10**11,s=1,marker ='d',c="lemonchiffon")

#Finalisation du graphique
plt.title("Trajectoires des planètes autour du Soleil")
plt.xlabel("Position x (m)")
plt.ylabel("Position y (m)")
ax.set_facecolor('midnightblue')

plt.show()