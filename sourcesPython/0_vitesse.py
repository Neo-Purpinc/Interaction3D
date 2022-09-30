from Vector import *


f = open("animations/0_vitesse.pv","w")

time_step = 0.04 
nb_points = 1
nb_frames = 75 
l_vitesse = 8./3.
vecteur_vitesse = Vector(l_vitesse, 0, 0)

# ecriture de l'entete
f.write("#PV==\n" + str(time_step) + "\n" + str(nb_points) + "\n" + str(nb_frames) + "\n")
f.write("Point 0:0\n")
f.write("====\n")


#-------------- animation
position_depart = Vector(-5,0,0)

# parcours des frames (images)
for frame in range(nb_frames):
    position_depart += vecteur_vitesse * time_step
    (x,y,z)=position_depart.get()
    f.write(str(x) + " " + str(y) + " " + str(z) + "\n")

f.close()

