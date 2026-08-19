import pandas as pd
import matplotlib.pyplot as plt

#collect number of particles
with open("particleInfo.csv", "r") as f:
    nP_line=f.readline()

    nP=int(nP_line.split(":")[1].strip())

df=pd.read_csv("particleInfo.csv", skiprows=3)

print("Loading Particle Trajectories")

target_particles=[]

num_chosen=0

#let user pick the particles they want to track 
while(num_chosen!=3):
    part_id=int(input(f"Chose a particle to track(1-{nP}):"))
    if part_id<=0 or part_id>nP:
        print("Chose an ID from 1-",nP)
    else:
        num_chosen+=1;
        target_particles.append(part_id)


#create slots in the dictionary that stores coordinates for the target particles
particle_paths={pid:{'x':[],'y':[]} for pid in target_particles}

##go through every row
for _, row in df.iterrows():
    #convert the particle id to a number and safely make it NaN if unsuccessful
    raw_pid=pd.to_numeric(row['Particle num'],errors='coerce')

    ##if the particle number is a number and is one of the target particles convert the number to an integer, the coordinates to number (with NaN error handling), and append it to the dictionary if the coordinates are numbers
    if not pd.isna(raw_pid) and int(raw_pid) in target_particles:
        p_id = int(raw_pid)
        px=pd.to_numeric(row['x'], errors='coerce')
        py=pd.to_numeric(row['y'], errors='coerce')

        if not (pd.isna(px) or pd.isna(py)):
            particle_paths[p_id]['x'].append(px)
            particle_paths[p_id]['y'].append(py)

###Debugging###
#print("Columns Pandas found:", df.columns.tolist())
#print("Data points collected:", {pid: len(coords['x']) for pid, coords in particle_paths.items()})


#create the graph
plt.figure(figsize=(8,6))

##plot the coordinates of every target particle
for particle,coords in particle_paths.items():
    plt.plot(coords['x'],coords['y'],label=f'Particle{particle}',linewidth=1.5)

plt.title("Sample Particle Trajectories")
plt.xlabel('X Position(m)')
plt.ylabel('Y Position(m)')
plt.legend()
plt.grid(True)

plt.show()

