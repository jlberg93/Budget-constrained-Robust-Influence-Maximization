
import random
import pandas as pd
import geopandas as gp
import sys
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import csv
import statistics as stat

data=gp.read_file("population_nodes.geojson")

edgelist=[]
adj={}
radj={}
P={}
rP={}
for i in range(len(data)):
	distance=data.distance(data.iloc[i]['geometry'])
	distance.drop([i],inplace=True)
	distance=distance*.621371/1000 #convert from meters to miles
	distance=distance[distance<1]
	if len(distance)>0:
		for j in distance.index:
			edgelist.append((i,j,1-distance[j]))
			if i in adj.keys():
				adj[i].append((i,j))
				P[i,j]=1-distance[j]
			else:
				adj[i]=[(i,j)]
				P[i,j]=1-distance[j]
			
			if j in radj.keys():
				radj[j].append((i,i))
				rP[j,i]=1-distance[j]
			else:
				radj[j]=[(i,j)]
				rP[j,i]=1-distance[j]

numV = len(data)
numE = len(adj)

f=open('oahupop.csv','r')
gains={}
for line in f:
    gains[int(line.split(' ')[0])]=int(line.split(' ')[1])
f.close()

#radj is used for the alpha_ji's 
#add empty lists to for nodes that don't have directed edges
for i in list(range(0,numV)):
	if i not in adj.keys():
		adj[i]=[]
	if i not in radj.keys():
		radj[i] = []
print("Graph has been built")

budgetlist = [1,2,3,4,5]
#budgetlist=[10]
'''
base_dir = output_dir + "/Alpha=1"
#base_dir = "C:/Users/Jeremy Berg/OneDrive - Naval Postgraduate School/01- Student/99 Thesis/Code/CIM-master-Robust/Results/Alpha=1"
seeded_nodes={}
seeded_discounts={}
inf_list={}
obj_vals={}

for budget in budgetlist:
	seeded_nodes[budget]=[]
	seeded_discounts[budget]=[]
	budget_file = open(base_dir + f"/B={budget}.txt")
	found=False
	for line in budget_file:
		splitline = line.split(sep='\t')
		if found:
			seeded_nodes[budget].append(int(splitline[0]))
			seeded_discounts[budget].append(float(splitline[1]))
		if splitline[0]=="CD":
			found=True

discounts={}
for budget in budgetlist:
	for i in list(range(0,numV)):
		if i in seeded_nodes[budget]:
			discounts[budget,i]=seeded_discounts[budget][seeded_nodes[budget].index(i)]
		else:
			discounts[budget,i]=0
'''

### NEED TO LOAD IN THE POPULATION VALUES
#LP Values
#LPseeds=[1,193,74,22,226]

#IC Values
LPseeds=[52,83,206,280,220]

seeded_nodes={}
seeded_discounts={}
inf_list={}
obj_vals={}

for budget in budgetlist:
	print(f'Budget {budget} true performance calculating')
	
	#the below chunk only runs if we want to calculate the IC performance of the budget distribution found above
	trial=0
	im_trials = 80000
	inf_list[budget]=[]
	seeds=[]

	
	for i in list(range(budget)):
		seeds.append(LPseeds[i])
	
	while trial < im_trials:
		#seed selection for trial
		vis=[]
		for i in list(range(0,numV)):
			vis.append(False)
			
		#propogation of seed set
		influenced=0
		q = []
		for seed in seeds:
			vis[seed]=True
			q.append(seed)
			influenced=influenced + gains[seed]
		
		while len(q) > 0:
			u = q.pop(0)
			for u,v in adj[u]:
				if vis[v]==False:
					runif=random.uniform(0,1)
					if runif < P[u,v]:
						vis[v]=True
						q.append(v)
						influenced=influenced+gains[v] ##add the population values here in terms of the list
			
		inf_list[budget].append(influenced)
		trial=trial+1
		if trial%10000==0:
			print(f'Trial {trial}')
	print(f'Budget {budget} performance: {sum(inf_list[budget])/len(inf_list[budget])}')