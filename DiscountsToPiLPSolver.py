
import random
import pandas
import sys
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import csv
import statistics as stat

#config file containing location of graph and discount function csv files
config_file = open("./data_files/email.config").readlines()
#config_file = open("./data_files/redditunif-negxsquared.config").readlines()
#config_file = open("./data_files/astrolaptop-negxsquared.config").readlines()

def ProbFunc(nodeType,c):
	#return 2*c-c*c
	
	if nodeType==1:
		return c*c
	elif nodeType==2:
		return c
	else:
		return (2-c)*c
	

#Set default values for filepaths and running variables
data_csv=""
func_csv=""
output_dir=""

#Note that we do not currently use start, end, DStep, or batches

data_counter=1
for line in config_file:
	if str(line)[0] =="#":
		next
	else:
		if data_counter==1:
			data_csv=line[:-1]
		elif data_counter==2:
			func_csv=line[:-1]
		elif data_counter==3:
			output_dir=line[:-1]
		data_counter+=1

#Graph Creation
temp=pandas.read_csv(data_csv,sep=",",header=None)

numV = int(temp[0][0].split()[0])
numE = int(temp[0][0].split()[1])

adj={}
radj={}
P={}
rP={}

for lines in temp[0]:
	line=lines.split()
	if len(line)<3:
		continue
	u=int(line[0])
	v=int(line[1])
	deg=int(line[2])
	#edge probabilites are calculated using the degree
	pp=1.0/deg
	
	if u in adj.keys():
		adj[u].append((u,v))
		P[u,v]=pp
	else:
		adj[u]=[(u,v)]
		P[u,v]=pp
	
	if v in radj.keys():
		radj[v].append((u,v))
		rP[v,u]=pp
	else:
		radj[v]=[(u,v)]
		rP[v,u]=pp

#radj is used for the alpha_ji's 
#add empty lists to for nodes that don't have directed edges
for i in list(range(0,numV)):
	if i not in adj.keys():
		adj[i]=[]
	if i not in radj.keys():
		radj[i] = []
print("Graph has been built")

#Discount adoption function assigner
temp=pandas.read_csv(func_csv,sep=",",header=None)

Type={}

for lines in temp[0]:
	line=lines.split()
	if int(line[1])==4: #carryover from previous code. Not sure this is needed
		line[1]=3
	Type[int(line[0])]=int(line[1])
	
LPSolver = SolverFactory("gurobi")

budgetlist = [10,20,30,40,50]
#budgetlist=[10]

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
#creates a full discount list from imported values
#turn dictionary into 2 dimensional key

for budget in budgetlist:
	print(f"============= Current Budget: {budget} =============")
	if budget == budgetlist[0]:
		pimodel = pyo.ConcreteModel()
		pimodel.pi = pyo.Var(pyo.Set(initialize=list(range(0,numV))),domain=pyo.NonNegativeReals, bounds=(0,1))
		
		pimodel.con1=pyo.ConstraintList()
		for i,j in list(P.keys()):
			pimodel.con1.add(expr=pimodel.pi[i] - pimodel.pi[j] <= 1 - P[i,j])
		
		pimodel.con2=pyo.ConstraintList()
		for i in list(range(0,numV)):
			pimodel.con2.add(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],discounts[budget,i]))
		
		pimodel.obj=pyo.Objective(expr=sum(pimodel.pi[i] for i in list(range(0,numV))), sense=pyo.minimize)
		LPSolver.solve(pimodel,tee=False)
		fullpi=pyo.value(pimodel.obj)
	else:
		for i in list(range(0,numV)):
			pimodel.con2[i+1].set_value(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],discounts[budget,i]))
		LPSolver.solve(pimodel, tee=False)
		fullpi=pyo.value(pimodel.obj)
	
	print(f"Budget {budget} full pi: {fullpi}")
