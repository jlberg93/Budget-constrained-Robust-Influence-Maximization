# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 07:57:00 2023

@author: Jeremy
"""


import random
import pandas as pd
import geopandas as gp
import sys
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import csv
import statistics as stat

data=gp.read_file("population_nodes.geojson")

def ProbFunc(nodeType,c):
	#return 2*c-c*c
	
	if nodeType==1:
		return c*c
	elif nodeType==2:
		return c
	else:
		return (2-c)*c
	

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
				radj[j].append((i,j))
				rP[j,i]=1-distance[j]
			else:
				radj[j]=[(i,j)]
				rP[j,i]=1-distance[j]

numV = len(data)
numE = len(adj)

f=open('oahupop.csv','r')
gains={}
largeV=0
for line in f:
	gains[int(line.split(' ')[0])]=int(line.split(' ')[1])
	largeV=largeV+int(line.split(' ')[1])
f.close()

#radj is used for the alpha_ji's 
#add empty lists to for nodes that don't have directed edges
for i in list(range(0,numV)):
	if i not in adj.keys():
		adj[i]=[]
	if i not in radj.keys():
		radj[i] = []
print("Graph has been built")

#Discount adoption function assigner
temp=pd.read_csv("./data_files/oahu_funcs.csv",sep=",",header=None)

Type={}

for lines in temp[0]:
	line=lines.split()
	if int(line[1])==4: #carryover from previous code. Not sure this is needed
		line[1]=3
	Type[int(line[0])]=int(line[1])
	

print("Functions have been assigned")

LPSolver = SolverFactory("gurobi")
NLPSolver = SolverFactory("ipopt")

budgetlist = [1,2,3,4,5]

seedlist={}

base_dir="./results/oahu/"

'''
cormickmodel = pyo.ConcreteModel()

#create the dual max problem for greek variables
cormickmodel.i = pyo.Set(initialize = list(range(0,numV)))
cormickmodel.beta = pyo.Var(cormickmodel.i, domain=pyo.NonNegativeReals)		
indexset = pyo.Set(initialize=list(P.keys()))
cormickmodel.alpha = pyo.Var(indexset, domain=pyo.NonNegativeReals, initialize=0)
cormickmodel.y=pyo.Var(cormickmodel.i, domain=pyo.Binary)
cormickmodel.w=pyo.Var(cormickmodel.i, domain=pyo.Reals)

#create a constraint list for the large greek constraint
cormickmodel.con = pyo.ConstraintList()
for i in list(range(0,numV)):
	cormickmodel.con.add(expr=(-sum(cormickmodel.alpha[i,j] for (i,j) in adj[i]) + sum(cormickmodel.alpha[j,i] for (j,i) in radj[i]) + cormickmodel.beta[i] <= gains[i]))

cormickmodel.con2 = pyo.ConstraintList()
for i in list(range(0,numV)):
	cormickmodel.con2.add(expr=(cormickmodel.beta[i] + largeV*cormickmodel.y[i] - largeV <= cormickmodel.w[i]))
	cormickmodel.con2.add(expr=(cormickmodel.w[i] <= cormickmodel.y[i]*largeV))
	cormickmodel.con2.add(expr=(cormickmodel.w[i] <= cormickmodel.beta[i]))

cormickmodel.obj = pyo.Objective(expr=(-sum((1-P[i,j])*cormickmodel.alpha[i,j] for (i,j) in P.keys()) + sum(cormickmodel.w[i] for i in cormickmodel.i)), sense=pyo.maximize)

cormickmodel.budget=pyo.Constraint(expr=(sum(cormickmodel.y[i] for i in cormickmodel.i)<=budgetlist[0]))

for budget in budgetlist:
	seedlist[budget]=[]
	cormickmodel.budget.set_value(expr=(sum(cormickmodel.y[i] for i in cormickmodel.i)<=budget))
	
	LPSolver.solve(cormickmodel, tee=True)
	seedlist[budget].append(pyo.value(cormickmodel.obj))
	for i in list(range(numV)):
		if pyo.value(cormickmodel.y[i]>0):
			seedlist[budget].append(i)

'''
seeded_nodes={}
seeded_discounts={}
inf_list={}
obj_vals={}
           
#for greedy starting points
seeded_nodes[1]=[1]
seeded_nodes[2]=[1,193]
seeded_nodes[3]=[1,193,74]
seeded_nodes[4]=[1,193,74,22]
seeded_nodes[5]=[1,193,74,22,226]

discounts={}
for budget in budgetlist:
	for i in list(range(0,numV)):
		if i in seeded_nodes[budget]:
			discounts[budget,i]=1
		else:
			discounts[budget,i]=0

'''
#IM as starting point
for budget in budgetlist:
	seeded_nodes[budget]=[]
	seeded_discounts[budget]=[]
	budget_file = open(f"./results/oahu/Alpha=1/B={budget}.txt")
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
#creates a full discount list from imported values
#turn dictionary into 2 dimensional key

for budget in budgetlist:
	print(f"Current Budget: {budget}")
	if budget == budgetlist[0]:
		#Greek model will solve for the dual variables, while the cmodel solves for discounts
		greekmodel = pyo.ConcreteModel()

		#create the dual max problem for greek variables
		greekmodel.i = pyo.Set(initialize = list(range(0,numV)))
		greekmodel.beta = pyo.Var(greekmodel.i, domain=pyo.NonNegativeReals)		
		indexset = pyo.Set(initialize=list(P.keys()))
		greekmodel.alpha = pyo.Var(indexset, domain=pyo.NonNegativeReals, initialize=0)
		
		#create a constraint list for the large greek constraint
		greekmodel.con = pyo.ConstraintList()
		for i in list(range(0,numV)):
			greekmodel.con.add(expr=(-sum(greekmodel.alpha[i,j] for (i,j) in adj[i]) + sum(greekmodel.alpha[j,i] for (j,i) in radj[i]) + greekmodel.beta[i] <= gains[i]))
		
		cmodel = pyo.ConcreteModel()
		cmodel.c = pyo.Var(pyo.Set(initialize = list(range(0,numV))), initialize=budget/numV, bounds=(0,1), domain=pyo.NonNegativeReals)
		for i in list(range(0,numV)):
			cmodel.c[i] = discounts[budget,i]
		cmodel.con = pyo.Constraint(expr = (sum(cmodel.c[i] for i in list(range(0,numV))) <= budget))
		
		pimodel = pyo.ConcreteModel()
		pimodel.pi = pyo.Var(greekmodel.i,domain=pyo.NonNegativeReals, bounds=(0,1))
		
		pimodel.con1=pyo.ConstraintList()
		for i,j in indexset:
			pimodel.con1.add(expr=pimodel.pi[i] - pimodel.pi[j] <= 1 - P[i,j])
		
		pimodel.con2=pyo.ConstraintList()
		for i in list(range(0,numV)):
			pimodel.con2.add(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],cmodel.c[i].value))
		
		greekmodel.obj = pyo.Objective(expr=(-sum((1-P[i,j])*greekmodel.alpha[i,j] for (i,j) in P.keys()) + sum(ProbFunc(Type[i],cmodel.c[i].value)*greekmodel.beta[i] for i in greekmodel.i)), sense=pyo.maximize)
		pimodel.obj=pyo.Objective(expr=sum(pimodel.pi[i] for i in list(range(0,numV))), sense=pyo.minimize)
		
		LPSolver.solve(greekmodel,tee=False)
		print(f'first greekmodel solved: {pyo.value(greekmodel.obj)}')
		
		cmodel.obj = pyo.Objective(expr=(-sum((1-P[i,j])*greekmodel.alpha[i,j].value for (i,j) in P.keys()) + sum(ProbFunc(Type[i],cmodel.c[i])*greekmodel.beta[i].value for i in greekmodel.i)), sense=pyo.maximize)
		
		NLPSolver.options['warm_start_init_point'] = 'yes'
		NLPSolver.options['warm_start_bound_push'] = 1e-6
		NLPSolver.options['warm_start_mult_bound_push'] = 1e-6
		NLPSolver.options['mu_init'] = 1e-6
		cmodel.ipopt_zL_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
		cmodel.ipopt_zU_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
		# Ipopt bound multipliers (sent to solver)
		cmodel.ipopt_zL_in = pyo.Suffix(direction=pyo.Suffix.EXPORT)
		cmodel.ipopt_zU_in = pyo.Suffix(direction=pyo.Suffix.EXPORT)
		# Obtain dual solutions from first solve and send to warm start
		cmodel.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT_EXPORT)
		cmodel.ipopt_zL_in.update(cmodel.ipopt_zL_out)
		cmodel.ipopt_zU_in.update(cmodel.ipopt_zU_out)
        
	else:
		for i in list(range(0,numV)):
			cmodel.c[i] = discounts[budget,i]
		
		cmodel.con.set_value(expr = (sum(cmodel.c[i] for i in list(range(0,numV))) <= budget))
		
		greekmodel.obj.set_value(expr=(-sum((1-P[i,j])*greekmodel.alpha[i,j] for (i,j) in P.keys()) + sum(ProbFunc(Type[i],cmodel.c[i].value)*greekmodel.beta[i] for i in greekmodel.i)))
		LPSolver.solve(greekmodel, tee=False)
		print(f'first greekmodel solved: {pyo.value(greekmodel.obj)}')
		
		cmodel.obj.set_value(expr=(-sum((1-P[i,j])*greekmodel.alpha[i,j].value for (i,j) in P.keys()) + sum(ProbFunc(Type[i],cmodel.c[i])*greekmodel.beta[i].value for i in greekmodel.i)))
		
		NLPSolver.options['warm_start_init_point'] = 'yes'
		NLPSolver.options['warm_start_bound_push'] = 1e-6
		NLPSolver.options['warm_start_mult_bound_push'] = 1e-6
		NLPSolver.options['mu_init'] = 1e-6
		cmodel.ipopt_zL_out.update(pyo.Suffix(direction=pyo.Suffix.IMPORT))
		cmodel.ipopt_zU_out.update(pyo.Suffix(direction=pyo.Suffix.IMPORT))
		# Ipopt bound multipliers (sent to solver)
		cmodel.ipopt_zL_in.update(pyo.Suffix(direction=pyo.Suffix.EXPORT))
		cmodel.ipopt_zU_in.update(pyo.Suffix(direction=pyo.Suffix.EXPORT))
		# Obtain dual solutions from first solve and send to warm start
		cmodel.dual.update(pyo.Suffix(direction=pyo.Suffix.IMPORT_EXPORT))
		cmodel.ipopt_zL_in.update(cmodel.ipopt_zL_out)
		cmodel.ipopt_zU_in.update(cmodel.ipopt_zU_out)
        
	
	NLPSolver.solve(cmodel,tee=False)
	print('discounts decided')
	
	#save the old objective value for comparison to next iteration
	old_objval = pyo.value(cmodel.obj)
	print(old_objval)
	
	diff = 100
	cutoff = 1e-3
	max_iterations=50
	current_iteration=0
	
	while (diff > cutoff or current_iteration > max_iterations):	
		current_iteration = current_iteration + 1
		print(f"==========NLP Iteration {current_iteration}==========")
		
		greekmodel.obj.set_value(expr=(-sum((1-P[i,j])*greekmodel.alpha[i,j] for (i,j) in P.keys()) + sum(ProbFunc(Type[i],cmodel.c[i].value)*greekmodel.beta[i] for i in greekmodel.i)))
		LPSolver.solve(greekmodel, tee=False)
		print('greekmodel solved')
		print(pyo.value(greekmodel.obj))
		
		cmodel.obj.set_value(expr=(-sum((1-P[i,j])*greekmodel.alpha[i,j].value for (i,j) in P.keys()) + sum(ProbFunc(Type[i],cmodel.c[i])*greekmodel.beta[i].value for i in greekmodel.i)))
		
		NLPSolver.options['warm_start_init_point'] = 'yes'
		NLPSolver.options['warm_start_bound_push'] = 1e-6
		NLPSolver.options['warm_start_mult_bound_push'] = 1e-6
		NLPSolver.options['mu_init'] = 1e-6
		cmodel.ipopt_zL_out.update(pyo.Suffix(direction=pyo.Suffix.IMPORT))
		cmodel.ipopt_zU_out.update(pyo.Suffix(direction=pyo.Suffix.IMPORT))
		# Ipopt bound multipliers (sent to solver)
		cmodel.ipopt_zL_in.update(pyo.Suffix(direction=pyo.Suffix.EXPORT))
		cmodel.ipopt_zU_in.update(pyo.Suffix(direction=pyo.Suffix.EXPORT))
		# Obtain dual solutions from first solve and send to warm start
		cmodel.dual.update(pyo.Suffix(direction=pyo.Suffix.IMPORT_EXPORT))
		cmodel.ipopt_zL_in.update(cmodel.ipopt_zL_out)
		cmodel.ipopt_zU_in.update(cmodel.ipopt_zU_out)
		
		NLPSolver.solve(cmodel,tee=False)
		print('discounts decided')
		
		new_objval=pyo.value(cmodel.obj)
		print(new_objval)
		diff = abs(new_objval-old_objval)
		old_objval = new_objval
	
	#update the constraints for the original pi LP for our final discount values. Constraint
	#numeration starts at 1 for some reason
	finaldiscount={}
	for i in list(range(0,numV)):
		pimodel.con2[i+1].set_value(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],cmodel.c[i].value))
		finaldiscount[i]=cmodel.c[i].value
		
	LPSolver.solve(pimodel,tee=False)
	print(f'Pi values solved. NLP Objective value: {pyo.value(pimodel.obj)}')
	
	obj_vals[budget]=pyo.value(pimodel.obj)
	'''
	f=open(output_dir + f"/B={budget} discounts.csv",'w',newline='')
	writer=csv.writer(f)
	writer.writerow(("Robust Obj val:",pyo.value(pimodel.obj)))
	for i in list(range(0,numV)):
		if cmodel.c[i].value>0:
			writer.writerow((i,cmodel.c[i].value))
	f.close()
	'''