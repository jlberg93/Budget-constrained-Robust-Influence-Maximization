
import statistics
import random
import pandas
import sys
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import csv

#config file containing location of graph and discount function csv files
config_file = open("wiki-Vote-desktop.config").readlines()
#config_file = open("wiki-Vote.config").readlines()

#Set default values for filepaths and running variables
data_csv=""
func_csv=""
output_dir=""
batches=100
start=1
end=5
DStep=0.05

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
		elif data_counter==4:
			batches=int(line[:-1])
		elif data_counter==5:
			start=int(line[:-1])
		elif data_counter==6:
			end=int(line[:-1])
		elif data_counter==7:
			DStep=float(line)
		else:
			print("Too many arguments given in the config file. Check config file for proper format")
			sys.exit()
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

list1=[]
list2=[]
list3=[]

for lines in temp[0]:
	line=lines.split()
	if int(line[1])==4: #carryover from previous code. Not sure this is needed
		line[1]=3
	Type[int(line[0])]=int(line[1])
	if int(line[1]==1):
		list1.append(int(line[0]))
	elif int(line[1]==2):
		list2.append(int(line[0]))
	else:
		list3.append(int(line[0]))

def ProbFunc(nodeType,c):
	if nodeType==1:
		return c*c
	elif nodeType==2:
		return c
	else:
		return (2-c)*c

print("Functions have been assigned")

NLPSolver = SolverFactory("ipopt")
LPSolver = SolverFactory("gurobi")

NLPSolver.warmstart=True

#budgetlist = [10,20,30,40,50]
budgetlist=[10]

base_dir = "C:/Users/Jeremy/OneDrive - Naval Postgraduate School/01- Student/99 Thesis/Code/CIM-master-Robust/Results/Alpha=1"
#base_dir = "C:/Users/Jeremy Berg/OneDrive - Naval Postgraduate School/01- Student/99 Thesis/Code/CIM-master-Robust/Results/Alpha=1"
seeded_nodes={}
seeded_discounts={}

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
	print(f"====Starting Budget = {budget}====")
	if budget == budgetlist[0]:
		#Greek model will solve for the dual variables, while the cmodel solves for discounts
		greekmodel = pyo.ConcreteModel()

		#create the dual max problem for greek variables
		greekmodel.i = pyo.Set(initialize = list(range(0,numV)))
		greekmodel.beta = pyo.Var(greekmodel.i, domain=pyo.NonNegativeReals)		
		indexset = pyo.Set(initialize=list(P.keys()))
		greekmodel.alpha = pyo.Var(indexset, domain=pyo.NonNegativeReals)
		
		#create a constraint list for the large greek constraint
		greekmodel.con = pyo.ConstraintList()
		for i in list(range(0,numV)):
			greekmodel.con.add(expr=(-sum(greekmodel.alpha[i,j] for (i,j) in adj[i]) + sum(greekmodel.alpha[j,i] for (j,i) in radj[i]) + greekmodel.beta[i] <= 1))
		
		pimodel = pyo.ConcreteModel()
		pimodel.pi = pyo.Var(greekmodel.i,domain=pyo.NonNegativeReals, bounds=(0,1))
		
		pimodel.con1=pyo.ConstraintList()
		for i,j in indexset:
			pimodel.con1.add(expr=pimodel.pi[i] - pimodel.pi[j] <= 1 - P[i,j])
		
		pimodel.con2=pyo.ConstraintList()
		for i in list(range(0,numV)):
			pimodel.con2.add(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],discounts[budget,i]))
		
		pimodel.obj=pyo.Objective(expr=sum(pimodel.pi[i] for i in list(range(0,numV))), sense=pyo.minimize)
		LPSolver.solve(pimodel,tee=False)
		print(f"Original IM pi value: {pyo.value(pimodel.obj)}")
		
		greekmodel.obj = pyo.Objective(expr=(-sum((1-P[i,j])*greekmodel.alpha[i,j] for (i,j) in P.keys()) + sum(ProbFunc(Type[i],discounts[budget,i])*greekmodel.beta[i] for i in greekmodel.i)), sense=pyo.maximize)
		LPSolver.solve(greekmodel,tee=False)
		print(f"Original LP obj value: {pyo.value(greekmodel.obj)}")
		
		greekmodel.c = pyo.Var(pyo.Set(initialize = list(range(0,numV))), initialize=budget/numV, bounds=(0,1), domain=pyo.NonNegativeReals)
		for i in list(range(0,numV)):
			greekmodel.c[i] = discounts[budget,i]
			
		greekmodel.con2 = pyo.Constraint(expr = (sum(greekmodel.c[i] for i in list(range(0,numV))) <= budget))
		
		greekmodel.obj.set_value(expr=(-sum((1-P[i,j])*greekmodel.alpha[i,j] for (i,j) in P.keys()) + sum((greekmodel.c[i]**2)*greekmodel.beta[i] for i in list1) 
								 + sum(greekmodel.c[i]*greekmodel.beta[i] for i in list2) + sum((2*greekmodel.c[i] - greekmodel.c[i]**2)*greekmodel.beta[i] for i in list3)))
		NLPSolver.solve(greekmodel,tee=True)
		print(f'Full model solved. Objective value: {pyo.value(greekmodel.obj)}')
		
		#update the constraints for the original pi LP for our final discount values. Constraint
		#numeration starts at 1 for some reason
		for i in list(range(0,numV)):
			pimodel.con2[i+1].set_value(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],greekmodel.c[i].value))
			
		LPSolver.solve(pimodel,tee=False)
		print(f'Pi values updated. NLP Objective value: {pyo.value(pimodel.obj)}')
		
		'''
		f=open(f"IMB={budget} discounts.csv",'w',newline='')
		writer=csv.writer(f)
		for i in list(range(0,numV)):
			writer.writerow((i,greekmodel.c[i].value))
		
		f.close()
		'''
	else:
		for i in list(range(0,numV)):
			greekmodel.c[i] = discounts[budget,i]
		
		for i in list(range(0,numV)):
			pimodel.con2[i+1].set_value(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],greekmodel.c[i].value))
		
		LPSolver.solve(pimodel,tee=False)
		print(f"Original IM pi value: {pyo.value(pimodel.obj)}")
		
		greekmodel.con2.set_value(expr = (sum(greekmodel.c[i] for i in list(range(0,numV))) <= budget))
		
		NLPSolver.solve(greekmodel, tee=True)
		print(f'Full model solved. Objective value: {pyo.value(greekmodel.obj)}')
			
		for i in list(range(0,numV)):
			pimodel.con2[i+1].set_value(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],greekmodel.c[i].value))
			
		LPSolver.solve(pimodel,tee=False)
		print(f'Pi values updated. NLP Objective value: {pyo.value(pimodel.obj)}')
		
		f=open(f"IMB={budget} discounts.csv",'w',newline='')
		writer=csv.writer(f)
		for i in list(range(0,numV)):
			writer.writerow((i,greekmodel.c[i].value))
		
		f.close()