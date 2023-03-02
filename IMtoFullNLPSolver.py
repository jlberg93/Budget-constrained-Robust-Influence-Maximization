
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
	

print("Functions have been assigned")

NLPSolver = SolverFactory("ipopt")
LPSolver = SolverFactory("gurobi")

budgetlist = [10,20,30,40,50]
#budgetlist=[10]

base_dir = output_dir + "/Alpha=1"
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
	print(f"Current Budget: {budget}")
	if budget == budgetlist[0]:
		#Greek model will solve for the dual variables, while the cmodel solves for discounts
		fullmodel = pyo.ConcreteModel()

		#create the dual max problem for greek variables
		fullmodel.i = pyo.Set(initialize = list(range(0,numV)))
		fullmodel.beta = pyo.Var(fullmodel.i, domain=pyo.NonNegativeReals)		
		indexset = pyo.Set(initialize=list(P.keys()))
		fullmodel.alpha = pyo.Var(indexset, domain=pyo.NonNegativeReals, initialize=0)
		
		#create a constraint list for the large greek constraint
		fullmodel.con = pyo.ConstraintList()
		for i in list(range(0,numV)):
			fullmodel.con.add(expr=(-sum(fullmodel.alpha[i,j] for (i,j) in adj[i]) + sum(fullmodel.alpha[j,i] for (j,i) in radj[i]) + fullmodel.beta[i] <= 1))
		
		fullmodel.c = pyo.Var(pyo.Set(initialize = list(range(0,numV))), initialize=budget/numV, bounds=(0,1), domain=pyo.NonNegativeReals)
		for i in list(range(0,numV)):
			fullmodel.c[i] = discounts[budget,i]
		fullmodel.con2 = pyo.Constraint(expr = (sum(fullmodel.c[i] for i in list(range(0,numV))) <= budget))
		
		pimodel = pyo.ConcreteModel()
		pimodel.pi = pyo.Var(fullmodel.i,domain=pyo.NonNegativeReals, bounds=(0,1))
		
		pimodel.con1=pyo.ConstraintList()
		for i,j in indexset:
			pimodel.con1.add(expr=pimodel.pi[i] - pimodel.pi[j] <= 1 - P[i,j])
		
		pimodel.con2=pyo.ConstraintList()
		for i in list(range(0,numV)):
			pimodel.con2.add(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],fullmodel.c[i].value))
		
		fullmodel.obj=pyo.Objective(expr=(-sum((1-P[i,j])*fullmodel.alpha[i,j] for (i,j) in P.keys()) + sum(ProbFunc(Type[i],fullmodel.c[i])*fullmodel.beta[i] for i in fullmodel.i)),sense=pyo.maximize)
		pimodel.obj=pyo.Objective(expr=sum(pimodel.pi[i] for i in list(range(0,numV))), sense=pyo.minimize)
		
		NLPSolver.options['warm_start_init_point'] = 'yes'
		NLPSolver.options['warm_start_bound_push'] = 1e-6
		NLPSolver.options['warm_start_mult_bound_push'] = 1e-6
		NLPSolver.options['mu_init'] = 1e-6
		fullmodel.ipopt_zL_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
		fullmodel.ipopt_zU_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
		# Ipopt bound multipliers (sent to solver)
		fullmodel.ipopt_zL_in = pyo.Suffix(direction=pyo.Suffix.EXPORT)
		fullmodel.ipopt_zU_in = pyo.Suffix(direction=pyo.Suffix.EXPORT)
		# Obtain dual solutions from first solve and send to warm start
		fullmodel.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT_EXPORT)
		fullmodel.ipopt_zL_in.update(fullmodel.ipopt_zL_out)
		fullmodel.ipopt_zU_in.update(fullmodel.ipopt_zU_out)
		
		
	else:
		for i in list(range(0,numV)):
			fullmodel.c[i] = discounts[budget,i]
		
		fullmodel.con2.set_value(expr = (sum(fullmodel.c[i] for i in list(range(0,numV))) <= budget))
		
		fullmodel.obj.set_value(expr=(-sum((1-P[i,j])*fullmodel.alpha[i,j] for (i,j) in P.keys()) + sum(ProbFunc(Type[i],fullmodel.c[i])*fullmodel.beta[i] for i in fullmodel.i)))
		
		NLPSolver.options['warm_start_init_point'] = 'yes'
		NLPSolver.options['warm_start_bound_push'] = 1e-6
		NLPSolver.options['warm_start_mult_bound_push'] = 1e-6
		NLPSolver.options['mu_init'] = 1e-6
		fullmodel.ipopt_zL_out.update(pyo.Suffix(direction=pyo.Suffix.IMPORT))
		fullmodel.ipopt_zU_out.update(pyo.Suffix(direction=pyo.Suffix.IMPORT))
		# Ipopt bound multipliers (sent to solver)
		fullmodel.ipopt_zL_in.update(pyo.Suffix(direction=pyo.Suffix.EXPORT))
		fullmodel.ipopt_zU_in.update(pyo.Suffix(direction=pyo.Suffix.EXPORT))
		# Obtain dual solutions from first solve and send to warm start
		fullmodel.dual.update(pyo.Suffix(direction=pyo.Suffix.IMPORT_EXPORT))
		fullmodel.ipopt_zL_in.update(fullmodel.ipopt_zL_out)
		fullmodel.ipopt_zU_in.update(fullmodel.ipopt_zU_out)
		
	
	NLPSolver.solve(fullmodel,tee=True)
	
	#update the constraints for the original pi LP for our final discount values. Constraint
	#numeration starts at 1 for some reason
	finaldiscount={}
	for i in list(range(0,numV)):
		pimodel.con2[i+1].set_value(expr= -pimodel.pi[i]<= -ProbFunc(Type[i],fullmodel.c[i].value))
		finaldiscount[i]=fullmodel.c[i].value
		
	LPSolver.solve(pimodel,tee=False)
	print(f'Pi values solved. NLP Objective value: {pyo.value(fullmodel.obj)}')
	print(f'LP Objective value: {pyo.value(pimodel.obj)}')
	
	obj_vals[budget]=pyo.value(pimodel.obj)
	'''
	f=open(output_dir + f"/B={budget} discounts.csv",'w',newline='')
	writer=csv.writer(f)
	writer.writerow(("Robust Obj val:",pyo.value(pimodel.obj)))
	for i in list(range(0,numV)):
		if fullmodel.c[i].value>0:
			writer.writerow((i,fullmodel.c[i].value))
	f.close()
	'''