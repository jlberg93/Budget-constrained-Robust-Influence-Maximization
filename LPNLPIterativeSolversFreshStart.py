
import statistics
import random
import pandas
import sys
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import csv
import statistics as stat

#config file containing location of graph and discount function csv files
#config_file = open("./data_files/wiki-Vote-negxsquared fresh.config").readlines()
#config_file = open("./data_files/redditunif-negxsquared.config").readlines()
config_file = open("./data_files/astro-negxsquared fresh.config").readlines()
print('negxsquared fresh start')

def ProbFunc(nodeType,c):
	return 2*c-c*c
	'''
	if nodeType==1:
		return c*c
	elif nodeType==2:
		return c
	else:
		return (2-c)*c
	'''

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

inf_list={}
obj_vals={}

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
			greekmodel.con.add(expr=(-sum(greekmodel.alpha[i,j] for (i,j) in adj[i]) + sum(greekmodel.alpha[j,i] for (j,i) in radj[i]) + greekmodel.beta[i] <= 1))
		
		cmodel = pyo.ConcreteModel()
		cmodel.c = pyo.Var(pyo.Set(initialize = list(range(0,numV))), initialize=0, bounds=(0,1), domain=pyo.NonNegativeReals)
		cmodel.con = pyo.Constraint(expr = (sum(cmodel.c[i] for i in list(range(0,numV))) <= budget))
		for i in list(range(0,budget)):
			cmodel.c[i]=1
		
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
		
	else:
		
		for i in list(range(0,numV)):
			cmodel.c[i] = budget/numV
		
		
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
	old_objval=pyo.value(cmodel.obj)
	print(old_objval)
	
	if budget == budgetlist[0]:
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
	
	#now lets throw the iterative solution into the full nlp for further optimization
	if budget==budgetlist[0]:
		cmodel.beta = pyo.Var(greekmodel.i, domain=pyo.NonNegativeReals)		
		cmodel.alpha = pyo.Var(indexset, domain=pyo.NonNegativeReals)
	
	for i in list(range(0,numV)):
		cmodel.beta[i]=greekmodel.beta[i].value
	
	for (i,j) in P.keys():
		cmodel.alpha[i,j]=greekmodel.alpha[i,j].value
	
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
	
	NLPSolver.options['max_iter']=200
	NLPSolver.options['max_cpu_time']=300
	
	if budget == budgetlist[0]:
		cmodel.con2 = pyo.ConstraintList()
		for i in list(range(0,numV)):
			cmodel.con2.add(expr=(-sum(cmodel.alpha[i,j] for (i,j) in adj[i]) + sum(cmodel.alpha[j,i] for (j,i) in radj[i]) + cmodel.beta[i] <= 1))
	else:
		cmodel.con2.activate()
		for i in list(range(0,numV)):
			cmodel.con2[i+1].set_value(expr=(-sum(cmodel.alpha[i,j] for (i,j) in adj[i]) + sum(cmodel.alpha[j,i] for (j,i) in radj[i]) + cmodel.beta[i] <= 1))
		
	cmodel.obj.set_value(expr=(-sum((1-P[i,j])*cmodel.alpha[i,j] for (i,j) in P.keys()) + sum(ProbFunc(Type[i],cmodel.c[i])*cmodel.beta[i] for i in greekmodel.i)))
	print("Solving full NLP")
	NLPSolver.solve(cmodel,tee=False)
	cmodel.con2.deactivate()
	
	NLPSolver.options['max_iter']=3000
	old_objval=pyo.value(cmodel.obj)
	
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
		
	LPSolver.solve(pimodel,tee=False)
	print(f'Pi values solved. NLP Objective value: {pyo.value(pimodel.obj)}')
	
	obj_vals[budget]=pyo.value(pimodel.obj)
	
	f=open(output_dir + f"/B={budget} discounts.csv",'w',newline='')
	writer=csv.writer(f)
	writer.writerow(("Robust Obj val:",pyo.value(pimodel.obj)))
	for i in list(range(0,numV)):
		if cmodel.c[i].value>0:
			writer.writerow((i,cmodel.c[i].value))
	f.close()
	