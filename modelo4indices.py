	import time
	import numpy as np
	import cplex
	from cplex import Cplex
	from cplex.exceptions import CplexError
	import sys
	import networkx as nx
	import matplotlib.pyplot as plt

# DATA READING ------------------------------------------------------------------------------------
def Resolver(site,N,D,P,K,V):

	C = N - D
	M = 99999
	before = time.clock()
	T_exec = 3600

	print("N = {} | D = {} | P = {} | K = {} | V = {}\n".format(N,D,P,K,V))

	content = []
	with open(site+"_cij.txt") as f:
	    for i in f.readline().split():
	        content.append(int(i))
	c = np.resize(content,(N,N))
	#print(c)

	content1 = []
	with open(site+"_dcp.txt") as f:
	    for i in f.readline().split():
	        content1.append(int(i))
	#print(content1)
	dem = np.resize(content1,(N,P))
	#print(dem)

	content2 = []
	with open(site+"_sdp.txt") as f:
	    for i in f.readline().split():
	        content2.append(int(i))
	s = np.resize(content2,(D,P))
	#print(s)

	content3 = []
	with open(site+"_vp.txt") as f:
	    for i in f.readline().split():
	        content3.append(int(i))
	v = content3
	#print (v)

	# MATHEMATICAL MODEL ---------------------------------------------------------------------------

	Model = cplex.Cplex()

	x_vars = np.array([[[["x("+str(i)+","+str(j)+","+str(k)+","+str(d)+")" for d in range(0,D)] for k in range(0,K)] for j in range(0,N)] for i in range(0,N)])
	x_varnames = x_vars.flatten()
	x_vartypes = 'B'*len(x_varnames)
	x_varlb = [0.0]*len(x_varnames)
	x_varub = [1.0]*len(x_varnames)
	x_varobj = []
	for i in range(N):
		for j in range(N):
			for k in range(K):
				for d in range(D):
					x_varobj.append(float(c[i,j]))

	Model.variables.add(obj = x_varobj, lb = x_varlb, ub = x_varub, types = x_vartypes, names = x_varnames)

	z_vars = np.array([[[["z("+str(i)+","+str(p)+","+str(k)+","+str(d)+")" for d in range(0,D)] for k in range(0,K)] for p in range(0,P)] for i in range(0,N)])
	z_varnames = z_vars.flatten()
	z_vartypes = 'I'*N*P*K*D
	z_varlb = [0.0]*N*P*K*D
	z_varub = [cplex.infinity]*N*P*K*D
	z_varobj = [0.0]*N*P*K*D

	Model.variables.add(obj = z_varobj, lb = z_varlb, ub = z_varub, types = z_vartypes, names = z_varnames)

	u_vars = np.array([[[["u("+str(i)+","+str(j)+","+str(k)+","+str(d)+")" for d in range(0,D)] for k in range(0,K)] for j in range(0,N)] for i in range(0,N)])
	u_varnames = u_vars.flatten()
	u_vartypes = 'I'*N*N*K*D
	u_varlb = [0.0]*N*N*K*D
	u_varub = [cplex.infinity]*N*N*K*D
	u_varobj = [0.0]*N*N*K*D

	Model.variables.add(obj = u_varobj, lb = u_varlb, ub = u_varub, types = u_vartypes, names = u_varnames)

	Model.objective.set_sense(Model.objective.sense.minimize)

	# RESTRICCIONES ---------------------------------------------------------------------------------------------------

	row3 = []
	val3 = []
	for d in range(D):
		for k in range(K):
			for j in range(D,N):
				row3.append(x_vars[d,j,k,d])
				val3.append(1.0)

	Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row3, val = val3)], senses = 'L', rhs = [K])

	for d in range(D):
		row4 = []
		val4 = []
		for k in range(K):
			for j in range(D,N):
				row4.append(x_vars[d,j,k,d])
				val4.append(1.0)
		Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row4, val= val4)], senses = 'L', rhs = [K])

	for k in range(K):
		row5 = []
		val5 = []
		for j in range(N):
			for d in range(D):
				row5.append(x_vars[d,j,k,d])
				val5.append(1.0)
		Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row5, val= val5)], senses = 'L', rhs = [1.0])

	for i in range(D):
		for d in range(D):
			row6 = []
			val6 = []
			if(i!=d):
				for k in range(K):
					for j in range(D,N):
						row6.append(x_vars[i,j,k,d])
						val6.append(1.0)
			Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row6, val= val6)], senses = 'E', rhs = [0.0])

	for j in range(D):
		for d in range(D):
			row7 = []
			val7 = []
			if(j!=d):
				for k in range(K):
					for i in range(D,N):
						row7.append(x_vars[i,j,k,d])
						val7.append(1.0)
			Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row7, val= val7)], senses = 'E', rhs = [0.0])


	for i in range(D):
		for j in range(D):
			row8 = []
			val8 = []
			for k in range(K):
				for d in range(D):
					row8.append(x_vars[j,i,k,d])
					val8.append(1.0)
			Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row8, val= val8)], senses = 'E', rhs = [0.0])


	for i in range(D,N):
		row9 = []
		val9 = []
		for k in range(K):
			for j in range(N):
				for d in range(D):
					row9.append(x_vars[j,i,k,d])
					val9.append(1.0)		
		Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row9, val= val9)], senses = 'G', rhs = [1.0])
		
	for i in range(D,N):
		for k in range(K):
			for d in range(D):
				row10 = []
				val10 = []
				for j in range(N):
					if(j!=i):
						row10.append(x_vars[i,j,k,d])
						val10.append(1.0)
						row10.append(x_vars[j,i,k,d])
						val10.append(-1.0)
				Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row10, val= val10)], senses = 'E', rhs = [0.0])

	for i in range(D,N):
		for k in range(K):
			for p in range(P):
				for d in range(D):
					row11 = [z_vars[i,p,k,d]]
					val11 = [1.0]
					for j in range(N):
						row11.append(x_vars[j,i,k,d])
						val11.append(-M)
					Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row11, val= val11)], senses = 'L', rhs = [0.0])

	for i in range(D,N):
		for p in range(P):
			row12 = []
			val12 = []
			for d in range(D):
				for k in range(K):
					row12.append(z_vars[i,p,k,d])
					val12.append(1.0)		
			Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row12, val= val12)], senses = 'E', rhs = [float(dem[i,p])])
	
	for k in range(K):
		row13 = []
		val13 = []
		for p in range(P):
			for i in range(D,N):
				for d in range(D):
					row13.append(z_vars[i,p,k,d])
					val13.append(float(v[p]))
		Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row13, val= val13)], senses = 'L', rhs = [V])

	for p in range(P):
		for d in range(D):
			row14 = []
			val14 = []
			for k in range(K):
				for i in range(D,N):
					row14.append(z_vars[i,p,k,d])
					val14.append(1.0)
			Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row14, val= val14)], senses = 'L', rhs = [float(s[d,p])])

	for i in range(N):
		for j in range(N):
			for k in range(K):
				for d in range(D):
					row15=[]
					val15=[]
					row15.append(u_vars[i,j,k,d])
					val15.append(1.0)
					row15.append(x_vars[i,j,k,d])
					val15.append(-V)
					Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row15, val= val15)], senses = 'L', rhs = [0.0])

	for i in range(D,N):
		for k in range(K):
			for d in range(D):
				row16 = []
				val16 = []
				for j in range(N):
					if (j!=i):
						row16.append(u_vars[j,i,k,d])
						val16.append(1.0)
				for j in range(D,N):
					if (j!=i):
						row16.append(u_vars[i,j,k,d])
						val16.append(-1.0)
				for p in range(P):
					row16.append(z_vars[i,p,k,d])
					val16.append(-1.0)
				Model.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = row16, val= val16)], senses = 'E', rhs = [0.0])


	Model.parameters.timelimit.set(float(T_exec))
	Model.parameters.workmem.set(9000.0)

	Model.solve()

	# SHOW SOLUTION FUNCTION -------------------------------------------------------------------------------------------

	def show_solution():
		print("\nObjective Function Value = {}".format(Model.solution.get_objective_value()))
		for i in range(0,N):
			for j in range(0,N):
				for k in range(K):
					for d in range(D):
						if(Model.solution.get_values("x("+str(i)+","+str(j)+","+str(k)+","+str(d)+")")!=0.0):
							print("x("+str(i)+","+str(j)+","+str(k)+","+str(d)+")"+" = "+str(Model.solution.get_values("x("+str(i)+","+str(j)+","+str(k)+","+str(d)+")")))
		print("")

	def show_solutionZ():
		for i in range(0,N):
			for p in range(0,P):
				for k in range(K):
					for d in range(D):
						if(Model.solution.get_values("z("+str(i)+","+str(p)+","+str(k)+","+str(d)+")")!=0.0):
							print("z("+str(i)+","+str(p)+","+str(k)+","+str(d)+")"+" = "+str(Model.solution.get_values("z("+str(i)+","+str(p)+","+str(k)+","+str(d)+")")))
		print("")

	def show_solutionU():
		for i in range(0,N):
			for j in range(0,N):
				for k in range(K):
					for d in range(D):
						if(Model.solution.get_values("u("+str(i)+","+str(j)+","+str(k)+","+str(d)+")")!=0.0):
							print("u("+str(i)+","+str(j)+","+str(k)+","+str(d)+")"+" = "+str(Model.solution.get_values("u("+str(i)+","+str(j)+","+str(k)+","+str(d)+")")))
		print("")			

	# FUNCTIONS --------------------------------------------------------------------------------------------------------

	def create_X():
		X = []
		for i in range(0,N):
			for j in range(0,N):
				if(atsp.solution.get_values("x("+str(i)+","+str(j)+")")==1.0):
					X.append([i,j])

		#print("X = {}".format(X))
		return X

	def find_subtours():
		G = nx.DiGraph(X)
		S = list(nx.simple_cycles(G))
				
		#print("S = {}".format(S))
		return S

	def add_cuts():
		for subtour in S:
			Z = len(subtour)-1
			new_ind = []
			for i in subtour:
				for j in subtour:
					if j!=i:
						new_ind.append("x("+str(i)+","+str(j)+")")
			new_val = [1.0]*(len(new_ind))
			new_line = [cplex.SparsePair(ind = new_ind, val = new_val)]
			atsp.linear_constraints.add(lin_expr = new_line, senses = 'L', rhs = [Z])

	# SOLVING ------------------------------------------------------------------------------------

	show_solution()
	show_solutionZ()
	show_solutionU()

	#Graficar Solución -------------------------------------------------------------------
	nodes = list(range(D,N))
	depots= list(range(D))
	path=[]
	for k in range(K):
		for i in range(0,N):
			for j in range(0,N):
				for d in range(D):
					if(Model.solution.get_values("x("+str(i)+","+str(j)+","+str(k)+","+str(d)+")")!=0.0):
						path.append((i,j))
						
	G=nx.DiGraph()
	#G.add_nodes_from(nodes)
	#G.add_nodes_from(depots, color='blue')
	G.add_edges_from(path)

	nx.draw(G, with_labels = True)
	plt.savefig("Solution.png")
	plt.show() 

Resolver("Ejemplo",7,3,2,3,65)
#                  N,D,P,K ,V


