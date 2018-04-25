/*********************************************
 * OPL 12.7.1.0 Model
 * Author: Gonzalo
 * Creation Date: 17-10-2017 at 16:21:52
 *********************************************/
int N=...;		//Number of nodes
int K=...;		//Number of vehicle
float V=...;	//Capacity of vehicle
int P=...;		//Number of product 
int D=...;		//Number of depot

range nodes=1..N;
range city=(D+1)..N;
range depot=1..D;
range vehicle=1..K;
range product=1..P;

tuple Vec {int i; int j;}
{Vec} Arcs={<i,j> | i in nodes , j in nodes};

tuple Vec_2 {int i; int j; int k;}
{Vec_2} Arcs_2={<i,j,k> | i in nodes , j in nodes, k in vehicle}; 

tuple Vec_3 {int i; int p; int k;}
{Vec_3} Arcs_3={<i,p,k> | i in nodes , p in product, k in vehicle};

tuple Vec_4 {int i; int p;}
{Vec_4} Arcs_4={<i,p> | i in city , p in product};

tuple Vec_6 {int d; int p;}
{Vec_6} Arcs_6={<d,p> | d in depot , p in product};

tuple Vec_7 {int i; int k;}
{Vec_7} Arcs_7={<i,k> | i in nodes , k in vehicle};

int d[Arcs_4]=...;		// Demand of costumer i by product p 
float c[Arcs]=...;		// Cost to travel from city i to j
float v[product]=...;	// Capacidad requerida del producto P
int i[Arcs_6]=...;      // Amount of product p in depot d
int y[depot]=...;		// vehicle by depot

//Model 
//Definition of Variable
dvar boolean x[Arcs_2]; 	// 1 if vehicle k, travel from i to j
dvar int+ z [Arcs_3];		// The amount of product p deliver to city i by vehicle k, assigned to depot d
dvar int+ u[Arcs_2];		// ranking in visit every city




//Objective Function
minimize sum (k in vehicle, i in nodes, j in nodes)c[<i,j>]*x[<i,j,k>]; //(1)

subject to{

//Ruteo

forall(d in depot)
  sum(k in vehicle, j in city)x[<d,j,k>]==y[d];  			//(3)No se pueden asignar a cada deposito más de K vehiculos

forall(k in vehicle)
  sum( j in nodes, d in depot)x[<d,j,k>]<=1;  //(4)Cada camion es asignado solo una vez

//forall(d in depot,k in vehicle, j in city)
//	x[<j,d,k>] == x[<d,j,k>];  //Fixed destination 	 

 forall(i in depot, j in depot)
   sum(k in vehicle)x[<j,i,k>]==0;   //(7)no hay rutas entre depositos
 
 forall(i in city)
   sum(k in vehicle, j in nodes)x[<j,i,k>]>=1;  //(8) Cada ciudad de ser visitada
  
 forall(i in city, k in vehicle)
   sum(j in nodes)x[<i,j,k>]-sum(j in nodes)x[<j,i,k>]==0; //(9)Camion que entra debe salir 
   

//Linked constrain  
forall (i in city, k in vehicle, p in product)
  z[<i,p,k>]<=10000*sum(j in nodes)x[<i,j,k>];          //(10)x puede existir solo si z existe

//Materiales
forall(i in city,p in product)
   sum(k in vehicle)z[<i,p,k>]==d[<i,p>];      //(11)Se debe cumplir la demanda

forall(k in vehicle)
    sum(p in product, i in city)z[<i,p,k>]*v[p]<=V;  //(12)No se puede superar el volumen
 
forall(p in product, d in depot)
	sum( k in vehicle)z[<d,p,k>]<=i[<d,p>];   //(13)Se debe respetar el inventario


//Variante GG-SEC
forall(i in nodes, j in nodes, k in vehicle)
  u[<i,j,k>]<=V*x[<i,j,k>];
  
forall (i in city, k in vehicle)
  sum(j in nodes:j!=i)u[<j,i,k>]-sum(j in city:j!=i)u[<i,j,k>]==sum(p in product)z[<i,p,k>];

 











/*
//Profe-SEC
forall(i in nodes, j in nodes, k in vehicle, d in depot:i!=j)
	u[<i,j,k,d>]<=D*x[<i,j,k,d>];

forall(i in city, k in vehicle, d in depot)
  u[<d,i,k,d>]==d*x[<d,i,k,d>];
  
forall(i in city, k in vehicle, d in depot)
  u[<i,d,k,d>]==d*x[<i,d,k,d>];
  
forall(i in city, k in vehicle, d in depot)
  sum(j in nodes:i!=j)u[<j,i,k,d>]-sum(j in nodes:i!=j)u[<i,j,k,d>]==0;
*/


//GG-SEC
/*
forall(i in nodes, j in nodes, k in vehicle, d in depot:i!=j)
  u[<i,j,k,d>]<=V*x[<i,j,k,d>];
  
//forall(i in nodes, j in nodes, k in vehicle, d in depot:i!=j)
  //0<=u[<i,j,k,d>];  
 
forall (i in city, k in vehicle, d in depot)
  sum(j in nodes:j!=i)u[<j,i,k,d>]-sum(j in city:j!=i)u[<i,j,k,d>]==1;
 */


/*

//MTZ-SEC  
forall(j in nodes, k in vehicle, d in depot)
  u[<j,k>]>=x[<d,j,k>];								//(14)

forall(i in city, k in vehicle, j in city)	//(15)
  u[<i,k>]+1 <= u[<j,k>]+1000*(1-x[<i,j,k>]);
  
forall(i in nodes, k in vehicle)						//(16)
  1<=u[<i,k>]<=N-D;
 */
}