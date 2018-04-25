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

tuple Vec_2 {int i; int j; int k; int d;}
{Vec_2} Arcs_2={<i,j,k,d> | i in nodes , j in nodes, k in vehicle, d in depot}; 

tuple Vec_3 {int i; int p; int k; int d;}
{Vec_3} Arcs_3={<i,p,k,d> | i in nodes , p in product, k in vehicle, d in depot};

tuple Vec_4 {int i; int p;}
{Vec_4} Arcs_4={<i,p> | i in city , p in product};

tuple Vec_6 {int d; int p;}
{Vec_6} Arcs_6={<d,p> | d in depot , p in product};

int d[Arcs_4]=...;		// Demand of costumer i by product p 
float c[Arcs]=...;		// Cost to travel from city i to j
float v[product]=...;	// Capacidad requerida del producto P
int i[Arcs_6]=...;      // Amount of product p in depot d

//Model 
//Definition of Variable
dvar boolean x[Arcs_2]; 	// 1 if vehicle k, assigned to depot d, travel from i to j
dvar int+ z [Arcs_3];		// The amount of product p deliver to city i by vehicle k, assigned to depot d
dvar int+ u[Arcs_2];		// ranking in visit every city



//Objective Function
minimize sum (k in vehicle, i in nodes, j in nodes, d in depot)c[<i,j>]*x[<i,j,k,d>]; //(1)

subject to{

sum(d in depot,k in vehicle, j in city)x[<d,j,k,d>]<=K;  //(2)Son asignados K vehiculos

forall(d in depot) 
  sum(k in vehicle, j in city)x[<d,j,k,d>]<=K;  			//(3)No se pueden asignar a cada deposito más de K vehiculos

forall(k in vehicle)
  sum( j in nodes, d in depot)x[<d,j,k,d>]<=1;  //(4)*OJO Esta debe ser enfocada solo en los camiones que asigne :/ 

forall(i in depot, d in depot:i!=d)
  sum(k in vehicle, j in city)x[<i,j,k,d>]==0; //(5)Si el camion k es asignado al deposito d, su x debe tener ese valor

forall(j in depot, d in depot:j!=d)
  sum(k in vehicle, i in city)x[<i,j,k,d>]==0; //(6)Tienen que volver a su deposito original

 //forall(k in vehicle)
   //sum ( j in city, d in depot)x[<d,j,k,d>]<=1;  //Cada vehiculo debe salir del deposito asignado
 
 forall(i in depot, j in depot)
   sum(k in vehicle, d in depot)x[<j,i,k,d>]==0;   //(7)no hay rutas entre depositos
 
 forall(i in city)
   sum(k in vehicle, j in nodes, d in depot)x[<j,i,k,d>]>=1;  //(8) Cada ciudad de ser visitada
 
 forall(i in city, k in vehicle, d in depot)
   sum(j in nodes)x[<i,j,k,d>]-sum(j in nodes)x[<j,i,k,d>]==0; //(9)Camion que entra debe salir  
  
forall (i in city, k in vehicle, p in product, d in depot)
  z[<i,p,k,d>]<=10000*sum(j in nodes)x[<i,j,k,d>];          //(10)x puede existir solo si z existe
  
forall(i in city,p in product)
   sum(d in depot, k in vehicle)z[<i,p,k,d>]==d[<i,p>];      //(11)Se debe cumplir la demanda


forall(k in vehicle)
    sum(p in product, i in city, d in depot)z[<i,p,k,d>]*v[p]<=V;  //(12)No se puede superar el volumen
  
forall(p in product, d in depot)
	sum( k in vehicle,i in city)z[<i,p,k,d>]<=i[<d,p>];   //(13)Se debe respetar el inventario

//Variante GG-SEC
forall(i in nodes, j in nodes, k in vehicle, d in depot)
  u[<i,j,k,d>]<=V*x[<i,j,k,d>];
  
forall (i in city, k in vehicle, d in depot)
  sum(j in nodes:j!=i)u[<j,i,k,d>]-sum(j in city:j!=i)u[<i,j,k,d>]==sum(p in product)z[<i,p,k,d>];
  
  
 
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
forall(d in depot, j in nodes, k in vehicle)
  u[<j,k>]>=x[<d,j,k,d>];								//(14)

forall(i in city, k in vehicle, j in city, d in depot)	//(15)
  u[<i,k>]+1 <= u[<j,k>]+1000*(1-x[<i,j,k,d>]);
  
forall(i in nodes, k in vehicle)						//(16)
  1<=u[<i,k>]<=N-D;*/
 
}