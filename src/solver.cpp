#include "solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h> 
#include <math.h>   
#include <algorithm>

void Solver::Init(unsigned N, float dt, float diff, float visc)
{
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;
	this->N = N;
    this->generators = true;
}

/*
float *u_prev, *v_prev, *dens_prev;
*/
/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/
void Solver::FreeData(void)
{
	//TODO: Libera los buffers de memoria.
	int flatDim = (N + 2)*(N + 2);
	if (u_prev != nullptr) u_prev = (float*)malloc(flatDim * sizeof(float));
	if (v_prev != nullptr) v_prev = (float*)malloc(flatDim * sizeof(float));
	if (dens_prev != nullptr) dens_prev = (float*)malloc(flatDim * sizeof(float));
	if (u != nullptr) u = (float*)malloc(flatDim * sizeof(float));
	if (v != nullptr) v = (float*)malloc(flatDim * sizeof(float));
	if (dens != nullptr) dens = (float*)malloc(flatDim * sizeof(float));
}

void Solver::ClearData(void)
{
	//TODO: Borra todo el contenido de los buffers
	int flatDim = (N + 2)*(N + 2);
	memset(u_prev, 0.0f, flatDim * sizeof(float));
	memset(v_prev, 0.0f, flatDim * sizeof(float));
	memset(dens_prev, 0.0f, flatDim * sizeof(float));
	memset(u, 0.0f, flatDim * sizeof(float));
	memset(v, 0.0f, flatDim * sizeof(float));
	memset(dens, 0.0f, flatDim * sizeof(float));
}

bool Solver::AllocateData(void)
{
	//TODO:
	//Reservamos memoria, en caso de fallo devlvemos false.
	//Antes de devolver true, hay que limpiar la memoria reservada con un ClearData().
	int flatDim = (N + 2)*(N + 2);
	u_prev = (float*)malloc(flatDim * sizeof(float));
	v_prev = (float*)malloc(flatDim * sizeof(float));
	dens_prev = (float*)malloc(flatDim * sizeof(float));
	u = (float*)malloc(flatDim * sizeof(float));
	v = (float*)malloc(flatDim * sizeof(float));
	dens = (float*)malloc(flatDim * sizeof(float));

	ClearData();

	return u_prev != nullptr && v_prev != nullptr && dens_prev != nullptr &&
		u != nullptr && v != nullptr && dens != nullptr;
}

void Solver::ClearPrevData() 
{
	//TODO: Borra el contenido de los buffers _prev
	int flatDim = (N + 2)*(N + 2);
	memset(u_prev, 0.0f, flatDim * sizeof(float));
	memset(v_prev, 0.0f, flatDim * sizeof(float));
	memset(dens_prev, 0.0f, flatDim * sizeof(float));
}

void Solver::AddDensity(unsigned x, unsigned y, float source)
{
	//TODO: Añade el valor de source al array de densidades. Sería interesante usar la macro: XY_TO_ARRAY
	dens_prev[XY_TO_ARRAY(x, y)] += source;
}

void Solver::AddVelocity(unsigned x, unsigned y, float forceX, float forceY)
{
	//TODO: Añade el valor de fuerza a sus respectivos arrays. Sería interesante usar la macro: XY_TO_ARRAY
	u_prev[XY_TO_ARRAY(x, y)] += forceX;
	v_prev[XY_TO_ARRAY(x, y)] += forceY;
}

void Solver::Solve()
{
	VelStep();
	DensStep();
}

void Solver::DensStep()
{
    if (generators) {
        dens[XY_TO_ARRAY(1, N / 2)] = 100;
        dens[XY_TO_ARRAY(N, 4)] = -10;
    }

	AddSource(dens, dens_prev);			//Adding input density (dens_prev) to final density (dens).
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Diffuse(0, dens, dens_prev);		//Writing result in dens because we made the swap before. bi = dens_prev. The initial trash in dens matrix, doesnt matter, because it converges anyways.
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Advect(0, dens, dens_prev, u, v);	//Advect phase, result in dens.
}

void Solver::VelStep()
{
    if (generators) {
        v[XY_TO_ARRAY(N, 3)] = 0.5;
        v[XY_TO_ARRAY(1, 2 * N / 3)] = -1.5;
        u[XY_TO_ARRAY(N, 3)] = -1;
        u[XY_TO_ARRAY(1, 2 * N / 3)] = 1;
    }

	AddSource(u, u_prev);
	AddSource(v, v_prev);
	SWAP (u_prev, u)			
	SWAP (v_prev, v)
	Diffuse(1, u, u_prev);  
	Diffuse(2, v, v_prev); 
	Project(u, v, u_prev, v_prev);		//Mass conserving.
	SWAP (u_prev, u)			
	SWAP (v_prev, v)
	Advect(1, u, u_prev, u_prev, v_prev);
	Advect(2, v, v_prev, u_prev, v_prev);
	Project(u, v, u_prev, v_prev);		//Mass conserving.
}

void Solver::AddSource(float * base, float * source)
{
	//TODO: Teniendo en cuenta dt (Delta Time), incrementar el array base con nuestro source. 
	//      Esto sirve tanto para añadir las nuevas densidades como las nuevas fuerzas.
	int i, j;
	FOR_EACH_CELL
		base[XY_TO_ARRAY(i, j)] += source[XY_TO_ARRAY(i, j)] * dt;
	END_FOR
}


void Solver::SetBounds(int b, float * x)
{
/*TODO:
Input b: 0, 1 or 2.
	0: borders = same value than the inner value.
	1: x axis borders inverted, y axis equal.
	2: y axis borders inverted, x axis equal.
	Corner values allways are mean value between associated edges.
*/
	float xSign = b == 1 ? -1 : 1;
	float ySign = b == 2 ? -1 : 1;

	for (int i = 1; i <= N; ++i) {
		x[XY_TO_ARRAY(i, 0)] = xSign * x[XY_TO_ARRAY(i, 1)];
		x[XY_TO_ARRAY(i, N + 1)] = xSign * x[XY_TO_ARRAY(i, N)];
		x[XY_TO_ARRAY(0, i)] = ySign * x[XY_TO_ARRAY(1, i)];
		x[XY_TO_ARRAY(N + 1, i)] = ySign * x[XY_TO_ARRAY(N, i)];
	}

	x[XY_TO_ARRAY(0, 0)] = (x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]) / 2;
	x[XY_TO_ARRAY(N + 1, 0)] = (x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]) / 2;
	x[XY_TO_ARRAY(0, N + 1)] = (x[XY_TO_ARRAY(1, N + 1)] + x[XY_TO_ARRAY(0, N)]) / 2;
	x[XY_TO_ARRAY(N + 1, N + 1)] = (x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]) / 2;

}

/*
https://www.youtube.com/watch?v=62_RUX_hrT4
https://es.wikipedia.org/wiki/M%C3%A9todo_de_Gauss-Seidel <- Solución de valores independientes.
Despreciando posibles valores de x no contiguos, se simplifica mucho. Mirar diapositivas y la solución de Gauss Seidel de términos independientes.
Gauss Seidel -> Matrix x and x0
*/
void Solver::LinSolve(int b, float * x, float * x0, float aij, float aii) {
    //TODO: Se recomienda usar FOR_EACH_CELL, END_FOR y XY_TO_ARRAY.
    int i, j;
	for (int k = 0; k < SOLVE_ITERATIONS; ++k) {
		FOR_EACH_CELL
			x[XY_TO_ARRAY(i, j)] = (-aij * (-x[XY_TO_ARRAY(i, j - 1)] - x[XY_TO_ARRAY(i - 1, j)] - x[XY_TO_ARRAY(i + 1, j)] - x[XY_TO_ARRAY(i, j + 1)]) + x0[XY_TO_ARRAY(i, j)]) / aii;
		END_FOR
		SetBounds(b, x);
	}
}

/*
Nuestra función de difusión solo debe resolver el sistema de ecuaciones simplificado a las celdas contiguas de la casilla que queremos resolver,
por lo que solo con la entrada de dos valores, debemos poder obtener el resultado.
*/
void Solver::Diffuse(int b, float * x, float * x0) {
//TODO: Solo necesitaremos pasar dos parámetros a nuestro resolutor de sistemas de ecuaciones de Gauss Seidel. Calculamos dichos valores y llamamos a la resolución del sistema.
    float aij = diff * dt * N * N;
    LinSolve(b, x, x0, aij, 1 + 4 * aij);
}

/*
d is overwrited with the initial d0 data and affected by the u & v vectorfield.
Hay que tener en cuenta que el centro de las casillas representa la posición entera dentro de la casilla, por lo que los bordes estan
en las posiciones x,5.
*/
void Solver::Advect(int b, float * d, float * d0, float * u, float * v) {
//TODO: Se aplica el campo vectorial realizando una interploación lineal entre las 4 casillas más cercanas donde caiga el nuevo valor.
	using namespace std;
	float dir = -dt * N, 
        currU, currV,
        deltaX, deltaY, 
        oX_f, oY_f;
	int i, j, 
        cellDest, 
        cell1, cell2, cell3, cell4,
        oX, oY;

	FOR_EACH_CELL
		cellDest = XY_TO_ARRAY(i, j);
        currU = u[cellDest];
        currV = v[cellDest];

        oX_f = i + (currU * dir);
        oY_f = j + (currV * dir);

        oX = (int)oX_f;
        oY = (int)oY_f;
        if (oX > 0 && oX < N && oY > 0 && oY < N) {
            cell1 = XY_TO_ARRAY(oX, oY);
            cell2 = XY_TO_ARRAY(oX + 1, oY);
            cell3 = XY_TO_ARRAY(oX, oY + 1);
            cell4 = XY_TO_ARRAY(oX + 1, oY + 1);

            deltaX = fabs(oX_f - oX);
            deltaY = fabs(oY_f - oY);

            d[cellDest] = d0[cell1] * (1 - deltaX) * (1 - deltaY) +
                          d0[cell2] * deltaX * (1 - deltaY) +
                          d0[cell3] * (1 - deltaX) * deltaY +
                          d0[cell4] * deltaX * deltaY;
        } else {
            d[cellDest] = 0;
        }
	END_FOR

	SetBounds(b, d);
}

/*
Se encarga de estabilizar el fluido y hacerlo conservativo de masa. Se usa solo en las matrices de velocidades.
No necesaria implementación por su complejidad.
*/
void Solver::Project(float * u, float * v, float * p, float * div)
{
	int i, j;

	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
		p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
	SetBounds(0, div);
	SetBounds(0, p);

	LinSolve(0, p, div, 1, 4);

	//Aproximamos: Laplaciano de q a su gradiente.
	FOR_EACH_CELL
		u[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
		v[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
	END_FOR
	SetBounds(1, u);
	SetBounds(2, v);
}