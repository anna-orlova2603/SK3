#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define M 500
#define N 1000

double u(double x, double y)
{
    return sqrt(4 + x*y);
}

double q(double x, double y)
{
	double t = x + y;
    return t >= 0 ? t : 0;
}

double k(double x, double y){
    return 4 + x + y;
}

/** phi_1(x) = u(4, y), R*/
double phi_1(double x, double y)
{
    return 2.0 * sqrt(1 + y);
}

/** phi_2(y) = u(0, y), L*/
double phi_2(double x, double y)
{
    return 2.0;
}

/** phi_3(x) = u(x, 3), T*/
double phi_3(double x, double y)
{
    return sqrt(4 + x*3);
}

/** phi_4(x) = u/dn(x, 3) + u, B*/
double phi_4(double x, double y)
{
    return -x/4 + 2; 
}

double F(double x, double y)
{
    return -(x + y)/(2*sqrt(x*y + 4)) + (y*y + x*x)*k(x, y)/(4*pow(x*y + 4, 1.5)) + q(x, y)*u(x, y);
}

int getColumns(int degree, int rank, int k)
{
    int numCol = 1;
    int columns[3];

    for(int i = 0; i < degree/2; i++)
    {
        numCol = numCol*2;
    }

    columns[0] = (N + 1) % numCol;
    columns[1] = (N + 1) / numCol;

    if ((rank + 1)%numCol)
    {
        columns[2] = (rank + 1) % numCol - 1;
    }
    else
    {
        columns[2] = numCol - 1;
    }
    
    if (k == 0)
        return columns[0];
    else if (k == 1)
        return columns[1];
    if (k == 2)
        return columns[2];
    return -1;
}

int getRows(int degree, int rank, int k)
{
    int numRow = 1, numCol = 1;
    int rows[3];

    for(int i = 0; i < degree/2; i++)
    {
        numCol = numCol * 2;
    }

    for(int i = 0; i < degree - degree/2; i++)
    {
        numRow = numRow * 2;
    }

    rows[0] = (M + 1) % numRow;
    rows[1] = (M + 1) / numRow;

    if ((rank + 1) % numCol)
    {
        rows[2] = (rank + 1) / numCol;
    }
    else
    {
        rows[2] = (rank + 1) / numCol - 1;
    }

    if (k == 0)
        return rows[0];
    else if (k == 1)
        return rows[1];
    if (k == 2)
        return rows[2];
    return -1;
}

double ro(int i, int j)
{
    if ((i == 0 && j == 0) || (i == M && j == 0) || (i == 0 && j == N) || (i == M && j == N)) 
        return 0.25;
    if (j == 0)
        return 0.5;
    return 1.0;
}

double iterPar(double u[], double v[], int i, int j, double h1, double h2)
{
    double res = 0.0;
    
    if ((i == 0 && j == 0) || (i == M && j == 0) || (i == 0 && j == N) || (i == M && j == N))
    {
        res = ro(i, j) * ro(i, j) * u[i * (N + 1) + j] * v[i * (N + 1) + j];
    } 
    else if(j == 0)
    {
        res = ro(i, j) * ro(i, j) * u[i * (N + 1) + j] * v[i * (N + 1) + j];
    }
    else if(i >= 1 && i <= M - 1 && j >= 1 && j <= N - 1)
    {
        res = ro(i, j) * ro(i, j) * u[i * (N + 1) + j] * v[i * (N + 1) + j];
    }

    return h1 * h2 * res;
}

double getB(int i, int j, double h1, double h2, int x_size, int y_size)
{
    double B = 0.0;
    double x = h1 * i, xl = h1 * (i - 1), xr = h1 * (i + 1);
    double y = h2 * j, yt = h2 * (j + 1), yb = h2 * (j - 1);
    
    if(i >= 2 && i <= M - 2 && j >= 1 && j <= N - 2)   //internal
    {
        B = F(x, y);
    }
    else if(i == 1 && j >= 1 && j <= N - 2)            //left side a1b1 -> a1b2 
    {
        B = F(x, y) + 2.0/(h1) * phi_2(xl, y);
    } 
    else if(i == M - 1 && j >= 1 && j <= N - 2)        //right side a2b1 -> a2b2
    {
        B = F(x, y) + 2.0/(h1) * phi_1(xr, y);
    }
    else if(i >= 2 && i <= M - 2 && j == N - 1)        //top side a1b2 -> a2b2 !!
    {
        B = F(x, y) + 2.0/(h2) * phi_3(xr, yt);//
    }
    else if(i >= 2 && i <= M - 2 && j == 0)            //bottom side a1b1 -> a2b1
    {
        B = F(x, y) + (2.0/h2) * phi_4(x, y);
    }
    else if(i == 1 && j == 0)                          //point 1, 0
    {
        B = F(x, y) + 2.0/h2 * phi_4(xl, y) + 2.0/h1 * phi_2(xl, y);
    }
    else if(i == 1 && j == N - 1)                          //point 1, N-1
    {
        B = F(x, y) + 2.0/h2 * phi_3(xl, yt) + 2.0/h1 * phi_2(xl, yt);
    }
    else if(i == M - 1 && j == 0)                          //point M-1, 0
    {
        B = F(x, y) + 2.0/h2 * phi_4(xr, y) + 2.0/h1 * phi_1(xr, y);
    }
    else if(i == M - 1 && j == N - 1)                          //point M-1, N-1
    {
        B = F(x, y) + 2.0/h1 * phi_1(xr, yt) + 2.0/h2 * phi_3(xr, yt);
    }

    return B; 
    
}

double getAw(double w[], int i, int j, double h1, double h2, int x_size, int y_size)
{
    double aw;
    double x = h1 * i, xl = h1 * (i - 1), xr = h1 * (i + 1);
    double y = h2 * j, yt = h2 * (j + 1), yb = h2 * (j - 1);


    if(i >= 2 && i <= M - 2 && j >= 1 && j <= N - 2)   //internal
    {
        aw = 1.0/(h1*h1) * (k(x + 0.5 * h1, y) * (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) 
        - k(x - 0.5 * h1, y) * (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]))
        + 1.0/(h2*h2) * (k(x, y + 0.5 * h2) * (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j]) 
        + k(x, y - 0.5 * h2) * (w[i * (N + 1) + j] - w[i * (N + 1) + j - 1])) 
        + q(x, y) * w[i * (N + 1) + j];
    }
    else if(i == 1 && j >= 1 && j <= N - 2)            //left side a1b1 -> a1b2 
    {
        aw = -1.0/(h1*h1) * (k(xr - 0.5 * h1, y) * (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) 
        - 1.0/(h1*h1) * k(x - 0.5 * h1, y) * w[i * (N + 1) + j])
        - 1.0/(h2*h2) * (k(x, y + 0.5 * h2) * (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j]) 
        - k(x, y - 0.5 * h2) * (w[i * (N + 1) + j] - w[i * (N + 1) + j - 1]))
        + q(x, y) * w[i * (N + 1) + j];
    } 
    else if(i == M - 1 && j >= 1 && j <= N - 2)        //right side a2b1 -> a2b2
    {
        aw = 1.0/(h1*h1) * (k(xl + 0.5 * h1, y)) * (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]) 
        + 1.0/(h1*h1) * (k(x + 0.5 * h1, y)) * w[i * (N + 1) + j]  
        - 1.0/(h2*h2) * (k(x, y + 0.5 * h2) * (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j]) 
         - k(x, y - 0.5 * h2) * (w[i * (N + 1) + j] - w[i * (N + 1) + j - 1]))
         + q(x, y) * w[i * (N + 1) + j];
    }
    else if(i >= 2 && i <= M - 2 && j == N - 1)       //top side a1b2 -> a2b2 !!
    {
        aw = -1.0/(h1*h1) * (k(x + 0.5 * h1, y) * (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) 
        - k(x - 0.5 * h1, y) * (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j])) 
        + 1.0/(h2*h2) * (k(x, yb + 0.5 * h2) * (w[i * (N + 1) + j] 
        - w[i * (N + 1) + j - 1]) + k(x, y + 0.5 * h2) * w[i * (N + 1) + j])
        + q(x, y) * w[i * (N + 1) + j];
    }
    else if(i == 1 && j == N - 1)                     //point a1b2 !! -1
    {
        aw = -1.0/(h1*h1) * k(xr - 0.5 * h1, y) * (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) 
        + 1.0/(h1*h1) * k(x - 0.5 * h1, y) * w[i * (N + 1) + j] 
        + 1.0/(h2*h2) * (k(x, yb + 0.5 * h2) * (w[i * (N + 1) + j] - w[i * (N + 1) + j - 1]) 
        + k(x, y + 0.5 * h2) * w[i * (N + 1) + j])
        + q(x, y) * w[i * (N + 1) + j] ;
    }
    else if(i == M - 1 && j == N - 1)                 //point a2b2 !!
    {
        aw = 1.0/(h1*h1) * (k(xl + 0.5 * h1, y) * (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]) 
        + k(x + 0.5 * h1, y) * w[i * (N + 1) + j]) 
        + 1.0/(h2*h2) * (k(x, yb + 0.5 * h2) * (w[i * (N + 1) + j] - w[i * (N + 1) + j - 1]) 
        + k(x, y + 0.5* h2) * w[i * (N + 1) + j])
        + q(x, y) * w[i * (N + 1) + j];
    }
    else if(i >= 2 && i <= M - 1 && j == 0)            //bottom side a1b1 -> a2b1
    {
        aw = -2.0/(h2*h2) * k(x, j + 0.5 * h2) * (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j]) 
        - 1.0/(h1*h1) * (k(x + 0.5 * h1, y) * (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) 
        - k(x - 0.5 * h1, y) * (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]))
         + (q(x, y) + 2/h2) * w[i * (N + 1) + j] ;
    }
    else if(i == 1 && j == 0)                          //point a1b1
    {
        aw = -1.0/(h1*h1) * k(xr - 0.5 * h1, y) * (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) 
        + 1.0/(h1*h1) * k(x - 0.5 * h1, y) * w[i * (N + 1) + j] //(y_size + 2)
        - 2.0/(h2*h2) * k(x, yt - 0.5 * h2) * (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j])
        + (q(x, y) + 2/h2) * w[i * (N + 1) + j];
    }
    else if(i == M - 1 && j == 0)                          //point a2b1
    {
        aw = 1.0/(h1*h1) * (k(xl + 0.5 * h1, y) * (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]) 
        + k(x + 0.5 * h1, y) * w[i * (N + 1) + j]) 
        - 2.0/(h2*h2) * k(x, yt - 0.5 * h2) * (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j])
        + (q(x, y) + (2/h2)) * w[i * (N + 1) + j];
    }
    else if (i == 0 || i == M || j == N)
    {
        aw = w[i * (N + 1) + j];
    }

    return aw;

}

int main(int argc, char *argv[]){

    int iter = 0,
        rank,
        size,
        root = 0, i, j;

    double epsilon = 1e-6,
           eps_local,
           eps_global = 0.0,
           tau, 
           tau_local[2] = {0.0, 0.0}, 
           tau_global[2] = {0.0, 0.0},
           maxTime,
           time,
           start;

    double A1 = 0.0,
           A2 = 4.0,
           B1 = 0.0,
           B2 = 3.0,
           h1 = (A2 - A1) / M,
           h2 = (B2 - B1) / N;

    double w[(M + 1) * (N + 1)],
           r2[(M + 1) * (N + 1)], 
           Ar2[(M + 1) * (N + 1)], 
           wk[(M + 1) * (N + 1)],
           error[(M + 1) * (N + 1)],
           r[(M + 1) * (N + 1)],
           Ar[(M + 1) * (N + 1)],
           fu[(M + 1) * (N + 1)];

    for(int i = 0; i <= M; i++)
        {
            for (int j = 0; j <= N; j++)
            {
                fu[i * (N + 1) + j] = u(i * A2 / M, j * B2 / N);
                w[i * (N + 1) + j] = 0.0;
                wk[i * (N + 1) + j] = 0.0;
                r[i * (N + 1) + j] = 0.0;
                r2[i * (N + 1) + j] = 0.0;
                Ar[i * (N + 1) + j] = 0.0;
                Ar2[i * (N + 1) + j] = 0.0;
                error[i * (N + 1) + j] = 0.0;
            }
        }
  
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    start = MPI_Wtime();

    int degree = 0;
    int a = size;

    int x_size[2];
    int y_size[2];
    int coord[2];
    int procRow[2];
    int procCol[2];

    while (a != 1) 
    {
        a /= 2;
        degree++;
    }

    procRow[0] = getRows(degree, rank, 0); procRow[1] = getRows(degree, rank, 1); coord[0] = getRows(degree, rank, 2);
    procCol[0] = getColumns(degree, rank, 0); procCol[1] = getColumns(degree, rank, 1); coord[1] = getColumns(degree, rank, 2);

    if(((coord[0] + 1) <= procRow[0]) && ((coord[1] + 1) <= procCol[0]))
    {
        x_size[0] = coord[0] * (procRow[1] + 1);
        x_size[1] = (coord[0] + 1) * (procRow[1] + 1);
        y_size[0] = coord[1] * (procCol[1] + 1);
        y_size[1] = (coord[1] + 1) * (procCol[1] + 1);
    }
    else if (((coord[0] + 1) > procRow[0]) && ((coord[1] + 1 <= procCol[0])))
    {
        x_size[0] = procRow[0] * (procRow[1] + 1) + (coord[0] - procRow[0]) * procRow[1];//procRow[0] + procRow[1] * coord[0]; //
        x_size[1] = procRow[0] * (procRow[1] + 1) + (coord[0] - procRow[0] + 1) * procRow[1];
        y_size[0] = coord[1] * (procCol[1] + 1);
        y_size[1] = (coord[1] + 1) * (procCol[1] + 1);
    }
    else if (((coord[0] + 1) <= procRow[0]) && ((coord[1] + 1) > procCol[0]))
    {
        x_size[0] = coord[0] * (procRow[1] + 1);
        x_size[1] = (coord[0] + 1) * (procRow[1] + 1);
        y_size[0] = procCol[0] * (procCol[1] + 1) + (coord[1] - procCol[0]) * procCol[1];
        y_size[1] = procCol[0] * (procCol[1] + 1) + (coord[1] - procCol[0] + 1) * procCol[1];
    }
    else
    {
        x_size[0] = procRow[0] * (procRow[1] + 1) + (coord[0] - procRow[0]) * procRow[1]; //
        x_size[1] = procRow[0] * (procRow[1] + 1) + (coord[0] - procRow[0] + 1) * procRow[1];
        y_size[0] = procCol[0] * (procCol[1] + 1) + (coord[1] - procCol[0]) * procCol[1];
        y_size[1] = procCol[0] * (procCol[1] + 1) + (coord[1] - procCol[0] + 1) * procCol[1];
    }

    do{
        iter++;

        //r = A*w - B
        for(int i = x_size[0]; i < x_size[1]; i++)
        {
            for (int j = y_size[0]; j < y_size[1]; j++)
            {
                r2[i * (N + 1) + j] = getAw(w, i, j, h1, h2, M, N) - getB(i, j, h1, h2, M, N);// 
            }
        }

        MPI_Allreduce(r2, r, (M + 1) * (N + 1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //A*r
        for(int i = x_size[0]; i < x_size[1]; i++)
        {
            for (int j = y_size[0]; j < y_size[1]; j++)
            {
                Ar2[i * (N + 1) + j] = getAw(r, i, j, h1, h2, M, N);
            }
        }

        MPI_Allreduce(Ar2, Ar, (M + 1) * (N + 1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        tau_local[0] = 0.0;
        tau_local[1] = 0.0;

        //r_k+1
        for(int i = x_size[0]; i < x_size[1]; i++)
        {
            for (int j = y_size[0]; j < y_size[1]; j++)
            {
                tau_local[0] += iterPar(Ar, r, i, j, h1, h2);
                tau_local[1] += iterPar(Ar, Ar, i, j, h1, h2);
            }
        }

        MPI_Allreduce(tau_local, tau_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        tau = tau_global[0] / tau_global[1];
        
        for(i = 0; i <= M; i++)
        {
            for(j = 0; j <= N; j++)
            {   
                wk[i * (N + 1) + j] = (w[i * (N + 1) + j] - tau * r[i * (N + 1) + j]);//
                if (i == 0 || i == M || j == N)
                {
                    error[i * (N + 1) + j] = 0.0;
                }
                else
                {
                    error[i * (N + 1) + j] = (wk[i * (N + 1) + j] - w[i * (N + 1) + j]);
                }
                w[i * (N + 1) + j] = wk[i * (N + 1) + j];
                r2[i * (N + 1) + j] = 0;
                Ar2[i * (N + 1) + j] = 0;
            }
        }

        eps_local = 0.0;

        for(int i = x_size[0]; i < x_size[1]; i++)
        {
            for (int j = y_size[0]; j < y_size[1]; j++)
            {
                eps_local += sqrt(iterPar(error, error, i, j, h1, h2));
            }
        }
        MPI_Allreduce(&eps_local, &eps_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    } while(eps_global > epsilon);
    
    
    time = MPI_Wtime() - start;
    MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    MPI_Finalize();

    if(rank == root)
    {
        printf("Iterations: %d, time: %f seconds\n", iter, maxTime);
        printf("Number of processes: %d, epsilon: %f, M: %d, N: %d\n", size, epsilon, M, N);
    }
    
    return 0;
}
