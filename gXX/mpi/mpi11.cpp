#include <iostream>
#include <vector>
#include <cstdlib>
#include <omp.h>
#include <cmath>
#include <mpi.h>
#include <algorithm>

#define N_SPECIES 9
using namespace std;

unsigned int seed;

// Initial seed
void init_r4uni(int input_seed)
{
    seed = input_seed + 987654321;
}

// Randomizer
float r4_uni()
{
    int seed_in = seed;

    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);

    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}

// Genaration of initial grid
void gen_initial_grid(int*** grid, long long N, float density, int input_seed, int id, int p)
{
    int start = id * floor(N / p) + min(id, ((int)N)%p);
    int end = (id+1)* floor(N / p) +  min(id+1, ((int)N)%p) -1;
    //cout << "aaa" << id << " " << start << " " << end << endl;

    init_r4uni(input_seed);
    for (int x = 1; x < N+1; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                if(r4_uni() < density) {
                    int species = (int)(r4_uni() * N_SPECIES) + 1;
                    if( x>=start+1 && x<=end+1) {
                        grid[x-start][y][z] = species;
                    }
                }                
            }
        }
    }
}

int** gen_plaquette(long long N) {
    int** grid = (int **) malloc(N * sizeof(int *));
    if(grid == NULL) {
        cout <<"Failed to allocate matrix\n";
        exit(1);
    }
    for(int y = 0; y < N; y++){
        grid[y] = (int *) malloc(N * sizeof(int));
        if(grid[y] == NULL) {
            cout <<"Failed to allocate matrix\n";
            exit(1);
        }
        for(int z = 0; z < N; z++){
            grid[y][z] = 0;
        }
    }
    return grid;
}

int*** gen_grid(long long N, int id, int p) {
    int*** grid;
    int size;
    double div = (float) N/p;
    // TODO NPROC TEM DE SER MENOR QUE N
    if (id >= N%p) size = floor(div);
    else size = ceil(div);
    //cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"<< id<<" " <<size<<" " <<(id == p-1) << (N%p != 0) << " "<< floor(N/p) << endl;
    if(size == 0) {
        cout <<"Failed to allocate matrix\n";
        exit(1);
    }
    grid = (int ***) malloc((size+2) * sizeof(int **));

    if(grid == NULL) {
        cout <<"Failed to allocate matrix\n";
        exit(1);
    }
    for(int x = 0; x < size+2; x++) {
        grid[x] = gen_plaquette(N);
        if(grid[x] == NULL) {
        cout <<"Failed to allocate matrix\n";
        exit(1);
    }
    }
    //cout << sizeof(grid)/sizeof(int **) << endl;
    return grid;
}



// Prints the grid
void print_grid(int*** grid, long long N) {
    for (int x = 0; x < N; x++) {
        cout << "Layer " << x << ":" << endl;
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                cout << grid[x][y][z] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}
// Prints the grid
void print_grid(int** grid, long long N) {
    for (int y = 0; y < N; y++) {
        for (int z = 0; z < N; z++) {
            cout << grid[y][z] << " ";
        }
        cout << endl;
    }
}
// Prints the results
void print_results(int** maximum){
    for(int i = 0; i < N_SPECIES; i++){
        cout << i+1 << " " << maximum[i][0] << " " << maximum[i][1] << endl;
    }
}

// Compare and modify function
void compare_and_modify(int* max, int** maximum, int epoch)
{
    for (int i = 0; i < N_SPECIES; i++) {
        if (max[i] > maximum[i][0]) {
            maximum[i][0] = max[i];
            maximum[i][1] = epoch;
        }
    }
}


/*
Duvidas:
- tem de se mudar todos os x,y,z para long long???



*/

// Fills vectors
void fill_vector(int* dimentions, int dimention, long long N) {
    if (dimention == 0) dimentions[0] = N-1;
    else dimentions[0] = dimention-1;

    dimentions[1] = dimention;

    if (dimention == N-1) dimentions[2] = 0;
    else dimentions[2] = dimention+1;
}

int alive_neighbour(int* max) {
    int m = 0;
    int index = 0;
    for(int i = 0; i < N_SPECIES; i++) {
        if (max[i] > m){
            m = max[i];
            index = i;
        }
    }
    return index+1;
}

int alive_or_dead(int*** grid, int* max, int sum, int x, int y, int z) {
    if(grid[x][y][z] != 0){
        sum--;
        max[grid[x][y][z]-1]--;
    } 
    if(grid[x][y][z] == 0 && sum >= 7 && sum <= 10){
        return alive_neighbour(max);
    } 
    else if(sum <= 4 || sum > 13) return 0;
    return grid[x][y][z];
}

// Visits node
int visit_node(int*** grid, long long N, int x, int y, int z) {
    int xi[3];
    int yi[3];
    int zi[3];
    fill_vector(xi,x,N);
    fill_vector(yi,y,N);
    fill_vector(zi,z,N);

    int max[N_SPECIES] = {0};
    int sum = 0;

    for (int xii : xi) {
        for (int yii : yi) {
            for (int zii : zi) {
                int species = grid[xii][yii][zii];
                if (species != 0){
                    //cout << endl << xii << " " << yii << " " << zii << endl;
                    max[species-1]++;
                    sum++;
                }    
            }
        }
    }
    return alive_or_dead(grid, max, sum, x, y ,z);
}

// Calculates the maximum of each species
void info_of_gen(int*** grid, int id, int p,int** maximum, int epoch, long long N)
{
    int size;
    double div = (float) N/p;
    // TODO NPROC TEM DE SER MENOR QUE N
    if (id >= N%p) size = floor(div);
    else size = ceil(div);
    int maxNew[N_SPECIES] = {0};
    int max[N_SPECIES] = {0};
    #pragma omp parallel for collapse(2) reduction(+:max)
    for (int x = 1; x < size+1; x++) {
        for (int y = 0; y < N; y++) {
            #pragma omp simd
            for (int z = 0; z < N; z++) {
                int species = grid[x][y][z];
                if (species != 0)
                    max[species-1]++;        
            }
        }
    }
    MPI_Reduce (max, maxNew, N_SPECIES, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(!id){
        compare_and_modify(maxNew, maximum, epoch);
    }}

// Simulates one generation
void gen_generation(int*** grid, int*** new_grid, int id, int p, int** maximum, MPI_Request** requests, int epoch, long long N)
{

    int size;
    double div = (float) N/p;
    // TODO NPROC TEM DE SER MENOR QUE N
    if (id >= N%p) size = floor(div);
    else size = ceil(div);
    int maxNew[N_SPECIES] = {0};
    int max[N_SPECIES] = {0};
    #pragma omp parallel for reduction(+:max)
    for (int x = 1; x < size+1; x++) {
        int xi = (x == 1) ? size : (x - 1);
        //cout << x << " " << xi << endl;
        for (int y = 0; y < N; y++) {
            #pragma omp simd
            for (int z = 0; z < N; z++) {
                int species = grid[xi][y][z];
                if (species != 0)
                    max[species-1]++;
                new_grid[xi][y][z] = visit_node(grid, N, xi, y, z);
            }         
        }
        if(x==1){
            int node = (id+1)%p;
            for (int i = 0; i< N;i++) {
                MPI_Isend(new_grid[size][i], N, MPI_INT, node, N+i, MPI_COMM_WORLD, &requests[1][i]);
                MPI_Irecv(new_grid[size+1][i], N, MPI_INT, node, i, MPI_COMM_WORLD, &requests[3][i]);    
            }
        } else if (x == 2) {
            int node =  (id + p - 1) % p;
            for (int i = 0; i< N;i++) {
                MPI_Isend(new_grid[1][i], N, MPI_INT, node, i, MPI_COMM_WORLD, &requests[0][i]);
                MPI_Irecv(new_grid[0][i], N, MPI_INT, node, N+i, MPI_COMM_WORLD, &requests[2][i]);        
            }
        }
            
    }
    MPI_Reduce (max, maxNew, N_SPECIES, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(!id){
        compare_and_modify(maxNew, maximum, epoch);
    }
    for (int i = 0; i< N;i++) {
        MPI_Wait(&requests[2][i], MPI_STATUS_IGNORE);
        MPI_Wait(&requests[3][i], MPI_STATUS_IGNORE);
    }
    for (int i = 0; i< N;i++) {
        MPI_Wait(&requests[0][i], MPI_STATUS_IGNORE);
        MPI_Wait(&requests[1][i], MPI_STATUS_IGNORE);
    }

}
// meter a contagem do primeiro fora do timer e aqui trocar o visit para o inicio e o max para depois (e usar o newgrid)

// Simulates all generations
void full_generation(int*** grid1, int*** grid2, int id, int p, int** maximum , MPI_Request** requests, int gens, long long N, float density, int seed){
    bool grid_to_use = true;

    for (int x = 0; x < gens; x++) {
        if(grid_to_use){
            gen_generation(grid1, grid2, id, p, maximum, requests, x, N);
            grid_to_use = false;
        }
        else {
            gen_generation(grid2, grid1, id, p, maximum, requests, x, N);
            grid_to_use = true;
        }
        //if (!id){
        //    cout << x << endl;
        //}
    }
    if(grid_to_use) info_of_gen(grid1, id, p, maximum, gens, N);
    else info_of_gen(grid2, id, p, maximum, gens, N);
}



int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        cout << "Usage: " << "<Generations> <N> <density> <seed>" << endl;
        return 1;
    }
    //cout << atoi(argv[1]) << " " << atoll(argv[2]) << " " << atof(argv[3]) << " " << atoi(argv[4]) << endl;
    MPI_Status status;
    double exec_time;
    int*** grid1 = NULL;
    int*** grid2 = NULL;
    int** maximum = NULL;
    //int** upper = NULL;
    //int** lower = NULL;
    int id, p, size, N;
    N = atoll(argv[2]);
    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    maximum = new int*[N_SPECIES];
    for (int i = 0; i < N_SPECIES; i++) {
        maximum[i] = new int[2];
        maximum[i][0] = 0;
        maximum[i][1] = 0;
    }

    grid1 = gen_grid(N, id, p);
    grid2 = gen_grid(N, id, p);
    gen_initial_grid(grid1, N, atof(argv[3]), atoi(argv[4]), id, p);
    
    double div = (float) N/p;
    // TODO NPROC TEM DE SER MENOR QUE N
    if (id >= N%p) size = floor(div);
    else size = ceil(div);
    MPI_Request** requests;
    requests = new MPI_Request*[4];
    for (int i = 0; i < 4; i++) {
        requests[i] = new MPI_Request[N];
    }
    //cout << (id-1)%p << " " << id << endl;
    for (int i = 0; i< N;i++) {
        MPI_Isend(grid1[1][i], N, MPI_INT, (id + p - 1) % p, i, MPI_COMM_WORLD, &requests[0][i]);
        MPI_Isend(grid1[size][i], N, MPI_INT, (id+1)%p, N+i, MPI_COMM_WORLD, &requests[1][i]);
        MPI_Irecv(grid1[0][i], N, MPI_INT, (id + p - 1) % p, N+i, MPI_COMM_WORLD, &requests[2][i]);        
        MPI_Irecv(grid1[size+1][i], N, MPI_INT, (id+1)%p, i, MPI_COMM_WORLD, &requests[3][i]);    
    }
    for (int i = 0; i< N;i++) {
        MPI_Wait(&requests[2][i], MPI_STATUS_IGNORE);
        MPI_Wait(&requests[3][i], MPI_STATUS_IGNORE);
    }
    for (int i = 0; i< N;i++) {
        MPI_Wait(&requests[0][i], MPI_STATUS_IGNORE);
        MPI_Wait(&requests[1][i], MPI_STATUS_IGNORE);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    exec_time = - MPI_Wtime();
    full_generation(grid1, grid2, id, p, maximum, requests, atoi(argv[1]), N, atof(argv[3]), atoi(argv[4]));
    exec_time += MPI_Wtime();
    if(!id){
        fprintf(stderr, "%.1fs\n", exec_time);
        fflush(stderr);
        print_results(maximum);
    }
    MPI_Finalize();
    return 0;
}