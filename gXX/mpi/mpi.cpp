#include <iostream>
#include <vector>
#include <cstdlib>
#include <omp.h>
#include <cmath>
#include <mpi.h>

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
void gen_initial_grid(int*** grid, long long N, float density, int input_seed)
{

    init_r4uni(input_seed);
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                if(r4_uni() < density) {
                    grid[x][y][z] = (int)(r4_uni() * N_SPECIES) + 1;
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
        if(grid[0] == NULL) {
            cout <<"Failed to allocate matrix\n";
            exit(1);
        }
    }
    return grid;
}

int*** gen_grid(long long N, int id, int p) {
    int*** grid;
    int size;

    if (id == p-1 && N%p != 0) size = floor(N/p);
    else size = ceil(N/p);

    grid = (int ***) malloc(size * sizeof(int **));

    if(grid == NULL) {
        cout <<"Failed to allocate matrix\n";
        exit(1);
    }
    for(int x = 0; x < size; x++) {
        grid[x] = gen_plaquette(N);
    }
    return grid;
}
void copyFromOrig(int*** orig, int*** grid, int id, int p, long long N) {
    int start = id * (N / p);
    int end;
    if (id == p - 1) end = N; 
    else end = start + (N / p);

    for (int x = start; x < end; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                grid[x - start][y][z] = orig[x][y][z];
            }
        }
    }
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
                    max[species-1]++;
                    sum++;
                }    
            }
        }
    }
    return alive_or_dead(grid, max, sum, x, y ,z);
}

// Calculates the maximum of each species
void info_of_gen(int*** grid, int** maximum, int epoch, long long N)
{
    int max[N_SPECIES] = {0};
    #pragma omp parallel for collapse(2) reduction(+:max)
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            #pragma omp simd
            for (int z = 0; z < N; z++) {
                int species = grid[x][y][z];
                if (species != 0)
                    max[species-1]++;        
            }
        }
    }
    compare_and_modify(max, maximum, epoch);
}

// Simulates one generation
void gen_generation(int*** grid, int*** new_grid, int** maximum, int epoch, long long N)
{
    int max[N_SPECIES] = {0};
    #pragma omp parallel for collapse(2) reduction(+:max)
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            #pragma omp simd
            for (int z = 0; z < N; z++) {
                int species = grid[x][y][z];
                if (species != 0)
                    max[species-1]++;
                new_grid[x][y][z] = visit_node(grid, N, x, y, z);
            }
        }
    }
    compare_and_modify(max, maximum, epoch);
}
// meter a contagem do primeiro fora do timer e aqui trocar o visit para o inicio e o max para depois (e usar o newgrid)

// Simulates all generations
void full_generation(int*** grid1, int*** grid2, int** maximum, int gens, long long N, float density, int seed){
    bool grid_to_use = true;

    for (int x = 0; x < gens; x++) {
        if(grid_to_use){
            gen_generation(grid1, grid2, maximum, x, N);
            grid_to_use = false;
        }
        else {
            gen_generation(grid2, grid1, maximum, x, N);
            grid_to_use = true;
        }
        cout << x << endl;
    }
    if(grid_to_use) info_of_gen(grid1, maximum, gens, N);
    else info_of_gen(grid2, maximum, gens, N);
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
    int*** orig = NULL;
    int*** grid1 = NULL;
    int*** grid2 = NULL;
    int** maximum = NULL;
    int** upper = NULL;
    int** lower = NULL;
    int id = 0;
    int p = 1;
    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);


    orig = gen_grid(atoll(argv[2]), 0, 1);
    maximum = new int*[N_SPECIES];
    for (int i = 0; i < N_SPECIES; i++) {
        maximum[i] = new int[2];
        maximum[i][0] = 0;
        maximum[i][1] = 0;
    }
    orig = gen_grid(atoll(argv[2]), 0, 1);
    gen_initial_grid(orig, atoll(argv[2]), atof(argv[3]), atoi(argv[4]));

    grid1 = gen_grid(atoll(argv[2]), id, p);
    grid2 = gen_grid(atoll(argv[2]), id, p);
    // maybe usar pares de plaquettes para se puder ir trocando tipo os grids
    upper = gen_plaquette(atoll(argv[2]));
    lower = gen_plaquette(atoll(argv[2]));
    copyFromOrig(orig, grid1, id, p, atoll(argv[2]));
    MPI_Barrier (MPI_COMM_WORLD);
    exec_time = - MPI_Wtime();
    full_generation(grid1, grid2, maximum, atoi(argv[1]), atoll(argv[2]), atof(argv[3]), atoi(argv[4]));
    exec_time += MPI_Wtime();
    MPI_Barrier (MPI_COMM_WORLD);
    if(!id){
        fprintf(stderr, "%.1fs\n", exec_time);
        print_results(maximum);
    }
    MPI_Finalize();
    return 0;
}