#include <iostream>
#include <vector>
#include <cstdlib>
#include <omp.h>

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
vector<vector<vector<int>>> gen_initial_grid(long long N, float density, int input_seed)
{
    vector<vector<vector<int>>> grid(N, vector<vector<int>>(N, vector<int>(N)));

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
    return grid;
}


// Prints the grid
void print_grid(vector<vector<vector<int>>> &grid, long long N) {
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
void print_results(vector<vector<int>> &maximum){
    for(int i = 0; i < N_SPECIES; i++){
        cout << i+1 << " " << maximum[i][0] << " " << maximum[i][1] << endl;
    }
}

// Compare and modify function
void compare_and_modify(vector<int> &max, vector<vector<int>> &maximum, int epoch)
{
    for (int i = 0; i < N_SPECIES; i++) {
        if (max[i] > maximum[i][0]) {
            maximum[i][0] = max[i];
            maximum[i][1] = epoch;
        }
    }
}

// Calculates the maximum of each species
void info_of_gen(vector<vector<vector<int>>> &grid, vector<vector<int>> &maximum, int epoch, long long N)
{
    vector<int> max(N_SPECIES, 0);
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                int species = grid[x][y][z];
                if (species != 0)
                    max[species-1]++;        
            }
        }
    }
    compare_and_modify(max, maximum, epoch);
}
/*
Duvidas:
- tem de se mudar todos os x,y,z para long long???



*/

// Fills vectors
void fill_vector(vector<int> &dimentions, int dimention, long long N) {
    if (dimention == 0) dimentions[0] = N-1;
    else dimentions[0] = dimention-1;

    dimentions[1] = dimention;

    if (dimention == N-1) dimentions[2] = 0;
    else dimentions[2] = dimention+1;
}

int alive_neighbour(vector<int> &max) {
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

int alive_or_dead(vector<vector<vector<int>>> &grid, vector<int> &max, int sum, int x, int y, int z) {
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
int visit_node(vector<vector<vector<int>>> &grid, long long N, int x, int y, int z) {
    vector<int> xi = vector<int>(3);
    vector<int> yi = vector<int>(3);
    vector<int> zi = vector<int>(3);
    fill_vector(xi,x,N);
    fill_vector(yi,y,N);
    fill_vector(zi,z,N);

    vector<int> max(N_SPECIES, 0);
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

// Simulates one generation
vector<vector<vector<int>>> gen_generation(vector<vector<vector<int>>> &grid, long long N)
{
    vector<vector<vector<int>>> new_grid(N, vector<vector<int>>(N, vector<int>(N)));
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                new_grid[x][y][z] = visit_node(grid, N, x, y, z);         
            }
        }
    }
    return new_grid;
}


// Simulates all generations
void full_generation(vector<vector<vector<int>>> &grid, vector<vector<int>> &maximum, long long gens, long long N, float density, int seed)
{
    for (int x = 0; x < gens; x++) {
        info_of_gen(grid, maximum, x, N);
        grid = gen_generation(grid, N);
        cout << x << endl;
    }
    info_of_gen(grid, maximum, gens, N);
}



int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        cout << "Usage: " << "<Generations> <N> <density> <seed>" << endl;
        return 1;
    }
    //cout << atoi(argv[1]) << " " << atoll(argv[2]) << " " << atof(argv[3]) << " " << atoi(argv[4]) << endl;
    double exec_time;
    vector<vector<vector<int>>> grid = gen_initial_grid(atoll(argv[2]), atof(argv[3]), atoi(argv[4]));
    vector<vector<int>> maximum(N_SPECIES, vector<int>(2,0));

    exec_time = -omp_get_wtime();
    full_generation(grid, maximum, atoi(argv[1]), atoll(argv[2]), atof(argv[3]), atoi(argv[4]));
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_results(maximum);

}