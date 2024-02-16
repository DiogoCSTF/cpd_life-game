#include <iostream>
#include <vector>
#include <cstdlib>

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

// Compare and modify function
void compare_and_modify(vector<int> &max, vector<int> &maximum)
{
    for (int i = 0; i < N_SPECIES; i++) {
        if (max[i] > maximum[i]) {
            maximum[i] = max[i];
        }
    }
}

// Calculates the maximum of each species
void info_of_gen(vector<vector<vector<int>>> &grid, vector<int> &maximum, long long N)
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
    compare_and_modify(max, maximum);
}
/*
Duvidas:
- tem de se mudar todos os x,y,z para long long???



*/


// Visits node
int visit_node(vector<vector<vector<int>>> &grid, int x, int y, int z) {
    return 0;
}

// Simulates one generation
vector<vector<vector<int>>> gen_generation(vector<vector<vector<int>>> &grid, long long N)
{
    vector<vector<vector<int>>> new_grid(N, vector<vector<int>>(N, vector<int>(N)));
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                new_grid[x][y][z] = visit_node(grid, x, y, z);         
            }
        }
    }
    return new_grid;
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

// Simulates all generations
vector<vector<vector<int>>> full_generation(long long gens, long long N, float density, int seed)
{
    vector<int> maximum(N_SPECIES, 0);
    vector<vector<vector<int>>> grid = gen_initial_grid(N, density, seed);
    print_grid(grid,N);
    for (int x = 0; x < gens; x++) {
        info_of_gen(grid, maximum, N);
        for (int i = 0; i < N_SPECIES; i++) {
            cout << maximum[i] << " ";
        }
        cout << endl;
        grid = gen_generation(grid, N);
    }
    return grid;
}



int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        cout << "Usage: " << "<Generations> <N> <density> <seed>" << endl;
        return 1;
    }
    cout << atoi(argv[1]) << " " << atoll(argv[2]) << " " << atof(argv[3]) << " " << atoi(argv[4]) << endl;

    full_generation(atoi(argv[1]), atoll(argv[2]), atof(argv[3]), atoi(argv[4]));
    return 0;
}