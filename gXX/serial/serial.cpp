#include <iostream>
#include <vector>
#include <cstdlib>

#define N_SPECIES 9
using namespace std;

unsigned int seed;

void init_r4uni(int input_seed)
{
    seed = input_seed + 987654321;
}

float r4_uni()
{
    int seed_in = seed;

    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);

    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}


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

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        cout << "Usage: " << "<Generations> <N> <density> <seed>" << endl;
        return 1;
    }
    int generations = atoi(argv[1]);
    long long N = atoll(argv[2]);
    float density = atof(argv[3]);
    int seed = atoi(argv[4]);
    cout << generations << " " << N << " " << density << " " << seed << endl;

    vector<vector<vector<int>>> grid = gen_initial_grid(N, density, seed);~
    /*
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
    */
    
    return 0;
}