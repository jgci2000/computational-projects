#include <iostream>
#include <cmath>
#include <map>

#include "read_NN.h"
#include "rng.h"

double dot_product(double *vector_a, double *vector_b, int size) {
    double product = 0;
    for (int i = 0; i < size; i++) {
        product += vector_a[i] * vector_b[i];
    }
    return product;
}

int main(int argc, char **argv) {
    // System parameters
    int L = 4;
    int N_atm = pow(L, 3);
    int J = 1;
    int NN = 6;
    int *NN_table = new int[N_atm * NN];
    read_NN_talbe("NN_tables/NN_3D_L" + std::to_string(L) + ".txt", NN_table);

    // Metropolis parameters
    int NT = 10;
    double temperatures[NT];
    printf("temperatures: ");
    for (int i = 0; i < NT; ++i) {
        temperatures[i] = i * (3.0 - 0.01) / (NT - 1) + 0.01;
        printf("%lf, ", temperatures[i]);
    }
    printf("\n");
    int steps = 50000;
    int measure = 40000;
    int n_measure = steps - measure;

    double *M = new double[NT];
    double *E = new double[NT];

    // Initialization of configuration and RNG
    RNG rng((int) time(NULL));
    RNG rng1(1);
    RNG rng2(2);

    double x1;
    double x2;

    std::map<int, double*> spin_lattice;
    for (int i = 0; i < N_atm; ++i) {
        do {
            x1 = rng1.rand_uniform() * 2 - 1;
            x2 = rng2.rand_uniform() * 2 - 1;
        } while(x1*x1 + x2*x2 >= 1);

        spin_lattice[i] = new double[3];
        spin_lattice[i][0] = 2 * x1 * sqrt(1 - x1*x1 - x2*x2);
        spin_lattice[i][1] = 2 * x2 * sqrt(1 - x1*x1 - x2*x2);
        spin_lattice[i][2] = 1 - 2 * (x1*x1 + x2*x2);
    }

    // Metropolis Sampling
    for (int t_idx = 0; t_idx < NT; ++t_idx) {
        double T = temperatures[t_idx];
        double M_tmp[n_measure];
        double E_tmp[n_measure];
        int measure_idx = 0;

        for (int t = 0; t < steps; ++t) {
            do {
                x1 = rng1.rand_uniform() * 2 - 1;
                x2 = rng2.rand_uniform() * 2 - 1;
            } while(x1*x1 + x2*x2 >= 1);

            double new_spin[3];
            new_spin[0] = 2 * x1 * sqrt(1 - x1*x1 - x2*x2);
            new_spin[1] = 2 * x2 * sqrt(1 - x1*x1 - x2*x2);
            new_spin[2] = 1 - 2 * (x1*x1 + x2*x2);

            int site = rng.rand() % N_atm;

            double sum_NNx = 0, sum_NNy = 0, sum_NNz = 0;
            for (int a = 0; a < NN; ++a) {
                sum_NNx += spin_lattice[NN_table[site * NN + a]][0];
                sum_NNy += spin_lattice[NN_table[site * NN + a]][1];
                sum_NNz += spin_lattice[NN_table[site * NN + a]][2];
            }

            double delta_E = - J * ((new_spin[0] - spin_lattice[site][0]) * sum_NNx + 
                                    (new_spin[1] - spin_lattice[site][1]) * sum_NNy + 
                                    (new_spin[2] - spin_lattice[site][2]) * sum_NNz);

            if (delta_E <= 0 || rng.rand_uniform() <= exp(- delta_E / T)) {
                spin_lattice[site][0] = new_spin[0];
                spin_lattice[site][1] = new_spin[1];
                spin_lattice[site][2] = new_spin[2];
            }

            if (t > measure) {
                // Magnetization
                double Mx = 0, My = 0, Mz = 0;
                for (int i = 0; i < N_atm; ++i) {
                    Mx += spin_lattice[i][0];
                    My += spin_lattice[i][1];
                    Mz += spin_lattice[i][2];
                }
                M_tmp[measure_idx] = sqrt(Mx*Mx + My*My + Mz*Mz) / N_atm;

                // Energy
                for (int i = 0; i < N_atm; ++i) {
                    for (int a = 0; a < NN; ++a) {
                        E_tmp[measure_idx] += dot_product(spin_lattice[site], spin_lattice[NN_table[i * NN + a]], 3);
                    }
                }
                E_tmp[measure_idx] = - J * E_tmp[measure_idx] / (2.0 * N_atm);

                measure_idx++;
            }
        }

        for (int i = 0; i < n_measure; ++i) {
            M[t_idx] += M_tmp[i] / n_measure;
            E[t_idx] += E_tmp[i] / n_measure;
        }

        printf("T = %d/%d \n", t_idx + 1, NT);
    }

    printf("M: ");
    for (int i = 0; i < NT; ++i) {
        printf("%lf, ", M[i]);
    }
    printf("\n");
    
    printf("E: ");
    for (int i = 0; i < NT; ++i) {
        printf("%lf, ", E[i]);
    }
    printf("\n");


    for (int i = 0; i < N_atm; ++i) {
        delete[] spin_lattice[i];
    }
    delete[] NN_table;
    delete[] M;
    delete[] E;

    return 0;
}

