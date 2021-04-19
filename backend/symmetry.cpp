#include <iostream>
#include "spglib.h"

void test_spg_get_symmetry(void)
{
    double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
    double position[][3] =
        {
            {0, 0, 0},
            {0.5, 0.5, 0.25},
            {0.3, 0.3, 0},
            {0.7, 0.7, 0},
            {0.2, 0.8, 0.25},
            {0.8, 0.2, 0.25},
            {0, 0, 0.5},
            {0.5, 0.5, 0.75},
            {0.3, 0.3, 0.5},
            {0.7, 0.7, 0.5},
            {0.2, 0.8, 0.75},
            {0.8, 0.2, 0.75}};
    int types[] = {1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2};
    int num_atom = 12;
    int max_size = 50;
    int i, j, size;
    int rotation[max_size][3][3];
    double translation[max_size][3];

    double origin_shift[3] = {0.1, 0.1, 0};
    for (i = 0; i < num_atom; i++)
    {
        for (j = 0; j < 3; j++)
        {
            position[i][j] += origin_shift[j];
        }
    }

    printf("*** Example of spg_get_symmetry (Rutile two unit cells) ***:\n");
    size = spg_get_symmetry(rotation,
                            translation,
                            max_size,
                            lattice,
                            position,
                            types,
                            num_atom,
                            1e-5);
    for (i = 0; i < size; i++)
    {
        printf("--- %d ---\n", i + 1);
        for (j = 0; j < 3; j++)
            printf("%2d %2d %2d\n", rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
        printf("%f %f %f\n",
               translation[i][0], translation[i][1], translation[i][2]);
    }
}

// static void test_spg_get_symmetry(void)
// {
//     double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
//     double position[][3] =
//         {
//             {0, 0, 0},
//             {0.5, 0.5, 0.25},
//             {0.3, 0.3, 0},
//             {0.7, 0.7, 0},
//             {0.2, 0.8, 0.25},
//             {0.8, 0.2, 0.25},
//             {0, 0, 0.5},
//             {0.5, 0.5, 0.75},
//             {0.3, 0.3, 0.5},
//             {0.7, 0.7, 0.5},
//             {0.2, 0.8, 0.75},
//             {0.8, 0.2, 0.75}};
//     int types[] = {1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2};
//     int num_atom = 12;
//     int max_size = 50;
//     int i, j, size;
//     int rotation[max_size][3][3];
//     double translation[max_size][3];

//     double origin_shift[3] = {0.1, 0.1, 0};
//     for (i = 0; i < num_atom; i++)
//     {
//         for (j = 0; j < 3; j++)
//         {
//             position[i][j] += origin_shift[j];
//         }
//     }

//     printf("*** Example of spg_get_symmetry (Rutile two unit cells) ***:\n");
//     size = spg_get_symmetry(rotation,
//                             translation,
//                             max_size,
//                             lattice,
//                             position,
//                             types,
//                             num_atom,
//                             1e-5);
//     for (i = 0; i < size; i++)
//     {
//         printf("--- %d ---\n", i + 1);
//         for (j = 0; j < 3; j++)
//             printf("%2d %2d %2d\n", rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
//         printf("%f %f %f\n",
//                translation[i][0], translation[i][1], translation[i][2]);
//     }
// }