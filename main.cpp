#include <iostream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <ctime>
#include <mpi.h>

using namespace std;
double A, B, C;
int N;
double startTime, endTime, parallelTimeTaken, timeSingle;

struct CalculateParameters
{
    double A;
    double B;
    double C;
};

double calculate_formula(double x)
{
    return (A * (x * x) + B * x + C);
}

double get_double_from_stdin(const char *message)
{
    double out;
    std::cout << message << std::endl;
    std::cin >> out;
    return out;
}

int get_int_from_stdin(const char *message)
{
    int out;
    std::cout << message << std::endl;
    std::cin >> out;
    return out;
}

int main()
{
    double sum, dx, start_point, end_point, result;
    double functionParameters[3] = {0};
    double functionLimits[2] = {0};

    int node, numOfNodes;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    cout << setprecision(3) << fixed;
    int lengths[4] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint displacements[3] = {0, sizeof(double), sizeof(double) * 2};

    MPI_Datatype mpiCalculateParametersDatatype;
    MPI_Type_create_struct(3, lengths, displacements, types, &mpiCalculateParametersDatatype);
    MPI_Type_commit(&mpiCalculateParametersDatatype);

    if (node == 0)
    {

        A = get_double_from_stdin("Specify a:");
        B = get_double_from_stdin("Specify b:");
        C = get_double_from_stdin("Specify c:");
        N = get_int_from_stdin("Specifiy N:");
        start_point = get_double_from_stdin("Specify start point:");
        end_point = get_double_from_stdin("Specify end point:");

        sum = 0;
        dx = (end_point - start_point) / N;

        startTime = MPI_Wtime();
        for (int i = 0; i <= N; i++)
        {
            sum += calculate_formula(start_point + i * dx);
        }

        sum *= dx;

        endTime = MPI_Wtime();

        timeSingle = endTime - startTime;
        printf("Single node result: %f\n", sum);
        printf("Single node time: %f second\n", timeSingle);
        struct CalculateParameters calculateStruct;
        //There is no MPI_Barrier here because only one node can reach this section
        startTime = MPI_Wtime();
        functionLimits[0] = start_point;
        functionLimits[1] = end_point;
        calculateStruct.A = A;
        calculateStruct.B = B;
        calculateStruct.C = C;
        for (int i = 1; i < numOfNodes; i++)
        {
            MPI_Send(&calculateStruct, 1, mpiCalculateParametersDatatype, i, 0, MPI_COMM_WORLD);
            MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(functionLimits, 2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        sum = 0;
        dx = (end_point - start_point) / N;

        for (int i = 0; i <= N / numOfNodes; i++)
        {
            sum += calculate_formula(start_point + i * dx);
        }

        sum *= dx;
    }
    else
    {
        struct CalculateParameters calcStruct;
        // MPI_Recv(functionParameters, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&calcStruct, 1, mpiCalculateParametersDatatype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // printf("What is inside? %f, %f, %f, %d\n", calcStruct.A, calcStruct.B, calcStruct.C, calcStruct.N);
        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(functionLimits, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        start_point = functionLimits[0];
        end_point = functionLimits[1];
        A = calcStruct.A;
        B = calcStruct.B;
        C = calcStruct.C;

        sum = 0;
        dx = (end_point - start_point) / N;

        for (int i = node * (N / numOfNodes); i <= node * N / numOfNodes + N / numOfNodes; i++)
        {
            sum += calculate_formula(start_point + i * dx);
        }
        int check = ((node * (N / numOfNodes)) + (N / numOfNodes));
        int check2 = node * N / numOfNodes + N / numOfNodes;
        printf("Check 1: %d Check 2: %d\n", check, check2);
        sum *= dx;
    }

    MPI_Reduce(&sum, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (node == 0)
    {
        endTime = MPI_Wtime();
        parallelTimeTaken = endTime - startTime;
        printf("Parallized result: %f\n", result);
        printf("Parallized time: %f seconds\n", parallelTimeTaken);
    }
    MPI_Type_free(&mpiCalculateParametersDatatype);
    MPI_Finalize();
    return 0;
}