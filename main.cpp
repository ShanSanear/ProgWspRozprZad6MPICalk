#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <vector>
#include <iostream>
#include <cmath>
int n = 1000;

struct ToCalculate {
    double a;
    double b;
    double c;
    double start;
    double end;
};

double get_double_from_stdin(const char* message) {
    double out;
    std::cout << message << std::endl;
    std::cin >> out;
    printf("Returning: %f\n", out);
    return out;
}
class IntegralCalculator {
private:
    double a, b, c;
    double starting_point, ending_point;
    double dx;
    
    double calculate_formula(double x) const;

public:
    IntegralCalculator(ToCalculate toBeCalculated);

    double calculate_rectangular_method();

    double calculate_trapezoidal_rule();

};


IntegralCalculator::IntegralCalculator(ToCalculate toBeCalculated) {
    this->a = toBeCalculated.a;
    this->b = toBeCalculated.b;
    this->c = toBeCalculated.c;
    this->starting_point = toBeCalculated.start;
    this->ending_point = toBeCalculated.end;
    this->dx = (ending_point - starting_point) / n;
}

double IntegralCalculator::calculate_formula(double x) const {
    return a * powl(x, 2) + b * x + c;
}

double IntegralCalculator::calculate_rectangular_method() {
    printf("Calculating using rectangular method\n");
    double calculated_integral = 0.0;
    
    for (int i = 1; i < n; i++) {
        calculated_integral += calculate_formula(starting_point + i * dx) * dx;
    }
    
    
    printf("Calculated integral for limits from %f to %f is: %.15f\n", starting_point, ending_point,
           calculated_integral);

    return calculated_integral;
}

double IntegralCalculator::calculate_trapezoidal_rule() {
    printf("Calculating using trapezoidal rule\n");
    double calculated_integral = 0.0;
    
    double value_for_start = calculate_formula(starting_point);
    double value_for_end = calculate_formula(ending_point);
    double average = (value_for_start + value_for_end) / 2;
    for (int i = 1; i < n; i++) {
        calculated_integral += (calculate_formula(starting_point + (i * dx)));
    }
    calculated_integral += average;
    calculated_integral *= dx;
    
    
    printf("Calculated integral for limits from %f to %f is: %.15f\n", starting_point, ending_point,
           calculated_integral);
    return calculated_integral;
}

int main(int argc, char *argv[])
{
    double a,b,c,start,end;
    
    MPI_Init(&argc, &argv);
    int lengths[5] = {1,1,1,1,1};
    MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint displacements[5] = {0, sizeof(double), sizeof(double) * 2, sizeof(double) * 3,
    sizeof(double) * 4};
    MPI_Datatype toCalculateMpiDatatype;
    MPI_Type_create_struct(5, lengths, displacements, types, &toCalculateMpiDatatype);
    MPI_Type_commit(&toCalculateMpiDatatype);
    int node, Nnodes, tag = 2;
    int finalized = 3;
    double data = 0;
    int rec, recv;

    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    MPI_Comm_size(MPI_COMM_WORLD, &Nnodes);
    if (node == 0) {
        a = get_double_from_stdin("Specify a:");
        b = get_double_from_stdin("Specify b:");
        c = get_double_from_stdin("Specify c:");
        start = get_double_from_stdin("Specify start point:");
        end = get_double_from_stdin("Specify end point:");
    }
    double delta = (end - start)/(Nnodes);
    double localSum = 0.0;
    double result = 0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    double time = MPI_Wtime();
    if (node == 0)
    {
        for (int i = 0; i < Nnodes; i++) {
            struct ToCalculate toBeCalculated;
            toBeCalculated.a = a;
            toBeCalculated.b = b;
            toBeCalculated.c = c;
            toBeCalculated.start = start + delta * i;
            toBeCalculated.end = start + delta * (i+1);
            double i_st = start + delta * i;
            double i_end = start + delta * (i+1);
            printf("Sending to: %d\n", i+1);
            printf("%d / %d \n", node, Nnodes);
            if (i+1 != Nnodes) {
            MPI_Send(&toBeCalculated, 1, toCalculateMpiDatatype, i+1, tag, MPI_COMM_WORLD);
            } else {
                IntegralCalculator calc(toBeCalculated);
        localSum += calc.calculate_trapezoidal_rule();
            }
        }
        printf("Ending with node 0\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (node > 0)
    {
        printf("Node: %d\n", node);
        struct ToCalculate receivedCalculate;
        MPI_Recv(&receivedCalculate, 1, toCalculateMpiDatatype, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        IntegralCalculator calc(receivedCalculate);
        localSum += calc.calculate_trapezoidal_rule();

    }
    MPI_Reduce(&localSum, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (node == 0) {
        time = MPI_Wtime() - time;
        printf("Calculating took: %f seconds\n", time);
        printf("Result: %f\n", result);
    }
    
    MPI_Type_free(&toCalculateMpiDatatype);
    MPI_Finalize();
    return 0;
}
