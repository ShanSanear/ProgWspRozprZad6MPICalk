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
    printf("Our values: %f, %f, %f, %f, %f, %d, %f\n",
    a, b, c, starting_point, ending_point, n, dx);
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
    MPI_Init(&argc, &argv);
    int lengths[5] = {1,1,1,1,1};
    MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint displacements[5] = {0, sizeof(double), sizeof(double) * 2, sizeof(double) * 3,
    sizeof(double) * 4};
    MPI_Datatype ToCalculateMpi;
    MPI_Type_create_struct(5, lengths, displacements, types, &ToCalculateMpi);
    MPI_Type_commit(&ToCalculateMpi);
    int node, Nnodes, tag = 2;
    int finalized = 3;
    double data = 0;
    int rec, recv;

    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    MPI_Comm_size(MPI_COMM_WORLD, &Nnodes);
    
    if (node == 0)
    {
        double a = 1.0;
        double b = 2.0;
        double c = 3.0;
        int n = 1000;
        double start = -5.0;
        double end = 5.0;
        double result = 0.0;
        double delta = (end - start)/Nnodes;
        for (int i = 0; i< Nnodes; i++) {
            struct ToCalculate toBeCalculated;
            toBeCalculated.a = a;
            toBeCalculated.b = b;
            toBeCalculated.c = c;
            toBeCalculated.start = start + delta * i;
            toBeCalculated.end = start + delta * (i+1);
            double i_st = start + delta * i;
            double i_end = start + delta * (i+1);
            IntegralCalculator calculator(toBeCalculated);
            result += calculator.calculate_trapezoidal_rule();
        }
        printf("Result: %f\n", result);
        
        }
    
    if (node > 0)
    {
        printf("Node: %d\n", node);
    }
    else
    {
        printf("Node in else: %d\n", node);
    }
    MPI_Finalize();
    return 0;
}
