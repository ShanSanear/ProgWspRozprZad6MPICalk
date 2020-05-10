#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <vector>
#include <iostream>
#include <cmath>

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
    int n;

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
    this->n = n;
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
    MPI_Init(&argc, &argv);
    int node, Nnodes, tag = 2;
    int finalized = 3;
    double data = 0;
    int rec, recv;

    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    MPI_Comm_size(MPI_COMM_WORLD, &Nnodes);
    struct ToCalculate toBeCalculated[Nnodes];
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
            toBeCalculated[i].a = a;
            toBeCalculated[i].b = b;
            toBeCalculated[i].c = c;
            toBeCalculated[i].end = start + delta * (i+1);
            toBeCalculated[i].start = start + delta * i;
            double i_st = start + delta * i;
            double i_end = start + delta * (i+1);
            IntegralCalculator calculator(toBeCalculated[i]);
            result += calculator.calculate_trapezoidal_rule();
        }
        printf("Result: %f\n", result);
        
        // for (int i = 1; i < 4; i++)
        // {
        //     MPI_Send(&data, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        // }
        // recv = 0;
        // for (int i = 1; i < 4; i++){
        //     MPI_Recv(&rec, 1, MPI_INT, i, finalized, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     recv += rec;
        // }
        // printf("FInished value is: %d", recv);
    }
    else if (node > 0)
    {
        printf("Node: %d\n", node);
        // MPI_Recv(&data, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD,
        //          MPI_STATUS_IGNORE);
        // printf("Wezel %d odebral %f od wezla 0\n", node, data);
        // printf("Wezel %d value: %d\n", node, f);
        // MPI_Send(&f, 1, MPI_INT, 0, finalized, MPI_COMM_WORLD);
    }
    else
    {
        printf("Node in else: %d\n", node);
    }
    MPI_Finalize();
    return 0;
}
