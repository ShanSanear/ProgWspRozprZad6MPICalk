#include <iostream>
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <stdio.h>
#include <ctime>
#include <mpi.h>
#include <plog/Log.h>
#include <plog/Appenders/ConsoleAppender.h>
#include <cmath>

double A, B, C;

struct CalculateParameters
{
    double A;
    double B;
    double C;
    int N_quadratic;
    double start_point;
    double end_point;
    int N_Pi;
};

double calculate_quadratic_formula(double x)
{
    return (A * (x * x) + B * x + C);
}

double calculate_pi_approx(double x) {
    return (1 / (1 + std::pow(x, 2)));
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
    int standardPrecision = 6;
    int piPrecision = 20;
    int N_quadratic, N_Pi;
    const int struct_size = 7;
    double sum_quadratic, sum_pi, dx, dt, start_point, end_point, result_quadratic, result_pi;
    double startTime, endTime, parallelTimeTaken, timeSingle;
    plog::RollingFileAppender<plog::CsvFormatter> fileAppender("Datalogger.txt", 1048576, 5);
    plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;
    plog::init(plog::info, &fileAppender).addAppender(&consoleAppender);
    int node, numOfNodes;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    int lengths[struct_size] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[struct_size] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER, MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER};
    MPI_Aint displacements[struct_size] = {
        offsetof(CalculateParameters, A),
        offsetof(CalculateParameters, B),
        offsetof(CalculateParameters, C),
        offsetof(CalculateParameters, N_quadratic),
        offsetof(CalculateParameters, start_point),
        offsetof(CalculateParameters, end_point),
        offsetof(CalculateParameters, N_Pi)
        };

    MPI_Datatype mpiCalculateParametersDatatype;
    MPI_Type_create_struct(struct_size, lengths, displacements, types, &mpiCalculateParametersDatatype);
    MPI_Type_commit(&mpiCalculateParametersDatatype);
    std::stringstream resultStream;
    resultStream << std::fixed << std::setprecision(standardPrecision);
    if (node == 0)
    {
        PLOG_INFO << "Getting input values";
        A = get_double_from_stdin("Specify a:");
        B = get_double_from_stdin("Specify b:");
        C = get_double_from_stdin("Specify c:");
        start_point = get_double_from_stdin("Specify start point for quadratic function:");
        end_point = get_double_from_stdin("Specify end point for quadratic function:");
        N_quadratic = get_int_from_stdin("Specify N for quadratic function:");
        N_Pi = get_int_from_stdin("Specify N for Pi:");

        sum_quadratic = 0;
        dx = (end_point - start_point) / N_quadratic;

        startTime = MPI_Wtime();
        PLOG_INFO << "Calculating integral on single node";
        for (int i = 0; i <= N_quadratic; i++)
        {
            sum_quadratic += calculate_quadratic_formula(start_point + i * dx);
        }

        sum_quadratic *= dx;

        
        PLOG_INFO << "Calculating Pi on single node";
        sum_pi = 0;
        dt = 1.0 / N_Pi;
        for (int i = 0; i <= N_Pi; i++) {
            sum_pi += calculate_pi_approx(dt * i);
        }
        sum_pi *= dt * 4;
        endTime = MPI_Wtime();
        timeSingle = endTime - startTime;
        resultStream << "Single node integral result: " << static_cast<double>(sum_quadratic);
        PLOG_INFO << resultStream.str();
        resultStream.str(std::string());
        resultStream << std::setprecision(piPrecision);
        resultStream << "Sinlge node pi: " << static_cast<double>(sum_pi);

        PLOG_INFO << resultStream.str();
        PLOG_INFO << "Single node time: " << timeSingle << " second(s)";
        resultStream << std::setprecision(standardPrecision);
        struct CalculateParameters calculateStruct;
        //There is no MPI_Barrier here because only one node can reach this section
        startTime = MPI_Wtime();
        calculateStruct.A = A;
        calculateStruct.B = B;
        calculateStruct.C = C;
        calculateStruct.N_quadratic = N_quadratic;
        calculateStruct.start_point = start_point;
        calculateStruct.end_point = end_point;
        calculateStruct.N_Pi = N_Pi;
        PLOG_INFO << "Sending data to other nodes";
        for (int i = 1; i < numOfNodes; i++)
        {
            MPI_Send(&calculateStruct, 1, mpiCalculateParametersDatatype, i, 0, MPI_COMM_WORLD);
            
        }

        sum_quadratic = 0;
        sum_pi = 0;
        PLOG_INFO << "Calculating integral, node: " << node;
        for (int i = 0; i <= N_quadratic / numOfNodes; i++)
        {
            sum_quadratic += calculate_quadratic_formula(start_point + i * dx);
        }

        sum_quadratic *= dx;

        PLOG_INFO << "Calculating Pi, node: " << node;
        for (int i = 0; i <= N_Pi / numOfNodes; i++) {
            sum_pi += calculate_pi_approx(dt * i);
        }
        sum_pi *= dt * 4;
    }
    else
    {
        struct CalculateParameters calcStruct;
        MPI_Recv(&calcStruct, 1, mpiCalculateParametersDatatype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        PLOG_INFO << "Calculating integral, node: " << node;
        
        A = calcStruct.A;
        B = calcStruct.B;
        C = calcStruct.C;
        N_quadratic = calcStruct.N_quadratic;
        N_Pi = calcStruct.N_Pi;
        start_point = calcStruct.start_point;
        end_point = calcStruct.end_point;

        sum_quadratic = 0;
        dx = (end_point - start_point) / N_quadratic;
        dt = 1.0 / N_Pi;

        for (int i = node * (N_quadratic / numOfNodes); i <= node * N_quadratic / numOfNodes + N_quadratic / numOfNodes; i++)
        {
            sum_quadratic += calculate_quadratic_formula(start_point + i * dx);
        }
        sum_quadratic *= dx;

        PLOG_INFO << "Calculating Pi, node: " << node;

        for (int i = node * (N_Pi / numOfNodes); i <= node * N_Pi / numOfNodes + N_Pi / numOfNodes; i++)
        {
            sum_pi += calculate_pi_approx(i * dt);
        }
        sum_pi *= dt * 4;
    }
    PLOG_INFO << "Finished processing, node: " << node;

    MPI_Reduce(&sum_quadratic, &result_quadratic, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_pi, &result_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (node == 0)
    {
        endTime = MPI_Wtime();
        parallelTimeTaken = endTime - startTime;
        resultStream.str(std::string());
        resultStream << "Parallized result: " << static_cast<double>(result_quadratic);
        PLOG_INFO << resultStream.str();
        resultStream.str(std::string());
        resultStream << std::setprecision(piPrecision);
        resultStream << "Parallized Pi result: " << static_cast<double>(result_pi);
        PLOG_INFO << resultStream.str();
        PLOG_INFO << "Parallized time: " << parallelTimeTaken << " second(s)";
    }
    MPI_Type_free(&mpiCalculateParametersDatatype);
    MPI_Finalize();
    return 0;
}