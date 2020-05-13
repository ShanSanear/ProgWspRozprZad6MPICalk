#include <iostream>
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <stdio.h>
#include <ctime>
#include <mpi.h>
#include <plog/Log.h>
#include <plog/Appenders/ConsoleAppender.h>

double A, B, C;

struct CalculateParameters
{
    double A;
    double B;
    double C;
    int N;
    double start_point;
    double end_point;
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
    int N;
    const int struct_size = 6;
    double sum, dx, start_point, end_point, result;
    double startTime, endTime, parallelTimeTaken, timeSingle;
    plog::RollingFileAppender<plog::CsvFormatter> fileAppender("Datalogger.txt", 1048576, 5);
    plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;
    plog::init(plog::info, &fileAppender).addAppender(&consoleAppender);
    int node, numOfNodes;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    int lengths[struct_size] = {1, 1, 1, 1, 1, 1};
    MPI_Datatype types[struct_size] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint displacements[struct_size] = {
        offsetof(CalculateParameters, A),
        offsetof(CalculateParameters, B),
        offsetof(CalculateParameters, C),
        offsetof(CalculateParameters, N),
        offsetof(CalculateParameters, start_point),
        offsetof(CalculateParameters, end_point)
        };

    MPI_Datatype mpiCalculateParametersDatatype;
    MPI_Type_create_struct(struct_size, lengths, displacements, types, &mpiCalculateParametersDatatype);
    MPI_Type_commit(&mpiCalculateParametersDatatype);
    std::stringstream resultStream;
    resultStream << std::fixed << std::setprecision(6);
    if (node == 0)
    {
        PLOG_INFO << "Getting input values";
        A = get_double_from_stdin("Specify a:");
        B = get_double_from_stdin("Specify b:");
        C = get_double_from_stdin("Specify c:");
        N = get_int_from_stdin("Specifiy N:");
        start_point = get_double_from_stdin("Specify start point:");
        end_point = get_double_from_stdin("Specify end point:");

        sum = 0;
        dx = (end_point - start_point) / N;

        startTime = MPI_Wtime();
        PLOG_INFO << "Processing on single node";
        for (int i = 0; i <= N; i++)
        {
            sum += calculate_formula(start_point + i * dx);
        }

        sum *= dx;

        endTime = MPI_Wtime();

        timeSingle = endTime - startTime;
        resultStream << "Single node result: " << static_cast<double>(sum);
        PLOG_INFO << resultStream.str();
        PLOG_INFO << "Single node time: " << timeSingle << " second(s)";
        struct CalculateParameters calculateStruct;
        //There is no MPI_Barrier here because only one node can reach this section
        startTime = MPI_Wtime();
        calculateStruct.A = A;
        calculateStruct.B = B;
        calculateStruct.C = C;
        calculateStruct.N = N;
        calculateStruct.start_point = start_point;
        calculateStruct.end_point = end_point;
        PLOG_INFO << "Sending data to other nodes";
        for (int i = 1; i < numOfNodes; i++)
        {
            MPI_Send(&calculateStruct, 1, mpiCalculateParametersDatatype, i, 0, MPI_COMM_WORLD);
            
        }

        sum = 0;
        dx = (end_point - start_point) / N;
        PLOG_INFO << "Processing data on node 0";
        for (int i = 0; i <= N / numOfNodes; i++)
        {
            sum += calculate_formula(start_point + i * dx);
        }

        sum *= dx;
    }
    else
    {
        struct CalculateParameters calcStruct;
        MPI_Recv(&calcStruct, 1, mpiCalculateParametersDatatype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        PLOG_INFO << "Processing data, node: " << node;
        
        A = calcStruct.A;
        B = calcStruct.B;
        C = calcStruct.C;
        N = calcStruct.N;
        start_point = calcStruct.start_point;
        end_point = calcStruct.end_point;

        sum = 0;
        dx = (end_point - start_point) / N;

        for (int i = node * (N / numOfNodes); i <= node * N / numOfNodes + N / numOfNodes; i++)
        {
            sum += calculate_formula(start_point + i * dx);
        }
        sum *= dx;
    }
    PLOG_INFO << "Finished processing, node: " << node;

    MPI_Reduce(&sum, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (node == 0)
    {
        endTime = MPI_Wtime();
        parallelTimeTaken = endTime - startTime;
        resultStream.str(std::string());
        resultStream << "Parallized result: " << static_cast<double>(result);
        PLOG_INFO << resultStream.str();
        PLOG_INFO << "Parallized time: " << parallelTimeTaken << " second(s)";
    }
    MPI_Type_free(&mpiCalculateParametersDatatype);
    MPI_Finalize();
    return 0;
}