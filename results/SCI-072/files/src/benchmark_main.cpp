// Main benchmark runner for MAPF experiments
#include "../include/mapf_types.h"
#include "../include/cbs.h"
#include "../include/eecbs.h"
#include "../include/lacam.h"
#include "../include/prioritized_planning.h"
#include "../include/lifelong_mapf.h"
#include "../include/distributed_mapf.h"
#include "../include/benchmark_generator.h"
#include <fstream>
#include <iomanip>
#include <string>

struct BenchmarkResult {
    std::string algorithm;
    std::string mapType;
    int numAgents;
    int mapSize;
    double runtime_ms;
    int totalCost;
    int makespan;
    bool solved;
};

void runScalabilityBenchmark(std::ofstream& csv) {
    std::cout << "=== Scalability Benchmark ===" << std::endl;
    std::vector<int> agentCounts = {5, 10, 20, 30, 50, 75, 100, 150, 200};
    int mapSize = 32;

    auto map = BenchmarkGenerator::generateMap(BenchmarkGenerator::RANDOM, mapSize, mapSize, 0.1, 42);

    for (int n : agentCounts) {
        auto inst = BenchmarkGenerator::generateInstance(map, n, 42);
        int actualAgents = inst.numAgents();
        std::cout << "  Agents: " << actualAgents << std::flush;

        // CBS (only for small instances)
        if (actualAgents <= 30) {
            auto sol = CBS::solve(inst, 30000);
            csv << "CBS,random," << actualAgents << "," << mapSize << ","
                << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                << sol.valid() << "\n";
            std::cout << " CBS:" << (sol.valid() ? std::to_string((int)sol.runtime_ms) + "ms" : "FAIL");
        } else {
            csv << "CBS,random," << actualAgents << "," << mapSize << ",-1,-1,-1,0\n";
        }

        // EECBS
        {
            auto sol = EECBS::solve(inst, 1.5, 30000);
            csv << "EECBS,random," << actualAgents << "," << mapSize << ","
                << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                << sol.valid() << "\n";
            std::cout << " EECBS:" << (sol.valid() ? std::to_string((int)sol.runtime_ms) + "ms" : "FAIL");
        }

        // Prioritized Planning
        {
            auto sol = PrioritizedPlanning::solve(inst, 30000);
            csv << "PP,random," << actualAgents << "," << mapSize << ","
                << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                << sol.valid() << "\n";
            std::cout << " PP:" << (sol.valid() ? std::to_string((int)sol.runtime_ms) + "ms" : "FAIL");
        }

        // LaCAM
        {
            auto sol = LaCAM::solve(inst, 30000);
            csv << "LaCAM,random," << actualAgents << "," << mapSize << ","
                << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                << sol.valid() << "\n";
            std::cout << " LaCAM:" << (sol.valid() ? std::to_string((int)sol.runtime_ms) + "ms" : "FAIL");
        }

        std::cout << std::endl;
    }
}

void runMapTypeBenchmark(std::ofstream& csv) {
    std::cout << "\n=== Map Type Benchmark ===" << std::endl;
    int mapSize = 32;
    int agents = 20;

    struct MapConfig {
        std::string name;
        BenchmarkGenerator::MapType type;
        double density;
    };
    std::vector<MapConfig> maps = {
        {"empty", BenchmarkGenerator::EMPTY, 0.0},
        {"random10", BenchmarkGenerator::RANDOM, 0.1},
        {"random20", BenchmarkGenerator::RANDOM, 0.2},
        {"warehouse", BenchmarkGenerator::WAREHOUSE, 0.0},
        {"maze", BenchmarkGenerator::MAZE, 0.0}
    };

    for (auto& mc : maps) {
        auto map = BenchmarkGenerator::generateMap(mc.type, mapSize, mapSize, mc.density, 42);
        auto inst = BenchmarkGenerator::generateInstance(map, agents, 42);
        int actual = inst.numAgents();
        std::cout << "  Map: " << mc.name << " agents: " << actual << std::flush;

        // CBS
        {
            auto sol = CBS::solve(inst, 30000);
            csv << "CBS," << mc.name << "," << actual << "," << mapSize << ","
                << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                << sol.valid() << "\n";
            std::cout << " CBS:" << (sol.valid() ? "OK" : "FAIL");
        }

        // EECBS
        {
            auto sol = EECBS::solve(inst, 1.5, 30000);
            csv << "EECBS," << mc.name << "," << actual << "," << mapSize << ","
                << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                << sol.valid() << "\n";
            std::cout << " EECBS:" << (sol.valid() ? "OK" : "FAIL");
        }

        // PP
        {
            auto sol = PrioritizedPlanning::solve(inst, 30000);
            csv << "PP," << mc.name << "," << actual << "," << mapSize << ","
                << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                << sol.valid() << "\n";
            std::cout << " PP:" << (sol.valid() ? "OK" : "FAIL");
        }

        // LaCAM
        {
            auto sol = LaCAM::solve(inst, 30000);
            csv << "LaCAM," << mc.name << "," << actual << "," << mapSize << ","
                << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                << sol.valid() << "\n";
            std::cout << " LaCAM:" << (sol.valid() ? "OK" : "FAIL");
        }
        std::cout << std::endl;
    }
}

void runSuboptimalityBenchmark(std::ofstream& csv2) {
    std::cout << "\n=== Suboptimality Analysis ===" << std::endl;
    int mapSize = 16;
    std::vector<double> wValues = {1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0};

    csv2 << "w,agents,runtime_ms,totalCost,optimal_cost,ratio\n";

    for (int agents : {5, 10, 15, 20}) {
        auto map = BenchmarkGenerator::generateMap(BenchmarkGenerator::RANDOM, mapSize, mapSize, 0.1, 42);
        auto inst = BenchmarkGenerator::generateInstance(map, agents, 42);
        int actual = inst.numAgents();

        // Get optimal cost with CBS
        auto optSol = CBS::solve(inst, 60000);
        int optCost = optSol.valid() ? optSol.totalCost() : -1;

        for (double w : wValues) {
            auto sol = EECBS::solve(inst, w, 30000);
            double ratio = (optCost > 0 && sol.valid()) ? (double)sol.totalCost() / optCost : -1;
            csv2 << w << "," << actual << "," << sol.runtime_ms << ","
                 << (sol.valid() ? sol.totalCost() : -1) << ","
                 << optCost << "," << ratio << "\n";
            std::cout << "  w=" << w << " agents=" << actual
                     << " cost=" << (sol.valid() ? sol.totalCost() : -1)
                     << " opt=" << optCost << " ratio=" << ratio << std::endl;
        }
    }
}

void runLifelongBenchmark(std::ofstream& csv3) {
    std::cout << "\n=== Lifelong MAPF Benchmark ===" << std::endl;
    csv3 << "agents,mapSize,tasks,completed,avgServiceTime,runtime_ms,replans\n";

    for (int agents : {5, 10, 20, 50, 100}) {
        int mapSize = 32;
        auto map = BenchmarkGenerator::generateMap(BenchmarkGenerator::WAREHOUSE, mapSize, mapSize, 0.0, 42);
        int tasks = agents * 10;
        auto m = LifelongMAPF::simulate(map, agents, tasks, 10, 20, 60000, 42);
        csv3 << agents << "," << mapSize << "," << tasks << ","
             << m.totalTasksCompleted << "," << m.avgServiceTime << ","
             << m.totalRuntime_ms << "," << m.replanCount << "\n";
        std::cout << "  Agents=" << agents << " completed=" << m.totalTasksCompleted
                 << "/" << tasks << " avgService=" << m.avgServiceTime << std::endl;
    }
}

void runDistributedBenchmark(std::ofstream& csv4) {
    std::cout << "\n=== Distributed MAPF Benchmark ===" << std::endl;
    csv4 << "agents,commRadius,dropRate,cost,makespan,conflicts,runtime_ms,avgMsgs\n";

    int mapSize = 32;
    auto map = BenchmarkGenerator::generateMap(BenchmarkGenerator::RANDOM, mapSize, mapSize, 0.1, 42);

    for (int agents : {10, 20, 50}) {
        auto inst = BenchmarkGenerator::generateInstance(map, agents, 42);
        int actual = inst.numAgents();

        for (double radius : {5.0, 10.0, 20.0, 50.0}) {
            for (double drop : {0.0, 0.1, 0.3}) {
                DistributedMAPF::CommConfig cfg = {radius, drop, 100};
                auto m = DistributedMAPF::solve(inst, cfg, 300, 42);
                csv4 << actual << "," << radius << "," << drop << ","
                     << m.totalCost << "," << m.makespan << "," << m.conflicts << ","
                     << m.runtime_ms << "," << m.avgMsgsPerStep << "\n";
                std::cout << "  a=" << actual << " r=" << radius << " d=" << drop
                         << " conflicts=" << m.conflicts << std::endl;
            }
        }
    }
}

void runWarehouseLargeScale(std::ofstream& csv5) {
    std::cout << "\n=== Large-Scale Warehouse Benchmark ===" << std::endl;
    csv5 << "agents,mapSize,algorithm,runtime_ms,totalCost,makespan,solved\n";

    for (int mapSize : {32, 64}) {
        auto map = BenchmarkGenerator::generateMap(BenchmarkGenerator::WAREHOUSE, mapSize, mapSize, 0.0, 42);
        for (int agents : {10, 50, 100, 200, 500}) {
            auto inst = BenchmarkGenerator::generateInstance(map, agents, 42);
            int actual = inst.numAgents();
            std::cout << "  Map " << mapSize << "x" << mapSize << " agents=" << actual << std::flush;

            // PP
            {
                auto sol = PrioritizedPlanning::solve(inst, 30000);
                csv5 << actual << "," << mapSize << ",PP,"
                     << sol.runtime_ms << "," << sol.totalCost() << "," << sol.makespan() << ","
                     << sol.valid() << "\n";
                std::cout << " PP:" << (sol.valid() ? std::to_string((int)sol.runtime_ms) + "ms" : "FAIL");
            }

            std::cout << std::endl;
        }
    }
}

int main() {
    std::cout << "MAPF Benchmark Suite" << std::endl;
    std::cout << "====================" << std::endl;

    {
        std::ofstream csv1("benchmarks/scalability.csv");
        csv1 << "algorithm,mapType,agents,mapSize,runtime_ms,totalCost,makespan,solved\n";
        runScalabilityBenchmark(csv1);
        runMapTypeBenchmark(csv1);
        csv1.close();
    }

    std::ofstream csv2("benchmarks/suboptimality.csv");
    runSuboptimalityBenchmark(csv2);
    csv2.close();

    std::ofstream csv3("benchmarks/lifelong.csv");
    runLifelongBenchmark(csv3);
    csv3.close();

    std::ofstream csv4("benchmarks/distributed.csv");
    runDistributedBenchmark(csv4);
    csv4.close();

    std::ofstream csv5("benchmarks/warehouse_large.csv");
    runWarehouseLargeScale(csv5);
    csv5.close();

    std::cout << "\nAll benchmarks complete. Results in benchmarks/" << std::endl;
    return 0;
}
