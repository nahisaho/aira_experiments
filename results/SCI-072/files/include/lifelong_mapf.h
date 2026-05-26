#pragma once
#include "mapf_types.h"
#include "prioritized_planning.h"
#include <random>

// Rolling-Horizon Collision Resolution for Lifelong MAPF
class LifelongMAPF {
public:
    struct Task {
        Location pickup, delivery;
    };

    struct Metrics {
        int totalTasksCompleted;
        double avgServiceTime;
        double totalRuntime_ms;
        int replanCount;
        std::vector<int> tasksCompletedPerWindow;
    };

    static Metrics simulate(const GridMap& map, int numAgents, int numTasks,
                            int windowSize = 5, int planHorizon = 20,
                            double timeLimitMs = 60000, unsigned seed = 42) {
        auto startTime = std::chrono::steady_clock::now();
        std::mt19937 rng(seed);
        Metrics metrics = {0, 0.0, 0.0, 0, {}};

        std::vector<Location> freeCells;
        for (int y = 0; y < map.height; y++)
            for (int x = 0; x < map.width; x++)
                if (!map.obstacles[y][x]) freeCells.push_back({x, y});

        if (freeCells.size() < (size_t)(numAgents * 2)) return metrics;

        auto randLoc = [&]() -> Location {
            return freeCells[rng() % freeCells.size()];
        };

        // Generate tasks
        std::vector<Task> tasks(numTasks);
        for (auto& t : tasks) {
            t.pickup = randLoc();
            t.delivery = randLoc();
            while (t.delivery == t.pickup) t.delivery = randLoc();
        }

        // Initialize agents at random free cells (unique)
        std::vector<Location> positions(numAgents);
        std::vector<Location> goals(numAgents);
        std::vector<int> currentTask(numAgents, -1);
        std::vector<bool> hasPickedUp(numAgents, false);
        int nextTask = 0;

        std::shuffle(freeCells.begin(), freeCells.end(), rng);
        for (int i = 0; i < numAgents; i++) {
            positions[i] = freeCells[i];
        }

        auto assignTask = [&](int agent) {
            if (nextTask < numTasks) {
                currentTask[agent] = nextTask++;
                hasPickedUp[agent] = false;
                goals[agent] = tasks[currentTask[agent]].pickup;
            } else {
                currentTask[agent] = -1;
                goals[agent] = positions[agent];
            }
        };

        for (int i = 0; i < numAgents; i++) assignTask(i);

        int timestep = 0;
        int maxTimesteps = numTasks * 100 / std::max(1, numAgents) + 200;

        while (metrics.totalTasksCompleted < numTasks && timestep < maxTimesteps) {
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - startTime).count();
            if (elapsed > timeLimitMs) break;

            // Build windowed MAPF instance
            MAPFInstance windowInst;
            windowInst.map = map;
            windowInst.starts = positions;
            windowInst.goals = goals;

            auto sol = PrioritizedPlanning::solve(windowInst, 3000);
            metrics.replanCount++;

            int windowCompleted = 0;
            int steps = windowSize;
            if (!sol.valid()) steps = 0;
            else {
                int minLen = 999999;
                for (auto& p : sol.paths) minLen = std::min(minLen, (int)p.locations.size());
                steps = std::min(steps, minLen - 1);
            }

            for (int t = 1; t <= steps; t++) {
                for (int i = 0; i < numAgents; i++) {
                    if (t < (int)sol.paths[i].locations.size())
                        positions[i] = sol.paths[i].locations[t];
                }
                timestep++;

                for (int i = 0; i < numAgents; i++) {
                    if (currentTask[i] < 0) continue;
                    if (!hasPickedUp[i] && positions[i] == tasks[currentTask[i]].pickup) {
                        hasPickedUp[i] = true;
                        goals[i] = tasks[currentTask[i]].delivery;
                    } else if (hasPickedUp[i] && positions[i] == tasks[currentTask[i]].delivery) {
                        metrics.totalTasksCompleted++;
                        windowCompleted++;
                        assignTask(i);
                    }
                }
            }
            metrics.tasksCompletedPerWindow.push_back(windowCompleted);

            if (steps <= 0) {
                // Try random restart
                for (int i = 0; i < numAgents; i++) {
                    if (currentTask[i] >= 0 && positions[i] == goals[i]) continue;
                }
                break;
            }
        }

        metrics.totalRuntime_ms = (double)std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - startTime).count() / 1000.0;
        metrics.avgServiceTime = (metrics.totalTasksCompleted > 0) ?
            (double)timestep / metrics.totalTasksCompleted : 0;
        return metrics;
    }
};
