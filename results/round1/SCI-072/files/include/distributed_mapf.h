#pragma once
#include "mapf_types.h"
#include <random>

// Distributed MAPF under communication constraints
class DistributedMAPF {
public:
    struct CommConfig {
        double commRadius;     // communication range
        double msgDropRate;    // probability of message loss
        int maxMsgPerStep;     // bandwidth limit
    };

    struct DistMetrics {
        int totalCost;
        int makespan;
        int conflicts;
        double runtime_ms;
        double avgMsgsPerStep;
    };

    static DistMetrics solve(const MAPFInstance& inst, const CommConfig& comm,
                             int maxSteps = 200, unsigned seed = 42) {
        auto startTime = std::chrono::steady_clock::now();
        std::mt19937 rng(seed);
        int n = inst.numAgents();

        // Each agent plans independently with local information
        std::vector<Location> positions = inst.starts;
        std::vector<std::vector<Location>> planned(n);

        // Initial individual A* plans
        for (int i = 0; i < n; i++) {
            auto path = ConstrainedAStar::search(inst.map, inst.starts[i], inst.goals[i], {}, i);
            if (!path.locations.empty()) planned[i] = path.locations;
            else planned[i] = {inst.starts[i]};
        }

        DistMetrics metrics = {0, 0, 0, 0.0, 0.0};
        int totalMsgs = 0;
        int totalConflicts = 0;

        for (int step = 0; step < maxSteps; step++) {
            bool allDone = true;
            for (int i = 0; i < n; i++) {
                if (positions[i] != inst.goals[i]) { allDone = false; break; }
            }
            if (allDone) break;

            // Communication phase: agents share positions with neighbors
            std::vector<std::vector<int>> neighbors(n);
            int msgsThisStep = 0;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double dist = std::sqrt(
                        std::pow(positions[i].x - positions[j].x, 2) +
                        std::pow(positions[i].y - positions[j].y, 2));
                    if (dist <= comm.commRadius) {
                        // Check message drop
                        std::uniform_real_distribution<> drop(0, 1);
                        if (drop(rng) > comm.msgDropRate && msgsThisStep < comm.maxMsgPerStep) {
                            neighbors[i].push_back(j);
                            neighbors[j].push_back(i);
                            msgsThisStep += 2;
                        }
                    }
                }
            }
            totalMsgs += msgsThisStep;

            // Local conflict resolution
            std::vector<Location> nextPos(n);
            for (int i = 0; i < n; i++) {
                int t = step + 1;
                if (t < (int)planned[i].size()) nextPos[i] = planned[i][t];
                else nextPos[i] = positions[i];
            }

            // Detect and resolve local conflicts
            for (int i = 0; i < n; i++) {
                for (int j : neighbors[i]) {
                    if (j <= i) continue;
                    if (nextPos[i] == nextPos[j]) {
                        totalConflicts++;
                        // Lower priority agent waits
                        int di = inst.map.heuristic(positions[i], inst.goals[i]);
                        int dj = inst.map.heuristic(positions[j], inst.goals[j]);
                        if (di < dj) nextPos[j] = positions[j];
                        else nextPos[i] = positions[i];
                    }
                }
            }

            // Move
            for (int i = 0; i < n; i++) positions[i] = nextPos[i];

            // Replan for agents that deviated
            for (int i = 0; i < n; i++) {
                int t = step + 1;
                if (t < (int)planned[i].size() && planned[i][t] != positions[i]) {
                    auto path = ConstrainedAStar::search(inst.map, positions[i], inst.goals[i], {}, i);
                    if (!path.locations.empty()) {
                        planned[i].clear();
                        planned[i] = path.locations;
                        // Reset plan indexing
                    }
                }
            }

            metrics.makespan = step + 1;
        }

        metrics.totalCost = 0;
        for (int i = 0; i < n; i++) {
            metrics.totalCost += inst.map.heuristic(inst.starts[i], positions[i]);
        }
        metrics.conflicts = totalConflicts;
        metrics.avgMsgsPerStep = (metrics.makespan > 0) ? (double)totalMsgs / metrics.makespan : 0;
        metrics.runtime_ms = (double)std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - startTime).count() / 1000.0;
        return metrics;
    }

private:
    // Simple A* without constraints for distributed use
    static Path simpleAStar(const GridMap& map, const Location& start, const Location& goal) {
        return ConstrainedAStar::search(map, start, goal, {}, 0);
    }
};
