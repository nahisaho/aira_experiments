#pragma once
#include "mapf_types.h"
#include "astar.h"

class EECBS {
public:
    static MAPFSolution solve(const MAPFInstance& inst, double w = 1.5, double timeLimitMs = 30000) {
        auto startTime = std::chrono::steady_clock::now();
        int n = inst.numAgents();

        struct CTNode {
            std::vector<Path> paths;
            std::vector<Constraint> constraints;
            int cost;
            int lowerBound;
            bool operator>(const CTNode& o) const {
                if (lowerBound != o.lowerBound) return lowerBound > o.lowerBound;
                return cost > o.cost;
            }
        };

        CTNode root;
        root.constraints = {};
        root.paths.resize(n);
        root.lowerBound = 0;
        for (int i = 0; i < n; i++) {
            root.paths[i] = ConstrainedAStar::search(inst.map, inst.starts[i], inst.goals[i], {}, i);
            if (root.paths[i].locations.empty()) return {{}, 0};
            root.lowerBound += inst.map.heuristic(inst.starts[i], inst.goals[i]);
        }
        root.cost = 0;
        for (auto& p : root.paths) root.cost += p.cost();

        using PQ = std::priority_queue<CTNode, std::vector<CTNode>, std::greater<CTNode>>;
        PQ open;
        PQ focal;
        open.push(root);
        focal.push(root);

        int bestLB = root.lowerBound;
        int expanded = 0;

        while (!focal.empty()) {
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - startTime).count();
            if (elapsed > timeLimitMs) break;

            auto node = focal.top(); focal.pop();
            expanded++;

            if (node.cost <= (int)(w * bestLB)) {
                auto conflict = findConflict(node.paths);
                if (conflict.agent1 == -1) {
                    MAPFSolution sol;
                    sol.paths = node.paths;
                    sol.runtime_ms = (double)std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::steady_clock::now() - startTime).count() / 1000.0;
                    return sol;
                }

                for (int i = 0; i < 2; i++) {
                    CTNode child = node;
                    Constraint c;
                    c.agent = (i == 0) ? conflict.agent1 : conflict.agent2;
                    c.time = conflict.time;
                    c.type = conflict.type;
                    c.loc1 = conflict.loc1;
                    c.loc2 = conflict.loc2;
                    if (conflict.type == Conflict::EDGE && i == 1) std::swap(c.loc1, c.loc2);
                    child.constraints.push_back(c);
                    child.paths[c.agent] = ConstrainedAStar::search(
                        inst.map, inst.starts[c.agent], inst.goals[c.agent],
                        child.constraints, c.agent);
                    if (child.paths[c.agent].locations.empty()) continue;
                    child.cost = 0;
                    for (auto& p : child.paths) child.cost += p.cost();
                    child.lowerBound = node.lowerBound;
                    open.push(child);
                    if (child.cost <= (int)(w * bestLB)) {
                        focal.push(child);
                    }
                }
            }

            // Update focal list
            if (!open.empty()) {
                int newLB = open.top().lowerBound;
                if (newLB > bestLB) {
                    bestLB = newLB;
                    PQ newFocal;
                    while (!open.empty()) {
                        auto n2 = open.top(); open.pop();
                        if (n2.cost <= (int)(w * bestLB)) newFocal.push(n2);
                        else { open.push(n2); break; }
                    }
                    focal = newFocal;
                }
            }
        }
        return {{}, 0};
    }

private:
    static Conflict findConflict(const std::vector<Path>& paths) {
        int n = (int)paths.size();
        int maxT = 0;
        for (auto& p : paths) maxT = std::max(maxT, (int)p.locations.size());
        for (int t = 0; t < maxT; t++) {
            for (int i = 0; i < n; i++) {
                Location li = (t < (int)paths[i].locations.size()) ? paths[i].locations[t] : paths[i].locations.back();
                for (int j = i + 1; j < n; j++) {
                    Location lj = (t < (int)paths[j].locations.size()) ? paths[j].locations[t] : paths[j].locations.back();
                    if (li == lj) return {i, j, li, li, t, Conflict::VERTEX};
                }
            }
        }
        return {-1, -1, {0,0}, {0,0}, 0, Conflict::VERTEX};
    }
};
