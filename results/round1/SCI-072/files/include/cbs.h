#pragma once
#include "mapf_types.h"
#include "astar.h"

class CBS {
public:
    struct CTNode {
        std::vector<Path> paths;
        std::vector<Constraint> constraints;
        int cost;
        bool operator>(const CTNode& o) const { return cost > o.cost; }
    };

    static MAPFSolution solve(const MAPFInstance& inst, double timeLimitMs = 30000) {
        auto startTime = std::chrono::steady_clock::now();
        int n = inst.numAgents();

        CTNode root;
        root.constraints = {};
        root.paths.resize(n);
        for (int i = 0; i < n; i++) {
            root.paths[i] = ConstrainedAStar::search(inst.map, inst.starts[i], inst.goals[i], {}, i);
            if (!root.paths[i].locations.size()) return {{}, 0};
        }
        root.cost = 0;
        for (auto& p : root.paths) root.cost += p.cost();

        using PQ = std::priority_queue<CTNode, std::vector<CTNode>, std::greater<CTNode>>;
        PQ open;
        open.push(root);
        int expanded = 0;

        while (!open.empty()) {
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - startTime).count();
            if (elapsed > timeLimitMs) break;

            auto node = open.top(); open.pop();
            expanded++;

            auto conflict = findFirstConflict(node.paths);
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
                if (conflict.type == Conflict::EDGE && i == 1) {
                    std::swap(c.loc1, c.loc2);
                }
                child.constraints.push_back(c);
                child.paths[c.agent] = ConstrainedAStar::search(
                    inst.map, inst.starts[c.agent], inst.goals[c.agent],
                    child.constraints, c.agent);
                if (child.paths[c.agent].locations.empty()) continue;
                child.cost = 0;
                for (auto& p : child.paths) child.cost += p.cost();
                open.push(child);
            }
        }
        return {{}, 0};
    }

private:
    static Conflict findFirstConflict(const std::vector<Path>& paths) {
        int n = (int)paths.size();
        int maxT = 0;
        for (auto& p : paths) maxT = std::max(maxT, (int)p.locations.size());

        for (int t = 0; t < maxT; t++) {
            for (int i = 0; i < n; i++) {
                Location li = getLocation(paths[i], t);
                for (int j = i + 1; j < n; j++) {
                    Location lj = getLocation(paths[j], t);
                    if (li == lj) return {i, j, li, li, t, Conflict::VERTEX};
                    if (t > 0) {
                        Location pi = getLocation(paths[i], t - 1);
                        Location pj = getLocation(paths[j], t - 1);
                        if (li == pj && lj == pi) return {i, j, pi, li, t, Conflict::EDGE};
                    }
                }
            }
        }
        return {-1, -1, {0,0}, {0,0}, 0, Conflict::VERTEX};
    }

    static Location getLocation(const Path& p, int t) {
        if (t < (int)p.locations.size()) return p.locations[t];
        return p.locations.back();
    }
};
