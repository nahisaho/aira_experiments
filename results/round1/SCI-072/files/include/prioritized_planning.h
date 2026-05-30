#pragma once
#include "mapf_types.h"
#include <queue>

class PrioritizedPlanning {
public:
    // Simple space-time A* that avoids occupied cells in a reservation table
    static Path stAStar(const GridMap& map, const Location& start, const Location& goal,
                        const std::unordered_map<int, std::unordered_set<int>>& reservedVertices,
                        int maxTime = 300) {
        struct Node {
            Location loc; int g, f, time;
            bool operator>(const Node& o) const { return f > o.f || (f == o.f && g > o.g); }
        };

        auto encode = [&](int x, int y, int t) -> long long {
            return (long long)t * map.width * map.height + (long long)y * map.width + x;
        };

        auto isReserved = [&](int x, int y, int t) -> bool {
            int key = y * map.width + x;
            auto it = reservedVertices.find(t);
            if (it != reservedVertices.end() && it->second.count(key)) return true;
            return false;
        };

        std::priority_queue<Node, std::vector<Node>, std::greater<Node>> open;
        std::unordered_set<long long> closed;
        std::unordered_map<long long, long long> parent;

        int h0 = map.heuristic(start, goal);
        open.push({start, 0, h0, 0});

        while (!open.empty()) {
            auto cur = open.top(); open.pop();
            long long curKey = encode(cur.loc.x, cur.loc.y, cur.time);
            if (closed.count(curKey)) continue;
            closed.insert(curKey);

            if (cur.loc == goal) {
                // Check no future reservations at goal
                bool safe = true;
                for (int ft = cur.time + 1; ft <= cur.time + 5; ft++) {
                    if (isReserved(goal.x, goal.y, ft)) { safe = false; break; }
                }
                if (safe) {
                    Path path;
                    long long k = curKey;
                    std::vector<Location> locs;
                    locs.push_back(cur.loc);
                    while (parent.count(k)) {
                        k = parent[k];
                        int t2 = (int)(k / (map.width * map.height));
                        int rem = (int)(k % (map.width * map.height));
                        locs.push_back({rem % map.width, rem / map.width});
                    }
                    std::reverse(locs.begin(), locs.end());
                    path.locations = locs;
                    return path;
                }
            }

            if (cur.time >= maxTime) continue;

            for (auto& next : map.neighbors(cur.loc)) {
                if (isReserved(next.x, next.y, cur.time + 1)) continue;
                long long nk = encode(next.x, next.y, cur.time + 1);
                if (closed.count(nk)) continue;
                parent[nk] = curKey;
                int ng = cur.g + 1;
                open.push({next, ng, ng + map.heuristic(next, goal), cur.time + 1});
            }
        }
        return Path{};
    }

    static MAPFSolution solve(const MAPFInstance& inst, double timeLimitMs = 30000) {
        auto startTime = std::chrono::steady_clock::now();
        int n = inst.numAgents();
        MAPFSolution sol;
        sol.paths.resize(n);

        // Reservation table: time -> set of (y*width+x)
        std::unordered_map<int, std::unordered_set<int>> reserved;

        for (int i = 0; i < n; i++) {
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - startTime).count();
            if (elapsed > timeLimitMs) return {{}, 0};

            sol.paths[i] = stAStar(inst.map, inst.starts[i], inst.goals[i], reserved);
            if (sol.paths[i].locations.empty()) return {{}, 0};

            // Reserve this path
            for (int t = 0; t < (int)sol.paths[i].locations.size(); t++) {
                auto& loc = sol.paths[i].locations[t];
                reserved[t].insert(loc.y * inst.map.width + loc.x);
            }
            // Reserve goal for extra timesteps
            auto& gloc = sol.paths[i].locations.back();
            int plen = (int)sol.paths[i].locations.size();
            for (int t = plen; t < plen + 30; t++) {
                reserved[t].insert(gloc.y * inst.map.width + gloc.x);
            }
        }

        sol.runtime_ms = (double)std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - startTime).count() / 1000.0;
        return sol;
    }
};
