#pragma once
#include "mapf_types.h"

class ConstrainedAStar {
public:
    static Path search(const GridMap& map, const Location& start, const Location& goal,
                       const std::vector<Constraint>& constraints, int agent, int maxTime = 500) {
        struct Node {
            Location loc;
            int g, f, time;
            bool operator>(const Node& o) const { return f > o.f; }
        };

        auto isConstrained = [&](const Location& from, const Location& to, int t) -> bool {
            for (auto& c : constraints) {
                if (c.agent != agent) continue;
                if (c.type == Conflict::VERTEX && c.loc1 == to && c.time == t) return true;
                if (c.type == Conflict::EDGE && c.loc1 == from && c.loc2 == to && c.time == t) return true;
            }
            return false;
        };

        using PQ = std::priority_queue<Node, std::vector<Node>, std::greater<Node>>;
        PQ open;
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> best; // [t][y][x]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, Location>>> parent;

        int h0 = map.heuristic(start, goal);
        open.push({start, 0, h0, 0});
        best[0][start.y][start.x] = 0;

        while (!open.empty()) {
            auto cur = open.top(); open.pop();
            if (cur.loc == goal && cur.time >= (int)constraints.size()) {
                // check no future constraints on goal
                bool hasFuture = false;
                for (auto& c : constraints) {
                    if (c.agent == agent && c.loc1 == goal && c.time > cur.time) {
                        hasFuture = true; break;
                    }
                }
                if (!hasFuture) {
                    Path path;
                    Location l = cur.loc;
                    int t = cur.time;
                    std::vector<Location> locs;
                    locs.push_back(l);
                    while (t > 0) {
                        l = parent[t][l.y][l.x];
                        locs.push_back(l);
                        t--;
                    }
                    std::reverse(locs.begin(), locs.end());
                    path.locations = locs;
                    return path;
                }
            }

            if (cur.time >= maxTime) continue;

            for (auto& next : map.neighbors(cur.loc)) {
                if (isConstrained(cur.loc, next, cur.time + 1)) continue;
                int ng = cur.g + 1;
                int nt = cur.time + 1;
                auto& b = best[nt][next.y][next.x];
                if (b == 0 || ng < b) {
                    if (!(nt == 0 && next == start)) b = ng;
                    else if (b == 0) b = ng;
                    else continue;
                    parent[nt][next.y][next.x] = cur.loc;
                    open.push({next, ng, ng + map.heuristic(next, goal), nt});
                }
            }
        }
        return Path{}; // no solution
    }
};
