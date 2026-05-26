#pragma once
#include "mapf_types.h"
#include <stack>

class LaCAM {
public:
    static MAPFSolution solve(const MAPFInstance& inst, double timeLimitMs = 30000) {
        auto startTime = std::chrono::steady_clock::now();
        int n = inst.numAgents();
        if (n == 0) return {{}, 0};

        // Config = joint position of all agents
        using Config = std::vector<Location>;

        auto configHash = [&](const Config& c) -> size_t {
            size_t h = 0;
            for (auto& l : c) {
                h ^= std::hash<int>()(l.x * 10007 + l.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        };

        auto isGoal = [&](const Config& c) -> bool {
            for (int i = 0; i < n; i++)
                if (c[i] != inst.goals[i]) return false;
            return true;
        };

        struct TreeNode {
            Config config;
            int parent;
            int depth;
        };

        std::vector<TreeNode> tree;
        std::unordered_set<size_t> visited;

        Config startCfg = inst.starts;
        tree.push_back({startCfg, -1, 0});
        visited.insert(configHash(startCfg));

        // Use BFS-like search with greedy expansion
        std::queue<int> bfsQ;
        bfsQ.push(0);

        int maxDepth = inst.map.width * inst.map.height;
        std::mt19937 rng(42);

        while (!bfsQ.empty()) {
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - startTime).count();
            if (elapsed > timeLimitMs) break;

            int idx = bfsQ.front(); bfsQ.pop();

            if (isGoal(tree[idx].config)) {
                MAPFSolution sol;
                sol.paths.resize(n);
                std::vector<int> seq;
                int cur = idx;
                while (cur != -1) { seq.push_back(cur); cur = tree[cur].parent; }
                std::reverse(seq.begin(), seq.end());
                for (int i = 0; i < n; i++)
                    for (int ci : seq)
                        sol.paths[i].locations.push_back(tree[ci].config[i]);
                sol.runtime_ms = (double)std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::steady_clock::now() - startTime).count() / 1000.0;
                return sol;
            }

            if (tree[idx].depth >= maxDepth) continue;

            // Try multiple successor configs
            for (int attempt = 0; attempt < 3; attempt++) {
                Config next(n);
                std::unordered_set<int> occupied;
                bool valid = true;

                // Determine agent order
                std::vector<int> order(n);
                for (int i = 0; i < n; i++) order[i] = i;
                if (attempt > 0) std::shuffle(order.begin(), order.end(), rng);
                else {
                    std::sort(order.begin(), order.end(), [&](int a, int b) {
                        return inst.map.heuristic(tree[idx].config[a], inst.goals[a]) >
                               inst.map.heuristic(tree[idx].config[b], inst.goals[b]);
                    });
                }

                for (int i : order) {
                    Location cur_loc = tree[idx].config[i];
                    auto nbrs = inst.map.neighbors(cur_loc);

                    std::sort(nbrs.begin(), nbrs.end(), [&](const Location& a, const Location& b) {
                        return inst.map.heuristic(a, inst.goals[i]) < inst.map.heuristic(b, inst.goals[i]);
                    });

                    bool placed = false;
                    for (auto& nb : nbrs) {
                        int key = nb.y * inst.map.width + nb.x;
                        if (!occupied.count(key)) {
                            next[i] = nb;
                            occupied.insert(key);
                            placed = true;
                            break;
                        }
                    }
                    if (!placed) {
                        int key = cur_loc.y * inst.map.width + cur_loc.x;
                        if (!occupied.count(key)) {
                            next[i] = cur_loc;
                            occupied.insert(key);
                        } else { valid = false; break; }
                    }
                }

                if (!valid) continue;
                size_t hv = configHash(next);
                if (visited.count(hv)) continue;
                visited.insert(hv);

                int childIdx = (int)tree.size();
                tree.push_back({next, idx, tree[idx].depth + 1});
                bfsQ.push(childIdx);
            }
        }

        return {{}, 0};
    }
};
