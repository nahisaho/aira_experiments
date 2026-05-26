#pragma once
#include "mapf_types.h"
#include <random>
#include <fstream>

class BenchmarkGenerator {
public:
    enum MapType { EMPTY, RANDOM, WAREHOUSE, MAZE };

    static GridMap generateMap(MapType type, int width, int height,
                               double obstacleDensity = 0.2, unsigned seed = 42) {
        std::mt19937 rng(seed);
        GridMap map(width, height);

        switch (type) {
        case EMPTY: break;
        case RANDOM:
            for (int y = 0; y < height; y++)
                for (int x = 0; x < width; x++)
                    if (std::uniform_real_distribution<>(0, 1)(rng) < obstacleDensity)
                        map.obstacles[y][x] = true;
            break;
        case WAREHOUSE:
            generateWarehouse(map, rng);
            break;
        case MAZE:
            generateMaze(map, rng);
            break;
        }
        return map;
    }

    static MAPFInstance generateInstance(const GridMap& map, int numAgents, unsigned seed = 42) {
        std::mt19937 rng(seed);
        MAPFInstance inst;
        inst.map = map;

        std::vector<Location> freeCells;
        for (int y = 0; y < map.height; y++)
            for (int x = 0; x < map.width; x++)
                if (!map.obstacles[y][x]) freeCells.push_back({x, y});

        std::shuffle(freeCells.begin(), freeCells.end(), rng);

        int n = std::min(numAgents, (int)freeCells.size() / 2);
        for (int i = 0; i < n; i++) {
            inst.starts.push_back(freeCells[i]);
            inst.goals.push_back(freeCells[n + i]);
        }
        return inst;
    }

    static void saveInstance(const MAPFInstance& inst, const std::string& filename) {
        std::ofstream f(filename);
        f << inst.map.width << " " << inst.map.height << "\n";
        for (int y = 0; y < inst.map.height; y++) {
            for (int x = 0; x < inst.map.width; x++)
                f << (inst.map.obstacles[y][x] ? '@' : '.');
            f << "\n";
        }
        f << inst.numAgents() << "\n";
        for (int i = 0; i < inst.numAgents(); i++)
            f << inst.starts[i].x << " " << inst.starts[i].y << " "
              << inst.goals[i].x << " " << inst.goals[i].y << "\n";
    }

private:
    static void generateWarehouse(GridMap& map, std::mt19937& rng) {
        // Create aisle-shelf pattern
        int shelfWidth = 2;
        int aisleWidth = 2;
        int blockWidth = shelfWidth + aisleWidth;

        for (int y = 2; y < map.height - 2; y++) {
            for (int x = 2; x < map.width - 2; x++) {
                int bx = (x - 2) % blockWidth;
                int by = y % (shelfWidth + aisleWidth);
                if (bx < shelfWidth && by < shelfWidth) {
                    // Cross-aisles every 8 rows
                    if (y % 8 >= 2) map.obstacles[y][x] = true;
                }
            }
        }
    }

    static void generateMaze(GridMap& map, std::mt19937& rng) {
        for (int y = 0; y < map.height; y++)
            for (int x = 0; x < map.width; x++)
                map.obstacles[y][x] = true;

        // DFS maze generation
        std::vector<std::vector<bool>> visited(map.height, std::vector<bool>(map.width, false));
        std::stack<Location> stack;
        Location start = {1, 1};
        map.obstacles[1][1] = false;
        visited[1][1] = true;
        stack.push(start);

        int dx[] = {0, 2, 0, -2};
        int dy[] = {2, 0, -2, 0};

        while (!stack.empty()) {
            Location cur = stack.top();
            std::vector<int> dirs;
            for (int d = 0; d < 4; d++) {
                int nx = cur.x + dx[d], ny = cur.y + dy[d];
                if (nx > 0 && nx < map.width - 1 && ny > 0 && ny < map.height - 1 && !visited[ny][nx])
                    dirs.push_back(d);
            }
            if (dirs.empty()) { stack.pop(); continue; }
            int d = dirs[rng() % dirs.size()];
            int mx = cur.x + dx[d] / 2, my = cur.y + dy[d] / 2;
            int nx = cur.x + dx[d], ny = cur.y + dy[d];
            map.obstacles[my][mx] = false;
            map.obstacles[ny][nx] = false;
            visited[ny][nx] = true;
            stack.push({nx, ny});
        }
    }
};
