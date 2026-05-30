#pragma once
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <chrono>
#include <random>
#include <cassert>
#include <iostream>
#include <sstream>
#include <memory>

struct Location {
    int x, y;
    bool operator==(const Location& o) const { return x == o.x && y == o.y; }
    bool operator!=(const Location& o) const { return !(*this == o); }
};

struct LocationHash {
    size_t operator()(const Location& l) const {
        return std::hash<int>()(l.x) ^ (std::hash<int>()(l.y) << 16);
    }
};

struct TimedLocation {
    Location loc;
    int time;
};

struct Conflict {
    int agent1, agent2;
    Location loc1, loc2;
    int time;
    enum Type { VERTEX, EDGE } type;
};

struct Constraint {
    int agent;
    Location loc1, loc2;
    int time;
    Conflict::Type type;
};

struct Path {
    std::vector<Location> locations;
    int cost() const { return locations.empty() ? 0 : (int)locations.size() - 1; }
};

struct GridMap {
    int width, height;
    std::vector<std::vector<bool>> obstacles;

    GridMap() : width(0), height(0) {}
    GridMap(int w, int h) : width(w), height(h), obstacles(h, std::vector<bool>(w, false)) {}

    bool isValid(int x, int y) const {
        return x >= 0 && x < width && y >= 0 && y < height && !obstacles[y][x];
    }

    std::vector<Location> neighbors(const Location& l) const {
        std::vector<Location> result;
        static const int dx[] = {0, 1, 0, -1, 0};
        static const int dy[] = {1, 0, -1, 0, 0};
        for (int i = 0; i < 5; i++) {
            int nx = l.x + dx[i], ny = l.y + dy[i];
            if (isValid(nx, ny)) result.push_back({nx, ny});
        }
        return result;
    }

    int heuristic(const Location& a, const Location& b) const {
        return std::abs(a.x - b.x) + std::abs(a.y - b.y);
    }
};

struct MAPFInstance {
    GridMap map;
    std::vector<Location> starts;
    std::vector<Location> goals;
    int numAgents() const { return (int)starts.size(); }
};

struct MAPFSolution {
    std::vector<Path> paths;
    double runtime_ms;
    int totalCost() const {
        int c = 0;
        for (auto& p : paths) c += p.cost();
        return c;
    }
    int makespan() const {
        int m = 0;
        for (auto& p : paths) m = std::max(m, p.cost());
        return m;
    }
    bool valid() const { return !paths.empty(); }
};
