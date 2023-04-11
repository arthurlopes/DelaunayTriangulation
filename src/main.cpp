#include <iostream>
#include <random>
#include <vector>
#include "delaunay.h"

using namespace delaunay;

template<typename XDistribution, typename YDistribution>
std::vector<float2> generateRandomPoints(int num_points, XDistribution& x_dist, YDistribution& y_dist) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<float2> points;
    points.reserve(num_points);

    for (int i = 0; i < num_points; ++i) {
        points.push_back({static_cast<float>(x_dist(gen)), static_cast<float>(y_dist(gen))});
    }

    return points;
}

void printPoints(const std::vector<float2>& points) {
    std::cout << "Points:\n";
    for (const auto& point : points) {
        std::cout << "(" << point.x << ", " << point.y << ")\n";
    }
}

void printTriangles(const std::vector<std::shared_ptr<Triangle>>& triangles) {
    std::cout << "Triangles:\n";
    int index = 0;
    for (const auto& triangle_ptr : triangles) {
        std::cout << "Triangle " << index << ":\n";
        std::cout << "  p1: (" << triangle_ptr->a->x << ", " << triangle_ptr->a->y << ")\n";
        std::cout << "  p2: (" << triangle_ptr->b->x << ", " << triangle_ptr->b->y << ")\n";
        std::cout << "  p3: (" << triangle_ptr->c->x << ", " << triangle_ptr->c->y << ")\n";
        std::cout << std::endl;
        ++index;
    }
}

int main() {
    int num_points = 10; // Set the desired number of random points

    // Configure the distributions
    std::uniform_real_distribution<> x_dist(0.0, 10.0);
    std::uniform_real_distribution<> y_dist(0.0, 10.0);

    // Alternatively, you can use different distributions
    // std::normal_distribution<> x_dist(5.0, 2.0);
    // std::normal_distribution<> y_dist(5.0, 2.0);

    std::vector<float2> points = generateRandomPoints(num_points, x_dist, y_dist);
    printPoints(points);

    std::vector<std::shared_ptr<Triangle>> triangles;

    RunDelaunayTriangulation(points, triangles);
    printTriangles(triangles);

    return 0;
}
