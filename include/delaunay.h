#include <string>
#include <vector>
#include <memory>
#include <list>
#include <stack>

namespace delaunay
{

typedef struct int2
{
    int x;
    int y;

    int2(int X, int Y)
        : x(X)
        , y(Y)
    {}
    int2()
        : x(0)
        , y(0)
    {}
} int2;

typedef struct float2
{
    float x;
    float y;
    int id;

    float2(float X, float Y, int id)
        : x(X)
        , y(Y)
        , id(id)
    {}
    float2(float X, float Y)
        : x(X)
        , y(Y)
        , id(0)
    {}
    float2()
        : x(0)
        , y(0)
        , id(0)
    {}
    float2(int2 a)
        : x(float(a.x))
        , y(float(a.y))
        , id(0)
    {}
    bool operator<(const float2& other) const {
        return x < other.x || (x == other.x && y < other.y);
    }
    bool operator==(const float2& other) const {
        return (x == other.x) && (y == other.y);
    }
    bool operator!=(const float2& other) const {
        return !(*this == other);
    }
} float2;
typedef struct double2
{
    double x;
    double y;

    double2(double X, double Y)
        : x(X)
        , y(Y)
    {}
    double2()
        : x(0)
        , y(0)
    {}
    double2(float2 a)
        : x(double(a.x))
        , y(double(a.y))
    {}
    double2(std::shared_ptr<float2> a)
        : x(double(a->x))
        , y(double(a->y))
    {}
    bool operator==(const double2& other) const {
        return (x == other.x) && (y == other.y);
    }
    bool operator!=(const double2& other) const {
        return !(*this == other);
    }
} double2;

typedef struct Triangle {
    // Vertex
    std::shared_ptr<float2> a;
    std::shared_ptr<float2> b;
    std::shared_ptr<float2> c;

    // Neighbor triangles
    std::shared_ptr<struct Triangle> t_a;
    std::shared_ptr<struct Triangle> t_b;
    std::shared_ptr<struct Triangle> t_c;

    // Neighbor edge (convex hull)
    std::shared_ptr<struct Edge> e_a;
    std::shared_ptr<struct Edge> e_b;
    std::shared_ptr<struct Edge> e_c;

    Triangle(float2 X, float2 Y, float2 Z)
        : a(std::make_shared<float2>(float2(X.x, X.y)))
        , b(std::make_shared<float2>(float2(Y.x, Y.y)))
        , c(std::make_shared<float2>(float2(Z.x, Z.y)))
        , t_a()
        , t_b()
        , t_c()
        , e_a()
        , e_b()
        , e_c()
    {}
    Triangle(std::shared_ptr<float2> X, std::shared_ptr<float2> Y, std::shared_ptr<float2> Z)
        : a(X)
        , b(Y)
        , c(Z)
        , t_a()
        , t_b()
        , t_c()
        , e_a()
        , e_b()
        , e_c()
    {}
    Triangle()
        : a()
        , b()
        , c()
        , t_a()
        , t_b()
        , t_c()
        , e_a()
        , e_b()
        , e_c()
    {}
    bool operator==(const Triangle& other) const {
        return (
            ((*a == *other.a) && (*b == *other.b) && (*c == *other.c)) ||
            ((*a == *other.b) && (*b == *other.c) && (*c == *other.a)) ||
            ((*a == *other.c) && (*b == *other.a) && (*c == *other.b))
            );
    }
    bool operator!=(const Triangle& other) const {
        return !(*this == other);
    }
} Triangle;

// Convex hull edge
typedef struct Edge {
    std::shared_ptr<float2> start, end;

    // triangle on the left side of the edge
    std::shared_ptr<Triangle> left;

    Edge(float2 X, float2 Y)
        : start(std::make_shared<float2>(float2(X.x, X.y)))
        , end(std::make_shared<float2>(float2(Y.x, Y.y)))
        , left()
    {}
    Edge(std::shared_ptr<float2> X, std::shared_ptr<float2> Y)
        : start(X)
        , end(Y)
        , left()
    {}
    Edge()
        : start()
        , end()
        , left()
    {}
    bool operator==(const Edge& other) const {
        return (start->x == other.start->x) &&
            (start->y == other.start->y) &&
            (end->x == other.end->x) &&
            (end->y == other.end->y);
    }
    bool operator!=(const Edge& other) const {
        return !(*this == other);
    }
} Edge;


void GlobalInit();
void GlobalTeardown();

double2 FindCircumcentre(std::shared_ptr<Triangle> const& t);
double Distance(double2 a, double2 b);
bool IsLegal(std::shared_ptr<Triangle> t1, std::shared_ptr<Triangle> t2);
void FlipEdge(std::shared_ptr<Triangle> t1, std::shared_ptr<Triangle> t2);
void Legalize(std::stack<std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>>>& stack);
void Legalize_recursive(std::shared_ptr<Triangle> t1, std::shared_ptr<Triangle> t2);
bool IsPointOnLeftOfEdge(Edge const& edge, std::shared_ptr<float2> const& point);
void InitializeDelaunay(std::vector<std::shared_ptr<float2>>& stars, std::list< std::shared_ptr<Edge>>& hull, std::vector<std::shared_ptr<Triangle>>& triangles);
void DelaunayTriangulation(std::vector<std::shared_ptr<float2>> const& stars, std::list< std::shared_ptr<Edge>>& hull, std::vector<std::shared_ptr<Triangle>>& triangles);

}