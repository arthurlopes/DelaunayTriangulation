#include "delaunay.h"
#include <stdexcept>
#include <list>
#include <iostream>
#include <stack>
#include <queue>
#include <unordered_set>
#include <cstdlib>
#include <fstream>

namespace delaunay {


double2 midpoint(const double2& p1, const double2& p2) {
    return { (p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0 };
}

double length(const double2& p1, const double2& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

double slope(const double2& p1, const double2& p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    return dx == 0.0 ? INFINITY : dy / dx;
}

void bisector(const double2& p1, const double2& p2, double& m, double& b) {
    double2 mid = midpoint(p1, p2);
    m = -1.0 / slope(p1, p2);
    b = mid.y - m * mid.x;
}

double2 Circumcentre(std::shared_ptr<Triangle> const& t) {
    double m1, b1, m2, b2;
    bisector(t->a, t->b, m1, b1);
    bisector(t->b, t->c, m2, b2);
    double cx = (b2 - b1) / (m1 - m2);
    double cy = m1 * cx + b1;

    return double2(cx, cy);
}

double2 FindCircumcentre_old(std::shared_ptr<Triangle> const& t) {
    double2 A = double2(t->a);
    double2 B = double2(t->b);
    double2 C = double2(t->c);

    double D = 1 / (2 * ((A.x * (B.y - C.y)) + B.x * (C.y - A.y) + C.x * (A.y - B.y)));
    double a = (A.x * A.x + A.y * A.y);
    double b = (B.x * B.x + B.y * B.y);
    double c = (C.x * C.x + C.y * C.y);
    double x = D * (a * (B.y - C.y) + b * (C.y - A.y) + c * (A.y - B.y));
    double y = D * (a * (C.x - B.x) + b * (A.x - C.x) + c * (B.x - A.x));

    return double2(x, y);
}

double Distance_old(double2 a, double2 b) {
    double dist = std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
    return dist;
}

double Distance(double2 a, double2 b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
}

bool CheckIfTrianglesShareEdge(std::shared_ptr<Triangle> t1, std::shared_ptr<Triangle> t2) {
    int count = 0;
    if (t1->a != t2->a && t1->a != t2->b && t1->a != t2->c) count++;
    if (t1->b != t2->a && t1->b != t2->b && t1->b != t2->c) count++;
    if (t1->c != t2->a && t1->c != t2->b && t1->c != t2->c) count++;

    if (count != 1) {
        return false;
    }

    return true;
}

bool IsLegal_old(std::shared_ptr<Triangle> t1, std::shared_ptr<Triangle> t2) {
    if (!CheckIfTrianglesShareEdge(t1, t2)) return true;

    double2 d;
    if (t1->a != t2->a && t1->a != t2->b && t1->a != t2->c) {
        d = *(t1->a);
    }
    else if (t1->b != t2->a && t1->b != t2->b && t1->b != t2->c) {
        d = *(t1->b);
    }
    else if (t1->c != t2->a && t1->c != t2->b && t1->c != t2->c) {
        d = *(t1->c);
    }

    double2 cc = Circumcentre(t2);
    double radius = Distance(cc, *(t2->a));
    double dist = Distance(cc, d);

    return dist >= radius;
}

bool IsLegal(std::shared_ptr<Triangle> t1, std::shared_ptr<Triangle> t2) {
    if (!CheckIfTrianglesShareEdge(t1, t2)) return true;

    double2 d;
    if (t1->a != t2->a && t1->a != t2->b && t1->a != t2->c) {
        d = *(t1->a);
    }
    else if (t1->b != t2->a && t1->b != t2->b && t1->b != t2->c) {
        d = *(t1->b);
    }
    else if (t1->c != t2->a && t1->c != t2->b && t1->c != t2->c) {
        d = *(t1->c);
    }

    double dx, dy, ex, ey, fx, fy, ap, bp, cp;
    dx = t2->a->x - d.x;
    dy = t2->a->y - d.y;
    ex = t2->b->x - d.x;
    ey = t2->b->y - d.y;
    fx = t2->c->x - d.x;
    fy = t2->c->y - d.y;

    ap = dx * dx + dy * dy;
    bp = ex * ex + ey * ey;
    cp = fx * fx + fy * fy;

    return (dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx)) < 0.0;
}

void FlipEdge(std::shared_ptr<Triangle> t1, std::shared_ptr<Triangle> t2) {
    std::shared_ptr<float2> a, b, c, d;
    std::shared_ptr<Triangle> t_a, t_b, t_c, t_d;
    std::shared_ptr<Edge> e_a, e_b, e_c, e_d;

    if (!CheckIfTrianglesShareEdge(t1, t2)) return;

    if (t1->a != t2->a && t1->a != t2->b && t1->a != t2->c) {
        b = t1->a;
        a = t1->b;
        t_b = t1->t_a;
        t_d = t1->t_c;
        e_b = t1->e_a;
        e_d = t1->e_c;
    }
    else if (t1->b != t2->a && t1->b != t2->b && t1->b != t2->c) {
        b = t1->b;
        a = t1->c;
        t_b = t1->t_b;
        t_d = t1->t_a;
        e_b = t1->e_b;
        e_d = t1->e_a;
    }
    else if (t1->c != t2->a && t1->c != t2->b && t1->c != t2->c) {
        b = t1->c;
        a = t1->a;
        t_b = t1->t_c;
        t_d = t1->t_b;
        e_b = t1->e_c;
        e_d = t1->e_b;
    }

    if (t2->a != t1->a && t2->a != t1->b && t2->a != t1->c) {
        c = t2->a;
        d = t2->b;
        t_c = t2->t_a;
        t_a = t2->t_c;
        e_c = t2->e_a;
        e_a = t2->e_c;
    }
    else if (t2->b != t1->a && t2->b != t1->b && t2->b != t1->c) {
        c = t2->b;
        d = t2->c;
        t_c = t2->t_b;
        t_a = t2->t_a;
        e_c = t2->e_b;
        e_a = t2->e_a;
    }
    else if (t2->c != t1->a && t2->c != t1->b && t2->c != t1->c) {
        c = t2->c;
        d = t2->a;
        t_c = t2->t_c;
        t_a = t2->t_b;
        e_c = t2->e_c;
        e_a = t2->e_b;
    }

    t1->a = b;
    t1->b = a;
    t1->c = c;
    t1->t_a = t_b;
    t1->t_b = t_a;
    t1->t_c = t2;
    t1->e_a = e_b;
    t1->e_b = e_a;
    t1->e_c = nullptr;

    if (t_a) {
        if (t_a->a == c) {
            t_a->t_a = t1;
        }
        else if (t_a->b == c) {
            t_a->t_b = t1;
        }
        else if (t_a->c == c) {
            t_a->t_c = t1;
        }
    }
    else {
        e_a->left = t1;
    }
    if (t_b) {
        if (t_b->a == a) {
            t_b->t_a = t1;
        }
        else if (t_b->b == a) {
            t_b->t_b = t1;
        }
        else if (t_b->c == a) {
            t_b->t_c = t1;
        }
    } else {
        e_b->left = t1;
    }

    t2->a = c;
    t2->b = d;
    t2->c = b;
    t2->t_a = t_c;
    t2->t_b = t_d;
    t2->t_c = t1;
    t2->e_a = e_c;
    t2->e_b = e_d;
    t2->e_c = nullptr;
    if (t_c) {
        if (t_c->a == d) {
            t_c->t_a = t2;
        }
        else if (t_c->b == d) {
            t_c->t_b = t2;
        }
        else if (t_c->c == d) {
            t_c->t_c = t2;
        }
    }
    else {
        e_c->left = t2;
    }
    if (t_d) {
        if (t_d->a == b) {
            t_d->t_a = t2;
        }
        else if (t_d->b == b) {
            t_d->t_b = t2;
        }
        else if (t_d->c == b) {
            t_d->t_c = t2;
        }
    }
    else {
        e_d->left = t2;
    }
}

void Legalize(std::stack<std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>>>& stack) {
    int count = 0;
    while (!stack.empty()) {
        count++;
        std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>> ts = stack.top();
        std::shared_ptr<Triangle> t1 = ts.first;
        std::shared_ptr<Triangle> t2 = ts.second;
        stack.pop();

        if (t1 && t2 && (!IsLegal(t1, t2))) {
            FlipEdge(t1, t2);

            if (t1->t_a && t1->t_a != t2 && !IsLegal(t1, t1->t_a))
                stack.push(std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>>(t1, t1->t_a));
            if (t1->t_b && t1->t_b != t2 && !IsLegal(t1, t1->t_b))
                stack.push(std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>>(t1, t1->t_b));
            if (t1->t_c && t1->t_c != t2 && !IsLegal(t1, t1->t_c))
                stack.push(std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>>(t1, t1->t_c));

            if (t2->t_a && t2->t_a != t1 && !IsLegal(t2, t2->t_a))
                stack.push(std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>>(t2, t2->t_a));
            if (t2->t_b && t2->t_b != t1 && !IsLegal(t2, t2->t_b))
                stack.push(std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>>(t2, t2->t_b));
            if (t2->t_c && t2->t_c != t1 && !IsLegal(t2, t2->t_c))
                stack.push(std::pair<std::shared_ptr<Triangle>, std::shared_ptr<Triangle>>(t2, t2->t_c));
        }
    }
}

void Legalize_recursive(std::shared_ptr<Triangle> tp1, std::shared_ptr<Triangle> tp2) {
    if ( tp1 && tp2 && (!IsLegal(tp1, tp2)) ) {
        FlipEdge(tp1, tp2);
        Triangle t1, t2;
        t1 = *tp1;
        t2 = *tp2;

        if (t1.t_a && t1.t_a != tp2)
            Legalize_recursive(tp1, t1.t_a);
        if (t1.t_b && t1.t_b != tp2)
            Legalize_recursive(tp1, t1.t_b);
        if (t1.t_c && t1.t_c != tp2)
            Legalize_recursive(tp1, t1.t_c);

        if (t2.t_a && t2.t_a != tp1)
            Legalize_recursive(tp2, t2.t_a);
        if (t2.t_b && t2.t_b != tp1)
            Legalize_recursive(tp2, t2.t_b);
        if (t2.t_c && t2.t_c != tp1)
            Legalize_recursive(tp2, t2.t_c);
    }
}

bool IsPointOnLeftOfEdge(Edge const& edge, std::shared_ptr<float2> const& point) {
    double v, px, py, edge_start_x, edge_start_y, edge_end_x, edge_end_y;
    px = point->x;
    py = point->y;
    edge_start_x = edge.start->x;
    edge_start_y = edge.start->y;
    edge_end_x = edge.end->x;
    edge_end_y = edge.end->y;

    v = (px - edge_start_x) * (edge_end_y - edge_start_y) - (edge_end_x - edge_start_x) * (py - edge_start_y);

    return v <= 0.f;
}

void InitializeDelaunay(std::vector<std::shared_ptr<float2>>& stars, std::list<std::shared_ptr<Edge>>& hull, std::vector<std::shared_ptr<Triangle>>& triangles) {
    std::sort(stars.begin(), stars.end(), [](const std::shared_ptr<float2>& a, const std::shared_ptr<float2>& b) {
        //return a->x < b->x;
        return a->x < b->x || (a->x == b->x && a->y < b->y);
    });

    for (int i = 0; i < stars.size(); i++)
        stars[i]->id = i;

    std::shared_ptr<float2> p0, p1, p2;
    p0 = stars[0];
    p1 = stars[1];
    p2 = stars[2];
    std::shared_ptr<Edge> p0p1 = std::make_shared<Edge>(p0, p1);
    if (IsPointOnLeftOfEdge(*p0p1, stars[2])) {
        std::shared_ptr<Triangle> t_p = std::make_shared<Triangle>(p0, p1, p2);
        triangles.emplace_back(t_p);

        std::shared_ptr<Edge> p1p2 = std::make_shared<Edge>(p1, p2);
        std::shared_ptr<Edge> p2p0 = std::make_shared<Edge>(p2, p0);

        p0p1->left = t_p;
        p1p2->left = t_p;
        p2p0->left = t_p;

        t_p->e_a = p0p1;
        t_p->e_b = p1p2;
        t_p->e_c = p2p0;

        hull.emplace_back(p0p1);
        hull.emplace_back(p1p2);
        hull.emplace_back(p2p0);
    }
    else {
        std::shared_ptr<Triangle> t_p = std::make_shared<Triangle>(p0, p2, p1);
        triangles.emplace_back(t_p);

        std::shared_ptr<Edge> p1p0 = std::make_shared<Edge>(p1, p0);
        std::shared_ptr<Edge> p0p2 = std::make_shared<Edge>(p0, p2);
        std::shared_ptr<Edge> p2p1 = std::make_shared<Edge>(p2, p1);

        p1p0->left = t_p;
        p0p2->left = t_p;
        p2p1->left = t_p;

        t_p->e_a = p0p2;
        t_p->e_b = p2p1;
        t_p->e_c = p1p0;

        hull.emplace_back(p1p0);
        hull.emplace_back(p0p2);
        hull.emplace_back(p2p1);
    }
}

void DelaunayTriangulation(std::vector<std::shared_ptr<float2>> const& points, std::list<std::shared_ptr<Edge>>& hull, std::vector<std::shared_ptr<Triangle>>& triangles) {
    for (auto i = 3; i < points.size(); i++) {
        std::list< std::shared_ptr<Edge>>::iterator it = hull.begin();
        while (it != hull.end()) {
            if (!IsPointOnLeftOfEdge(**it, points[i])) {
                std::shared_ptr<float2> e = points[i];

                std::shared_ptr<Triangle> t = std::make_shared<Triangle>((*it)->start, e, (*it)->end);
                triangles.emplace_back(t);

                t->t_c = (*it)->left;
                if (*((*it)->left->a) == *((*it)->start)) {
                    (*it)->left->t_a = t;
                }
                else if (*((*it)->left->b) == *((*it)->start)) {
                    (*it)->left->t_b = t;
                }
                else if (*((*it)->left->c) == *((*it)->start)) {
                    (*it)->left->t_c = t;
                }

                std::shared_ptr<Edge> edge_low = std::make_shared<Edge>((*it)->start, e);
                edge_low->left = t;
                std::shared_ptr<Edge> edge_high = std::make_shared<Edge>(e, (*it)->end);
                edge_high->left = t;

                Edge rev_edge_low = Edge(e, (*it)->start);
                Edge rev_edge_high = Edge((*it)->end, e);

                t->e_a = edge_low;
                t->e_b = edge_high;

                auto next = std::next(it);
                if (next == hull.end()) {
                    next = hull.begin();
                }
                if (**next == rev_edge_high) {
                    hull.erase(next);
                }
                else {
                    hull.insert(std::next(it), edge_high);
                }
                std::shared_ptr<Edge> prev;
                if (it == hull.begin())
                    prev = hull.back();
                else
                    prev = *std::prev(it);
                bool flag = false;
                if (*prev == rev_edge_low) {
                    if (t->a == (*it)->start) {
                        t->t_a = prev->left;
                        t->e_a = nullptr;
                    }
                    else if (t->b == (*it)->start) {
                        t->t_b = prev->left;
                        t->e_b = nullptr;
                    }
                    else if (t->c == (*it)->start) {
                        t->t_c = prev->left;
                        t->e_c = nullptr;
                    }

                    if (prev->left->a == e) {
                        prev->left->t_a = t;
                        prev->left->e_a = nullptr;
                    }
                    else if (prev->left->b == e) {
                        prev->left->t_b = t;
                        prev->left->e_b = nullptr;
                    }
                    else if (prev->left->c == e) {
                        prev->left->t_c = t;
                        prev->left->e_c = nullptr;
                    }
                    flag = true;
                    if (it == hull.begin())
                        hull.pop_back();
                    else
                        hull.erase(std::prev(it));
                }
                else {
                    hull.insert(it, edge_low);
                }

                if ((*it)->left->e_a && *((*it)->left->e_a) == **it) {
                    (*it)->left->e_a = nullptr;
                } else if ((*it)->left->e_b && *((*it)->left->e_b) == **it) {
                    (*it)->left->e_b = nullptr;
                } else if ((*it)->left->e_c && *((*it)->left->e_c) == **it) {
                    (*it)->left->e_c = nullptr;
                }

                hull.erase(it++);

                Legalize_recursive(t, t->t_c);
                if (flag && it != hull.end()) Legalize_recursive(t, (*it)->left);
            }
            if (it != hull.end()) ++it;
        }
    }
}

void RunDelaunayTriangulation(std::vector<float2> const& points, std::vector<std::shared_ptr<Triangle>> triangles) {
    if (points.size() < 3) {
        return;
    }

    std::vector<std::shared_ptr<float2>> _points;
    for (const auto& s : points) {
        _points.push_back(std::make_shared<float2>(s));
    }
    std::list< std::shared_ptr<Edge>> hull;

    // Initialize, creates first triangle and init hull
    InitializeDelaunay(_points, hull, triangles);
    DelaunayTriangulation(_points, hull, triangles);

    // // Free memory
    // for (auto& ptr : triangles) {
    //     if (ptr->t_a)
    //         ptr->t_a.reset();
    //     if (ptr->t_b)
    //         ptr->t_b.reset();
    //     if (ptr->t_c)
    //         ptr->t_c.reset();
    //     ptr.reset();
    // }
    for (auto& ptr : hull) {
        if (ptr->left)
            ptr->left.reset();
        ptr.reset();
    }
    for (auto& ptr : _points) {
        ptr.reset();
    }
}

}