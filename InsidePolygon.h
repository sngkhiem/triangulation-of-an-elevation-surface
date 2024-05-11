#ifndef INSIDEPOLYGON_H
#define INSIDEPOLYGON_H

#include "qhull_tools.h"
#include<OnSegment.h>
#include<Point.h>

using namespace std;

template<class P>
bool inPolygon(vector<P> &p, P a, bool strict) {
    int cnt = 0, n = p.size();
    for (int i = 0; i < n; i++) {
        P q = p[(i + 1) % n];
        if (onSegment(p[i], q, a)) return !strict;
        cnt ^= ((a.y<p[i].y) - (a.y<q.y)) * a.cross(p[i], q) > 0;
    }
    return cnt;
}

#endif
