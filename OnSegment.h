#ifndef ONSEGMENT_H
#define ONSEGMENT_H

#include "qhull_tools.h"
#include<Point.h>

using namespace std;

template<class P> bool onSegment(P s, P e, P p) {
    return p.cross(s, e) == 0 && (s - p).dot(e - p) <= 0;
}

#endif
