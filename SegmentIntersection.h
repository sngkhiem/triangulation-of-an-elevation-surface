#ifndef SEGMENTINTERSECTION_H
#define SEGMENTINTERSECTION_H

#include <Point.h>
#include <OnSegment.h>

using namespace std;

template<class P> bool segInter(P a, P b, P c, P d)
{
	auto oa = c.cross(d, a), ob = c.cross(d, b), oc = a.cross(b, c), od = a.cross(b, d);
	if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
		return true;
	return false;
}

#endif
