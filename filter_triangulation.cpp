#include "filter_triangulation.h"
#include<Point.h>
#include<OnSegment.h>
#include<SegmentDistance.h>
#include<InsidePolygon.h>
#include "qhull_tools.h"
#include <Eigen/Dense>
#include <vcg/complex/algorithms/convex_hull.h>

using namespace std;
using namespace vcg;

const coordT EPS = 1e-9;

struct MyPoint {
	coordT x, y;
};

struct MyEdge
{
	coordT x1, y1;
	coordT x2, y2;
	// adj vertex (which form with this edge to be a triangle)
	coordT adjX, adjY;

	 // Define equality operator for MyEdge struct
	bool operator==(const MyEdge& other) const {
		return (x1 == other.x1 && y1 == other.y1 && x2 == other.x2 && y2 == other.y2) ||
			   (x1 == other.x2 && y1 == other.y2 && x2 == other.x1 && y2 == other.y1);
	}
};

namespace std {
template<>
struct hash<MyEdge>
{
	size_t operator()(const MyEdge& edge) const
	{
		size_t hash_val = 0;
		hash_val ^= hash<coordT> {}(edge.x1) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
		hash_val ^= hash<coordT> {}(edge.y1) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
		hash_val ^= hash<coordT> {}(edge.x2) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
		hash_val ^= hash<coordT> {}(edge.y2) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
		return hash_val;
	}
};
}

QhullPlugin::QhullPlugin()
{
	typeList = {
		FP_TRIANGULATION
	};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QhullPlugin::~QhullPlugin()
{
}

QString QhullPlugin::pluginName() const
{
	return "FilterTriangulation";
}

QString QhullPlugin::filterName(ActionIDType f) const
{
	switch (f) {
	case FP_TRIANGULATION: return QString("Triangulation: Elevation surface structured by a sparse 3D grid");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::pythonFilterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_TRIANGULATION: return QString("Method to triangulate an elevation surface structured in a sparse 3D grid by using a fast search algorithm based on the 2D Delaunay triangulation");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_TRIANGULATION:
		return QString("Method to triangulate an elevation surface structured in a sparse 3D grid by using a fast search algorithm based on the 2D Delaunay triangulation");
	default: assert(0);
	}
	return QString("Error: Unknown Filter");
}

QhullPlugin::FilterClass QhullPlugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_TRIANGULATION: return FilterClass(FilterPlugin::Remeshing);
	default: assert(0);
	}
	return FilterClass(0);
}

RichParameterList QhullPlugin::initParameterList(const QAction* action, const MeshModel& m)
{
	RichParameterList parlst;
	switch (ID(action)) {
	case FP_TRIANGULATION:
	default: break; // do not add any parameter for the other filters
	}
	return parlst;
}

/*
Flat the 2d array of coordinate to 1d array by, one point has three element respectively to its axis
Thus, we need to allocate numpoints * dimension * sizeof(coordT) to save the coordinate of all point
*/
coordT* readpointsFromMesh(int* numpoints, int* dimension, MeshModel& m)
{
	coordT *points, *coords;

	coords = points = (coordT*) malloc((*numpoints) * (*dimension) * sizeof(coordT));

	int cnt = 0;
	CMeshO::VertexIterator vi;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
		if (!(*vi).IsD()) {
			for (int ii = 0; ii < *dimension; ++ii)
				*(coords++) = (*vi).P()[ii];
			++cnt;
		}
	assert(cnt == m.cm.vn);

	return (points);
}

bool isInCircumcircle(coordT pX1, coordT pY1, coordT pX2, coordT pY2, coordT pX3, coordT pY3, coordT pX, coordT pY)
{
	coordT centerX, centerY;
	coordT a1, b1, d1, a2, b2, d2;
	coordT radius, checkDistance;

	// Make the below calculation easier
	a1 = -2.0 * pX1 + 2.0 * pX2;
	b1 = -2.0 * pY1 + 2.0 * pY2;
	d1 = -pX1 * pX1 - pY1 * pY1 + pX2 * pX2 + pY2 * pY2;
	a2 = -2.0 * pX1 + 2.0 * pX3;
	b2 = -2.0 * pY1 + 2.0 * pY3;
	d2 = -pX1 * pX1 - pY1 * pY1 + pX3 * pX3 + pY3 * pY3;

	// Using determinant to calculate the coordinate x, and y of center point
	centerX = (d1 * b2 - d2 * b1) / (a1 * b2 - a2 * b1);
	centerY = (d2 * a1 - d1 * a2) / (a1 * b2 - a2 * b1);
	radius  = sqrt((centerX - pX1) * (centerX - pX1) + (centerY - pY1) * (centerY - pY1));

	// Distance between check point and center point
	checkDistance = sqrt((centerX - pX) * (centerX - pX) + (centerY - pY) * (centerY - pY));
	return checkDistance < radius;
}

// Check if the given point is in the right side of the given edge
/*
	 - Assume the edge is AB and the point is C
	 - Calculate cross product of vector AB * AC
	 - Using determinant we have AB * AC = (x_b-x_a) * (y_c - y_a) - (x_c - x_a) * (y_b - y_a)
	 - If the cross product is > 0, then this point is in the right side of the given edge
	 - If the cross product is < 0, then this point is in the left side of the given edge
	 - If the cross product = 0, then this point and this edge is collinear
*/
int isRight(MyEdge edge, coordT pX, coordT pY) 
{
	if ((edge.x2 - edge.x1) * (pY - edge.y1) - (pX - edge.x1) * (edge.y2 - edge.y1) > 0) // In the right
		return 1; 
	else if ((edge.x2 - edge.x1) * (pY - edge.y1) - (pX - edge.x1) * (edge.y2 - edge.y1) < 0) // In the left
		return 0;
	else if ((edge.x2 - edge.x1) * (pY - edge.y1) - (pX - edge.x1) * (edge.y2 - edge.y1) == 0) // collinear
		return 2;
}

bool compareEqual(coordT a, coordT b)
{
	return abs(a-b) < EPS;
}


// Check if a given point is lied in a given rectangle
bool isIn(MyPoint point, MyPoint corner1, MyPoint corner2, MyPoint corner3, MyPoint corner4)
{
	vector<P> v = { P{corner1.x, corner1.y}, P{corner2.x, corner2.y}, P{corner3.x, corner3.y}, P {corner4.x, corner4.y} };
	return inPolygon(v, P{point.x, point.y}, false);
}

std::map<std::string, QVariant> QhullPlugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* /* cb*/)
{
	qhT  qh_qh = {};
	qhT* qh    = &qh_qh;

	switch (ID(filter)) {
	case FP_TRIANGULATION: {
		MeshModel  &m  = *md.mm();
		MeshModel &nm  = *md.addNewMesh("", "Triangulation 2d");
		MeshModel &fm  = *md.addNewMesh("", "Final Triangulation");
		int dim = 3;
		int numpoints = m.cm.vn;
		int row = sqrt(numpoints);
		int col = sqrt(numpoints);
		coordT *points;
		Eigen::MatrixXi sparseMatrix(row, col);
		Eigen::MatrixXd zAxis(row, col);
		sparseMatrix = Eigen::MatrixXi::Zero(row, col);
		zAxis = Eigen::MatrixXd::Zero(row, col);
		points = readpointsFromMesh(&numpoints, &dim, m);
		
		//Convert 3d surface to 2d grid by eliminating the z coordinate
		for (int i = 0; i < numpoints; i++) {
			Point3d curVertex = {points[i * dim], points[i * dim + 1], 0};
			coordT curX = points[i * dim];
			coordT curY = points[i * dim + 1];
			// Use sparse matrix as a lookup table
			sparseMatrix((int)curX, (int)curY) = 1;
			zAxis((int)curX, (int)curY) = points[i * dim + 2];
			tri::Allocator<CMeshO>::AddVertex(nm.cm, curVertex);
		}
		/* Create seed triangle */
		// (0, 0) (0, 1) and (1, 0)
		Point3d p0 = {0, 0, 0};
		Point3d p1 = {0, 1, 0};
		Point3d p2 = {1, 0, 0};
		// Edges of this seed triangle
		MyEdge seedE3 = {1, 0, 0, 1, 0, 0};
		tri::Allocator<CMeshO>::AddFace(nm.cm, p0, p1, p2);

		/* Triangulating a surface */
		queue<MyEdge> EdgePool;
		unordered_set<MyEdge>   visited;
		//vector<MyEdge> visited;
		visited.insert(seedE3);
		// Add edge of seed triangle into edges pool
		EdgePool.push(seedE3);
		while (!EdgePool.empty()) {
			// Take the front edge
			MyEdge curEdge = EdgePool.front();
			EdgePool.pop();

			// Find all valid neighbor points of this edge
			bool   curRight             = 0;
			coordT directionVectorX     = curEdge.x2 - curEdge.x1;
			coordT directionVectorY     = curEdge.y2 - curEdge.y1;
			coordT perpendicularVectorX = -directionVectorY;
			coordT perpendicularVectorY = directionVectorX;
			coordT unitVectorX =
				perpendicularVectorX / sqrt(
										   perpendicularVectorX * perpendicularVectorX +
										   perpendicularVectorY * perpendicularVectorY);
			coordT unitVectorY =
				perpendicularVectorY / sqrt(
										   perpendicularVectorX * perpendicularVectorX +
										   perpendicularVectorY * perpendicularVectorY);
			coordT extendVectorX, extendVectorY;
			coordT k = sqrt(
				(curEdge.x1 - curEdge.x2) * (curEdge.x1 - curEdge.x2) +
				(curEdge.y1 - curEdge.y2) * (curEdge.y1 - curEdge.y2));
			if (isRight(curEdge, curEdge.adjX, curEdge.adjY)) {
				extendVectorX = -unitVectorX;
				extendVectorY = -unitVectorY;
			}
			else {
				curRight      = 1;
				extendVectorX = unitVectorX;
				extendVectorY = unitVectorY;
			}
			// Find the other two corners of the rectangle Re(e_i)
			coordT corner1X = curEdge.x1 + extendVectorX * k;
			coordT corner1Y = curEdge.y1 + extendVectorY * k;
			coordT corner2X = curEdge.x2 + extendVectorX * k;
			coordT corner2Y = curEdge.y2 + extendVectorY * k;

			// Calculate above and below point of this rectangle
			coordT minX = min({curEdge.x1, curEdge.x2, corner1X, corner2X});
			coordT maxX = max({curEdge.x1, curEdge.x2, corner1X, corner2X});
			coordT minY = min({curEdge.y1, curEdge.y2, corner1Y, corner2Y});
			coordT maxY = max({curEdge.y1, curEdge.y2, corner1Y, corner2Y});

			vector<MyPoint> neighborPoints;
			MyPoint         corner1 = {curEdge.x1, curEdge.y1};
			MyPoint         corner2 = {corner1X, corner1Y};
			MyPoint         corner3 = {corner2X, corner2Y};
			MyPoint         corner4 = {curEdge.x2, curEdge.y2};
			coordT          i       = minX;
			coordT          j       = minY;
			while (true) {
				MyPoint curPoint = {i, j};
				bool    isEdge   = false;
				if (i == corner1.x && j == corner1.y)
					isEdge = true;
				if (i == corner4.x && j == corner4.y)
					isEdge = true;
				if (!isEdge && i >= 0.0 && j >= 0.0 && i <= sqrt(numpoints) - 1 && j <= sqrt(numpoints) - 1 &&
					isIn(curPoint, corner1, corner2, corner3, corner4) &&
					sparseMatrix((int) i, (int) j) == 1) {
					neighborPoints.push_back({i, j});
				}
				if (compareEqual(j, maxY)) {
					i += 1.0, j = minY;
					if (i > maxX)
						break;
				}
				else
					j += 1.0;
			}

			// Find the best point to form new triangle
			coordT bestCompact = -1e18;
			coordT bestPX = -1.0, bestPY = -1.0;
			for (int k = 0; k < neighborPoints.size(); k++) {
				coordT curX = neighborPoints[k].x;
				coordT curY = neighborPoints[k].y;
				// Check if this point achieve delaunay criteration
				bool isFalse = false;
				for (int ii = 0; ii < neighborPoints.size(); ii++) {
					if (isInCircumcircle(curEdge.x1, curEdge.y1, curEdge.x2, curEdge.y2, curX, curY, neighborPoints[ii].x, neighborPoints[ii].y)) {
						isFalse = true;
						break;
					}
				}

				if (!isFalse) {
					coordT l0 = sqrt(
						(curEdge.x1 - curEdge.x2) * (curEdge.x1 - curEdge.x2) +
						(curEdge.y1 - curEdge.y2) * (curEdge.y1 - curEdge.y2));
					coordT l1 = sqrt(
						(curEdge.x1 - curX) * (curEdge.x1 - curX) +
						(curEdge.y1 - curY) * (curEdge.y1 - curY));
					coordT l2 = sqrt(
						(curEdge.x2 - curX) * (curEdge.x2 - curX) +
						(curEdge.y2 - curY) * (curEdge.y2 - curY));
					coordT p          = (l0 + l1 + l2) / 2.0;
					coordT A          = sqrt(p * (p - l0) * (p - l1) * (p - l2));
					coordT curCompact = (4.0 * sqrt(3.0) * A) / (l0 * l0 + l1 * l1 + l2 * l2);
					if (curCompact > bestCompact) {
						bestCompact = curCompact;
						bestPX      = curX;
						bestPY      = curY;
					}
				}
			}

			// Add new triangle
			if (bestPX != -1.0 && bestPY != -1.0) {
				Point3d p0 = {curEdge.x1, curEdge.y1, 0.0};
				Point3d p1 = {curEdge.x2, curEdge.y2, 0.0};
				Point3d p2 = {bestPX, bestPY, 0.0};
				tri::Allocator<CMeshO>::AddFace(nm.cm, p0, p1, p2);
				// Add new edge
				MyEdge newE0 = {curEdge.x1, curEdge.y1, bestPX, bestPY, curEdge.x2, curEdge.y2};
				MyEdge newE1 = {curEdge.x2, curEdge.y2, bestPX, bestPY, curEdge.x1, curEdge.y1};

				if (visited.find(newE0) == visited.end()) {
					EdgePool.push(newE0);
					visited.insert(newE0);
				}

				if (visited.find(newE1) == visited.end()) {
					EdgePool.push(newE1);
					visited.insert(newE1);
				}
			}
		}

		// Delete duplicate vertices and faces
		tri::Clean<CMeshO>::RemoveDuplicateVertex(nm.cm);
		tri::Clean<CMeshO>::RemoveDuplicateFace(nm.cm);

		// Restore z
		CMeshO::FaceIterator fi;
		for (fi = nm.cm.face.begin(); fi != nm.cm.face.end(); fi++) {
			if (!(*fi).IsD()) {
				coordT cur1X, cur1Y, cur2X, cur2Y, cur3X, cur3Y, cur1Z, cur2Z, cur3Z;
				cur1X = (*fi).P(0)[0];
				cur1Y = (*fi).P(0)[1];
				cur1Z = zAxis(cur1X, cur1Y);
				cur2X = (*fi).P(1)[0];
				cur2Y = (*fi).P(1)[1];
				cur2Z = zAxis(cur2X, cur2Y);
				cur3X = (*fi).P(2)[0];
				cur3Y = (*fi).P(2)[1];
				cur3Z = zAxis(cur3X, cur3Y);

				Point3d p0 = {cur1X, cur1Y, cur1Z};
				Point3d p1 = {cur2X, cur2Y, cur2Z};
				Point3d p2 = {cur3X, cur3Y, cur3Z};

				tri::Allocator<CMeshO>::AddFace(fm.cm, p0, p1, p2);
			}
		}

		// Delete duplicate vertices and faces
		tri::Clean<CMeshO>::RemoveDuplicateVertex(fm.cm);
		tri::Clean<CMeshO>::RemoveDuplicateFace(fm.cm);
	} break;
	default: wrongActionCalled(filter);
	}
	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
