# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


# Only build if we have qhull
if(TARGET external-qhull)

	set(SOURCES filter_triangulation.cpp qhull_tools.cpp)

	set(HEADERS filter_triangulation.h qhull_tools.h Point.h OnSegment.h SegmentDistance.h InsidePolygon.h SegmentIntersection.h)

	add_meshlab_plugin(filter_triangulation ${SOURCES} ${HEADERS})

	target_link_libraries(filter_triangulation PRIVATE external-qhull)

else()
	message(STATUS "Skipping filter_qhull - missing qhull")
endif()
