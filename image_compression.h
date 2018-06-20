#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_1d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_2d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\color.h"

template <class primary = int, class position = int>
class moment {
	primary prim;
	position pos;

	moment(primary p_prim, position p_pos)
	{
		prim = p_prim;
		pos = p_pos;
	}

	point_2d cv_point()
	{
		return point_2d(prim, center(pos));
	}
};

template <class vert_type = point_2d, class area_type = double>
class area_accum {

public:

	area_type area;
	vert_type point;

	inline void init(vert_type first, area_type p_area = (area_type)0)
	{
		area = p_area;
		point = first;
	}

	area_accum ()
	{
	}

	area_accum (vert_type first, area_type p_area = (area_type)0)
	{
		init(first, p_area);
	}

	// defines which direction is considered positive or negative area
	// this will give a net positive area to a polygon parsed in the positive angle direction (counter-clockwise)
	// v1 is before v2 in the positive angle direction
	static area_type new_area(vert_type v1, vert_type v2)
	{
		return trapizoid_area(v2, v1);
	}

	// a1 is before a2 in the positive angle direction
	// this should not have any impact on the sign of the area of the polygon 
	static area_type diff_area(area_type a1, area_type a2)
	{
		return a2 - a1;
	}
		
	// return area under two area accumulators that have joined (triagle-shaped area)
	static area_type area_under(area_accum p1a, area_accum p1b, area_accum p2a, area_accum p2b)
	{
		return (diff_area(p1a.area, p1b.area) + diff_area(p2a.area, p2b.area)) - new_area(p1a.point, p2b.point);
	}

	// returns height under two area accumulators that have joined (triangle-shaped area)
	static double ave_height_under(area_accum p1a, area_accum p1b, area_accum p2a, area_accum p2b)
	{
		return area_under(p1a, p1b, p2a, p2b) / p1a.point.dist(p2b.point);
	}

	/*
	area_type new_area(area_type p_start)
	{
		return new_area(p_start, point);
	}

	area_type diff_area(vert_type p_start)
	{
		return diff_area(p_start, area);
	}
	*/

	area_accum operator + (vert_type v)
	{
		return area_accum(v, area + new_area(point, v));
	}

	area_accum operator += (vert_type v)
	{
		return *this = *this + v;
	}

	area_accum operator - (vert_type v)
	{
		return area_accum(point, area - new_area(point, v));
	}

	area_accum operator -= (vert_type v)
	{
		return *this = *this - v;
	}

	void inc_ind(float p_dep)
	{
		point.x++;
		point.y = p_dep;
		area += double(p_dep);
	}

	// returns accumulated area under connecting line segment from this to next
	area_type area_under (area_accum p_start)
	{
		return diff_area(p_start.area, area) - new_area(p_start.point, point);
	}

	// returns accumulated area under connecting line segment from this to next
	area_type area_under2 (area_accum p_start)
	{
		//show_text(diff_area(p_start.area, area), "diff area", 10);
		//show_text(new_area(p_start.point, point), "new area", 11);

		return diff_area(p_start.area, area) + new_area(p_start.point, point);
	}

	// returns average height of area under connecting line segment from this to next
	double ave_height(area_accum p_start, area_type p_offset = 0.0)
	{
		return (p_offset + area_under(p_start)) / point.dist(p_start.point);
	}

	double ave_height2(area_accum p_start)
	{
		return area_under2(p_start) / point.dist(p_start.point);
	}

	double ave_height(area_accum p_prev, area_accum p_next)
	{
		return area_under(p_prev, p_next) / p_prev.point.dist(p_next.point);
	}

	inline void add_col_x2 (vert_type p)
	{
		area += trapizoid_area_x2(p, point);
		point = p;
	}

	inline void add_col (vert_type p)
	{
		area += trapizoid_area(p, point);
		point = p;
	}

	inline void sub_col_x2(vert_type p)
	{
		area -= trapizoid_area_x2(p, point);
		point = p;
	}

	inline void sub_col(vert_type p)
	{
		area -= trapizoid_area(p, point);
		point = p;
	}

	void iterate(vert_type v)
	{
		add_col(v);
	}

	void iterate_rev(vert_type v)
	{
		sub_col(v);
	}
};

