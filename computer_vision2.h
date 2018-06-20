#ifndef INC_COMPUTER_VISION2
#define INC_COMPUTER_VISION2

#include <math.h>
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_1d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_2d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_3d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\my_stl.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\graphics_debug.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\graphics_2d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\color.h"
#include <vector>

// interface for reading sequential samples from source
template <class sample_type = rgb>
class sample_reader {
	
	sample_type* source;
	quantum_2d source_dims;
	quantum_2d bot_lef, top_rig;
	unsigned int row_start;
	quantum_2d pos;

	inline unsigned int calc_source_index(int x, int y)
	{
		return x + y * source_dims.x;
	}

public:

	sample_reader (
		void* p_source, quantum_2d p_source_dims,
		quantum_2d p_bot_lef = VECT_2D_INT_ORIG, quantum_2d p_top_rig = p_source_dims - quantum_2d(1, 1))
	{
		source = (sample_type*)p_source;
		source_dims = p_source_dims;
		bot_lef = p_bot_lef;
		top_rig = p_top_rig;
	}

	~sample_reader ()
	{
	}

	quantum_2d get_source_pos()
	{
		return pos;
	}

	bool x_cont()
	{
		return pos.x <= top_rig.x;
	}

	bool y_cont()
	{
		return pos.y <= top_rig.y;
	}

	inline sample_type get_source_sample()
	{
		return source[row_start + pos.x].reverse();
	}

	void init_x()
	{		
		pos.x = bot_lef.x;
	}

	void init_y()
	{
		pos.y = bot_lef.y;
		row_start = calc_source_index(0, bot_lef.y);// + source_dims.x;
	}

	void advance_x()
	{
		pos.x++;
	}

	void advance_y()
	{
		pos.y++;
		row_start += source_dims.x;
	}

	point_2d div_x()
	{
		return get_source_pos().cv_float() + point_2d(0.0f, -0.5f);
	}

	point_2d div_y()
	{
		return get_source_pos().cv_float() + point_2d(0.5f, -1.0f);
	}

	sample_type map(int x, int y)
	{
		return source[x + y * source_dims.x];
	}
};


template <class element = point_2d, class area_type = double>
class corner {
public:
	area_accum<element, area_type> acc;

	corner () 
	{
	}

	corner (area_accum<element, area_type> p_acc)
	{
		acc = p_acc;
	}
};


template <class element = point_2d, class area_type = double>
class corner_detector {
	
	int process(float p_threshhold)
	{
		int res = test(p_threshhold);

		if (res != 0) {
			ref = curr;
		}

		return res;
	}

public:
	// current and reference accumulators
	area_accum<element, area_type> prev, curr, ref;

	void init(element p_first)
	{
		prev = curr = ref = area_accum<element, area_type>(p_first);
	}

	corner_detector()
	{
	}

	corner_detector(element p_first)
	{
		init(p_first);
	}

	int update(element p_next, float p_threshhold)
	{
		prev = curr;
		curr += p_next;
		return process(p_threshhold);
	}

	int update_no_reset(element p_next, float p_threshhold)
	{
		prev = curr;
		curr += p_next;

	}

	int update_rev(element p_next)
	{
		prev = curr;
		curr -= p_next;
		return process();
	}

	area_accum<element, area_type> get_vert()
	{
		return prev;
		//return curr;
	}

	area_type last_area()
	{
		return area_accum<element, area_type>::diff_area(ref.area, curr.area);
	}

	double ave_height()
	{
		return curr.ave_height(ref);
	}

	double ave_height(corner_detector<element> p_next)
	{
		double area = last_area() - p_next.last_area() - area_accum<element, area_type>::new_area(ref.point, p_next.ref.point);
		return area / ref.point.dist(p_next.ref.point);
	}

	int test(double h_bar, double p_threshhold)
	{
		return side_of_threshhold(h_bar, p_threshhold);
	}

	int test(double p_threshhold)
	{
		return test(ave_height(), p_threshhold);
	}

	int test(corner_detector<element> p_next, double p_threshhold = 0.5)
	{
		return test(ave_height(p_next), p_threshhold);
	}
};


template <class element = point_2d, class area_type = double>
class perim_tracer {

	point_2dd sum;
	int ct;
	corner_detector<element, area_type> det;
	std::list<element> samp;

	void init_first(element p_start)
	{
		det = corner_detector<element>(p_start);
		first = p_start;
		area_offset = 0.0;
	}

	void trace_fwd(point_2d p)
	{
		samp.push_back(p);

		if (det.update(p, 0.5) != 0) {
			log_corner_fwd();
		}

		sum += p;
		ct++;
	}

	void trace_rev(point_2d p)
	{
		samp.push_front(p);

		if (det.update(p, 0.5) != 0) {
			log_corner_rev();
		}

		sum += p;
		ct++;
	}

public:

	perim_tracer* next;
	perim_tracer* other_end;
	std::list<corner<>> corners;
	element first;
	bool forward;
	double area_offset;

	perim_tracer(element p_start)
	{
		init_first(p_start);
	}

	void init_left(perim_tracer* p_right)
	{
		forward = false; // reverse
		next = other_end = p_right;
		area_offset = 0.0;
	}

	void init_right(perim_tracer* p_left)
	{
		forward = true; // forward
		next = NULL;
		other_end = p_left;
		area_offset = 0.0;
	}

	void log_corner_fwd(area_accum<element, area_type> acc)
	{
		corners.push_back(corner<>(acc));
	}

	void log_corner_fwd()
	{
		log_corner_fwd(det.get_vert());
	}

	void log_corner_rev(area_accum<element, area_type> acc)
	{
		corners.push_front(corner<>(acc));
	}

	void log_corner_rev()
	{
		log_corner_rev(det.get_vert());
	}

	element last_element()
	{
		return det.curr.point;
	}

	element last_ref()
	{
		return det.ref.point;
	}

	area_type get_area()
	{
		return det.curr.area;
	}

	area_type get_area(perim_tracer* p_right)
	{
		return get_area() + p_right->get_area();
	}

	area_accum<element, area_type> get_orig()
	{
		return area_accum<element, area_type>(first);
	}

	void trace_left_fwd(int p_x, int p_y)
	{
		trace_fwd(point_2d(p_x, p_y).left());
	}

	void trace_left_rev(int p_x, int p_y)
	{
		trace_rev(point_2d(p_x, p_y).left());
	}

	void trace_right_fwd(int p_x, int p_y)
	{
		trace_fwd(point_2d(p_x, p_y).right());
	}

	void trace_right_rev(int p_x, int p_y)
	{
		trace_rev(point_2d(p_x, p_y).right());
	}

	void trace_ref_rev(int p_x, int p_y)
	{
		trace_rev(point_2d(p_x, p_y).top());
	}

	void trace_wrk_edge_fwd_pos(int p_lef_x, int p_rig_x, int p_y)
	{
		for (int x = p_lef_x; x <= p_rig_x; x++) {
			trace_fwd(point_2d(x, p_y).bottom());
		}
	}

	void trace_wrk_edge_fwd_neg(int p_lef_x, int p_rig_x, int p_y)
	{
		for (int x = p_rig_x; x >= p_lef_x; x--) {
			trace_fwd(point_2d(x, p_y).bottom());
		}
	}

	void trace_wrk_edge_rev_pos(int p_lef_x, int p_rig_x, int p_y)
	{
		for (int x = p_lef_x; x < p_rig_x; x++) {
			trace_rev(point_2d(x, p_y).bottom());
		}
	}

	void trace_wrk_edge_rev_neg(int p_lef_x, int p_rig_x, int p_y)
	{
		for (int x = p_rig_x; x >= p_lef_x; x--) {
			trace_rev(point_2d(x, p_y).bottom());
			//show_point(point_2d(x, p_y).center(), RGB_AQUA, 3, false);
		}
	}

	void trace_ref_edge_fwd_pos(int p_lef_x, int p_rig_x, int p_y)
	{
		for (int x = p_lef_x; x <= p_rig_x; x++) {
			trace_fwd(point_2d(x, p_y).top());
		}
	}

	void trace_ref_edge_fwd_neg(int p_lef_x, int p_rig_x, int p_y)
	{
		for (int x = p_rig_x; x >= p_lef_x; x--) {
			trace_fwd(point_2d(x, p_y).top());
		}
	}

	void trace_ref_edge_rev_pos(int p_lef_x, int p_rig_x, int p_y)
	{
		for (int x = p_lef_x; x <= p_rig_x; x++) {
			trace_rev(point_2d(x, p_y).top());
		}
	}

	void trace_ref_edge_rev_neg(int p_lef_x, int p_rig_x, int p_y)
	{
		for (int x = p_rig_x; x >= p_rig_x; x--) {
			trace_rev(point_2d(x, p_y).top());
		}
	}

	// appends p_next to end of this
	bool append(perim_tracer* p_next)
	{
		if (other_end == p_next) {
			return true;
		} else {
			perim_tracer* tmp = other_end;
			other_end->other_end = p_next->other_end;
			p_next->other_end->other_end = tmp;
			next = p_next;
			return false;
		}
	}

	int join_term(perim_tracer* p_next, double p_threshhold = 0.5)
	{
		double h_bar = area_accum<>::ave_height_under(det.ref, det.curr, p_next->det.curr, p_next->det.ref);
		return side_of_threshhold(h_bar, p_threshhold);
	}

	// calculates and stores the difference of area between
	// conjoined tracers.  stores offset on the previous (not next)
	void get_area_offset(perim_tracer* p_next)
	{
		area_offset = get_area() - p_next->get_area();
	}

	int join_spawn(perim_tracer* p_next, double p_threshhold = 0.5)
	{
		double h_bar = area_accum<>::ave_height_under(corners.back().acc, get_orig(), get_orig(), p_next->corners.front().acc);
		return side_of_threshhold(h_bar, p_threshhold);
	}

	// assumes that this is tracer of beginning of circuit
	void report(rgb c = RGB_GREEN)
	{
		std::vector<corner<>> corn;

		// show samples for whole circuit in order
		for (perim_tracer* t = this; t != NULL; t = t->next) {
			for (std::list<element>::iterator i = t->samp.begin(); i != t->samp.end(); ++i) {
				//show_point(*i, c, 3, true);
			}
		}

		area_accum<> prev = POINT_2DD_UNDEF_VALUE;
		line_2d_dydx prev_l = LINE_2D_DYDX_UNDEF;
		line_2d_dydx of;
		line_2d_dydx seg_l;

		// show corners for whole circuit in order
		for (perim_tracer* t = this; t != NULL; t = t->next) {

			for (std::list<corner<>>::iterator i = t->corners.begin(); i != t->corners.end(); ++i) {
				if (!prev.point.is_undef()) {
					show_seg(prev.point, i->acc.point, c, false);
					//segment_2d<> seg = segment_2d<>(prev.point, i->acc.point).offset(prev.ave_height(i->acc) * 1.0);
					//of = seg.line_dydx();
					//seg_l = segment_2d<>(prev.point, i->acc.point).offset(prev.ave_height(i->acc.point) * 1.0).line_dydx();
					//point_2d in = seg_l.intersect(prev_l);
					//show_point(in);

					//show_text("h bar", prev.ave_height(i->acc), 22);
					//show_text("offset", t->area_offset, 23);
					//show_seg(seg.a, seg.b, RGB_RED, true);
					//corn.push_back(*i);
				}
				prev = i->acc;
				//prev_l = of;
			}
		}

		if (!prev.point.is_undef() && corners.size() > 0) {
			show_seg(prev.point, corners.front().acc.point, c, false);
			//segment_2d<> seg = segment_2d<>(prev.point, corners.front().acc.point).offset(prev.ave_height(corners.front().acc, area_offset));
			//show_text("h bar", prev.ave_height(corners.front().acc, area_offset) * 2.0, 22);
			//show_text("offset", t->area_offset, 23);
			//show_seg(seg.a, seg.b, RGB_RED, true);
		}

		//for (int i = 0; i < corn.size(); ++i) {

		//}
	}
};


// segments parsed horizontally within tolerance of cannon
template <class inp_sample_type = rgb, class wrk_sample_type = pixel<float>>
class group_seg {
public:
	typedef segment_2d<int> seg_type;

	seg_type seg;
	wrk_sample_type cannon;
	inp_sample_type raw_cannon;
	pixel<int> sum;
	bool connected;
	perim_tracer<>* left;
	perim_tracer<>* right;

	group_seg()
	{
	}

	group_seg(int p_x, inp_sample_type p_sum, wrk_sample_type p_cannon, inp_sample_type p_raw_cannon)
	{
		seg.a = p_x;
		cannon = p_cannon;
		raw_cannon = p_raw_cannon;
		sum = p_sum;
		connected = false;
	}

	bool test(wrk_sample_type p_sample, float p_threshhold)
	{
		return cannon.euc_dist_thresh(p_sample, p_threshhold);
	}
	
	bool conx_seg(group_seg ref)
	{
		return cannon == ref.cannon;
	}

	wrk_sample_type mean_sample()
	{
		return (sum * (1.0 / double(seg.b - seg.a + 1)));
	}

	// spawn case: parse top and sides of working
	void spawn_case(int p_y)
	{
		// create new left and right tracers
		point_2d orig = point_2d(seg.a, p_y);
		left = new perim_tracer<>(orig);
		right = new perim_tracer<>(orig);
		left->init_left(right);
		right->init_right(left);
		
		// trace this section of perim
		left->trace_left_rev(seg.a, p_y);
		right->trace_wrk_edge_fwd_pos(seg.a, seg.b, p_y);
		right->trace_right_fwd(seg.b, p_y);
	}

	void left_cont_case(std::list<group_seg<>>::iterator ref, int p_y)
	{
		left = ref->left;

		if (seg.a < ref->seg.a) {
			left->trace_wrk_edge_rev_neg(seg.a, ref->seg.a - 1, p_y);
		} else  {
			left->trace_wrk_edge_rev_pos(ref->seg.a, seg.a, p_y);
		}

		left->trace_left_rev(seg.a, p_y);
	}

	void right_cont_case_a(std::list<group_seg<>>::iterator ref, int p_y)
	{
		// tracer right
		right = ref->right;
		right->trace_wrk_edge_fwd_pos(ref->seg.b + 1, seg.b, p_y);
		right->trace_right_fwd(seg.b, p_y);
	}

	void right_cont_case_b(std::list<group_seg<>>::iterator ref, int p_y)
	{
		// trace right
		right = ref->right;
		right->trace_wrk_edge_fwd_neg(seg.b + 1, ref->seg.b, p_y);
		right->trace_right_fwd(seg.b, p_y);
	}

	void diverge_case(std::list<group_seg<>>::iterator p_right, int p_y)
	{
		// create tracers for left and right
		point_2d orig = point_2d(p_right->seg.a - 1, p_y);
		p_right->left = new perim_tracer<>(orig);
		right = new perim_tracer<>(orig);
		p_right->left->init_left(right);
		right->init_right(p_right->left);

		// trace this section of new perim
		right->trace_right_fwd(p_right->seg.a - 1, p_y);
		right->trace_wrk_edge_fwd_neg(seg.b + 1, p_right->seg.a - 1, p_y);
		right->trace_right_fwd(seg.b, p_y);
	}

	void converge_case(std::list<group_seg<>>::iterator p_right, int p_y)
	{
		// trace this section of perim
		right->trace_wrk_edge_fwd_pos(seg.b + 1, p_right->seg.a - 1, p_y);
		right->trace_left_fwd(p_right->seg.a, p_y - 1);

		if (right->join_term(p_right->left) != 0) {
			p_right->left->get_area_offset(right);
			right->log_corner_fwd();
		}

		if (right->append(p_right->left)) {
			// inny circuit is complete, parse results

			// find missing corners
			for (perim_tracer<>* t = p_right->left; t != NULL; t = t->next) {
				//show_point(t->corners.back().acc.point, RGB_BLUE);
				if (!t->forward) {
					if (t->join_spawn(t->next)) {
						t->log_corner_fwd(t->get_orig());
						//show_point(t->first, RGB_YELLOW, 5);
					} else {
						//show_point(t->first, RGB_PURPLE, 5);
					}
				}
			}

			//p_right->left->report(RGB_AQUA); 
		} 
	}

	// terminal case: parse bottom of reference and check for complete circuit
	void term_case(int p_y)
	{
		left->trace_ref_edge_rev_pos(seg.a, seg.b, p_y);
		left->trace_right_rev(seg.b, p_y);

		if (right->join_term(left) != 0) {
			left->get_area_offset(right);
			left->log_corner_rev();
			//show_point(left->last_element(), RGB_RED, 5);
		} else {
			//show_point(left->last_element(), RGB_GREEN, 5);
		}

		if (right->append(left)) {
			// outy circuit complete
			//mess("left area", left->get_area());
			//mess("right area", right->get_area());

			// find missing corners
			for (perim_tracer<>* t = left; t != NULL; t = t->next) {
				//show_point(t->corners.back().acc.point, RGB_BLUE);
				if (!t->forward) {
					if (t->join_spawn(t->next)) {
						t->log_corner_fwd(t->get_orig());
						//show_point(t->first, RGB_YELLOW, 5);
					} else {
						//show_point(t->first, RGB_PURPLE, 5);
					}
				}
			}

			if (fabs(left->get_area(right)) > 0.0) {
				left->report(RGB_GREEN);
			}

			//mess("area left", left->get_area());
			//mess("area right", right->get_area());

		} 
	}
};


template <class element = float>
class area_tracer {
public:
	area_accum<> head, tail;

	void init(point_2d p_p)
	{
		head.init(p_p);
		tail.init(p_p);
	}

	void init(element p_ind, element p_dep)
	{
		init((float)p_ind, (float)p_dep);
	}

	area_tracer(element p_ind, element p_dep)
	{
		init(p_ind, p_ptr);
	}

	void operator ++ ()
	{

	}
};


template <class element = float>
class worm {
public:

	area_accum<> head;
	area_accum<> tail;
	int prev_m_case;

	void init(element p_ind, element p_dep)
	{
		point_2d p = point_2d((float(p_ind), float(p_dep)));
		head.init(p);
		tail.init(p);
		prev_m_case = -2;
	}

	worm()
	{
	}

	worm(element p_ind, element p_dep)
	{
		init(p_ind, p_dep);
	}

	double ave_height()
	{
		return tail.ave_height2(head);
	}

	int test(double p_threshhold)
	{
		return side_of_threshhold(ave_height(), p_threshhold);
	}

	void update_head(element p_dep)
	{
		head.inc_ind(p_dep);
	}

	void update_tail(element p_dep)
	{
		tail.inc_ind(p_dep);
	}

	int head_ind()
	{
		return (int)head.point.x;
	}

	int tail_ind()
	{
		return (int)tail.point.x;
	}

	double slope()
	{
		return head.point.slope(tail.point);
	}

	double alt_slope()
	{
		return head.point.alt_slope(tail.point);
	}

	int slope_case(double p_threshhold)
	{
		return side_of_threshhold(slope(), p_threshhold);
	}

	int delta_slope_case(double p_threshhold)
	{
		int m_case = slope_case(p_threshhold);
		int ret = m_case;

		if (prev_m_case == m_case) {
			ret = 0;
		}

		prev_m_case = m_case;
		return ret;
	}
};


// containter of pre-parsed segments 
template <class inp_sample_type = rgb, class wrk_sample_type = pixel<float>, class contain = std::list<group_seg<inp_sample_type, wrk_sample_type>>>
class group_seg_list {

	void init_seg(int p_x)
	{
		wrk_seg = group_seg<inp_sample_type, wrk_sample_type>(seg_start(p_x), 0, sample, raw_sample);
	}

	wrk_sample_type conv_inp_sample(inp_sample_type p_sample)
	{
		return p_sample.color_percept3();
	}

public:

	contain* working;
	contain* reference;
	group_seg<> wrk_seg;
	(typename contain)::iterator ref_iter;
	wrk_sample_type sample;
	inp_sample_type raw_sample;
	int start_x;
	int end_x;
	int y;
	wrk_sample_type* samp_hist;
	inp_sample_type* raw_hist;

	group_seg_list(int p_start_x, int p_end_x, int p_y)
	{
		start_x = p_start_x;
		end_x = p_end_x;
		y = p_y;
		working = new contain;
		reference = new contain;
		samp_hist = new wrk_sample_type[end_x - start_x + 1];
		raw_hist = new inp_sample_type[end_x - start_x + 1];
	}

	~group_seg_list()
	{
		delete working;
		delete reference;
		delete samp_hist;
		delete raw_hist;
	}

	inline int seg_start(int p_x)
	{
		return p_x;
	}

	inline int seg_end(int p_x)
	{
		return p_x - 1;
	}

	void set_sample(inp_sample_type p_sample, int p_x)
	{
		wrk_seg.sum += raw_sample;
		sample = conv_inp_sample(p_sample);
		raw_sample = p_sample;
		samp_hist[p_x] = sample;
		raw_hist[p_x] = raw_sample;
	}
		
	void init_row(int p_y, inp_sample_type p_sample)
	{
		raw_sample = 0;//p_sample;
		sample = 0;//conv_inp_sample(p_sample);
		init_seg(start_x);
		y = p_y;

		swap();

		// clear working list about to be used
		working->clear();
		ref_iter = reference->begin();
	}

	void cut_seg(int p_x)
	{
		// p_x is the first alias of the new segment
		push_working_prot(seg_end(p_x));
		init_seg(p_x);
	}

	void copy_cannon()
	{
		wrk_seg.cannon = ref_iter->cannon;
		wrk_seg.raw_cannon = ref_iter->raw_cannon;
		wrk_seg.connected = true;
	}

	void copy_cannon_con()
	{
		if (!wrk_seg.connected) {
			copy_cannon();
		}
	}

	void end_row(int p_x)
	{
		//wrk_seg.seg.b = p_x - 1;
		//working->push_back(wrk_seg);
		push_working(p_x - 1);
	}

	void swap()
	{
		// swap working and reference lists
		contain* temp = working;
		working = reference;
		reference = temp;
	}

	bool test_working(float p_threshhold)
	{
		return wrk_seg.test(sample, p_threshhold);
	}

	bool update_ref(quantum_2d p_p)
	{
		bool moved = false;

		if (!reference->empty()) {
			while (ref_iter->seg.b < p_p.x) {
				++ref_iter;
				moved = true;
			}
		}

		return moved;
	}

	bool test_reference(quantum_2d p_p, float p_threshhold)
	{
		return !reference->empty() && ref_iter->test(sample, p_threshhold);
	}

	void push_working(int x)
	{
		wrk_seg.seg.b = x;
		working->push_back(wrk_seg);
	}

	void push_working_prot(int x)
	{
		if (wrk_seg.seg.a <= x) {
			push_working(x);
		}
	}

	void resize_working(int x, pixel<int> accum_sub)
	{
		if (working->back().seg.a > x) {
			working->pop_back();
		} else {
			working->back().sum -= accum_sub;
			working->back().seg.b = x;
		}
	}

	void backtrack(int limit, int x, float p_threshhold)
	{
		int first_x = x;
		pixel<int> color_shift_wrk = 0;
		pixel<int> color_shift_his = 0;
		int good = 0;
		int bad = 0;

		while (
			!(working->back().connected && (x - 1) <= working->back().seg.b) 
			&& (x - 1) >= limit
			&& ref_iter->test(samp_hist[x - 1 - start_x], p_threshhold)) {

			if (x - 1 == working->back().seg.a - 1) {
				working->pop_back();
			}

			if (x - 1 == working->back().seg.b) {
				color_shift_his = 0;
				bad = 0;
			} 

			color_shift_wrk += raw_hist[x - 1 - start_x].cv_int();
			color_shift_his += raw_hist[x - 1 - start_x].cv_int();
			x--;
			bad++;
		}

		if (x >= wrk_seg.seg.a) {
			wrk_seg.sum -= color_shift_wrk;
			push_working_prot(x - 1);
		} else {
			good = 0;
			for (int i = working->back().seg.b; i >= x; i--) {
				good++;
			}

			if (abs(good - bad) > 0) {
				show_point(point_2d(x, y).center(), RGB_RED, 8, true);
			}

			resize_working(x - 1, color_shift_his);
		}

		wrk_seg.sum = color_shift_wrk;
		wrk_seg.seg.a = x;
	}

	void parse1(sample_reader<inp_sample_type>& samp, float threshhold)
	{
		if (!test_working(threshhold)) { // test horizontally
			cut_seg(samp.get_source_pos().x);
		}

		update_ref(samp.get_source_pos());
	
		if (test_reference(samp.get_source_pos(), threshhold)) { // test vertically
			copy_cannon_con();
		}
	}

	void parse2(sample_reader<inp_sample_type>& samp, float threshhold)
	{
		if (!test_working(threshhold)) { // horizontally DISconnected
			cut_seg(samp.get_source_pos().x);
		}

		if (!wrk_seg.connected) { // is not connected to any references
			update_ref(samp.get_source_pos());

			if (test_reference(samp.get_source_pos(), threshhold)) { // vertically connected here
				backtrack(start_x, samp.get_source_pos().x, threshhold);
				copy_cannon();
			}
		}
	}

	void parse_row(sample_reader<inp_sample_type>& samp, float threshhold)
	{
		init_row(samp.get_source_pos().y, samp.get_source_sample());

		for (samp.init_x(); samp.x_cont(); samp.advance_x()) {
			set_sample(samp.get_source_sample(), samp.get_source_pos().x);
			parse2(samp, threshhold);
		}

		end_row(samp.get_source_pos().x);
	}

	inline int wrk_y()
	{
		return y;
	}

	inline int ref_y()
	{
		return y - 1;
	}

	void color_segs(sample_reader<inp_sample_type>& samp, quantum_2d bot_lef, quantum_2d top_rig)
	{
		int prev_b = -1;

		for (std::list<group_seg<>>::iterator i = working->begin(); i != working->end(); ++i) {
			color_seg_horz(i->seg.a, i->seg.b, samp.get_source_pos().y, i->mean_sample());
			//color_seg_horz(i->seg.a, i->seg.b, samp.get_source_pos().y, i->raw_cannon);
			//color_seg_horz(i->seg.a, i->seg.b, samp.get_source_pos().y, i->raw_cannon.false_color());
			//color_seg_horz(i->seg.a, i->seg.b, samp.get_source_pos().y, RGB_AQUA);

			if (i->seg.a > i->seg.b) { // if segment has non-positive length
				show_point(point_2d(i->seg.a, samp.get_source_pos().y).center(), RGB_BLUE, 3, false);
				show_point(point_2d(i->seg.b, samp.get_source_pos().y).center(), RGB_BLUE, 3, true);
			}

			if (prev_b == -1) {
				if (i->seg.a != bot_lef.x) { // if first segment doesn't butt up against left of window
					show_point(point_2d(i->seg.b, samp.get_source_pos().y).center(), RGB_RED, 3, true);
				}
			} else {
				if (i->seg.a != prev_b + 1) { // if segment doesn't butt up against previous segment
					show_point(point_2d(i->seg.a, samp.get_source_pos().y).center(), RGB_ORANGE, 3, false);
					show_point(point_2d(prev_b, samp.get_source_pos().y).center(), RGB_ORANGE, 2, true);
				}
			} 

			prev_b = i->seg.b;
		}

		if (prev_b != top_rig.x) {
			show_point(point_2d(prev_b, samp.get_source_pos().y).center(), RGB_YELLOW, 3, true);
		}
	}

	void trace_first(sample_reader<inp_sample_type>& samp, float threshhold)
	{
		samp.init_y();

		parse_row(samp, threshhold);

		// generate spawn case for every first segment
		for (std::list<group_seg<>>::iterator i = working->begin(); i != working->end(); i++) {
			i->spawn_case(wrk_y());
		}

		samp.advance_y();
	}

	void trace_last()
	{
		for (std::list<group_seg<>>::iterator i = working->begin(); i != working->end(); i++) {
			i->term_case(wrk_y());
		}
	}

	void trace_segs()
	{
		std::list<group_seg<>>::iterator wrk = working->begin();
		std::list<group_seg<>>::iterator ref = reference->begin();
		std::list<group_seg<>>::iterator wrk_con = working->end();
		std::list<group_seg<>>::iterator ref_con = reference->end();

		if (reference->empty()) {
			// generate spawn case for every first segment
			while (wrk != working->end()) {
				wrk->spawn_case(wrk_y());
				wrk++;
			} 
		} else {
			while (wrk != working->end() && ref != reference->end()) {
				// if working and reference connect in this iteration
				if (wrk->conx_seg(*ref)) {
					// if working has at least 1 other connection (not including this iteration)
					if (wrk_con != working->end()) {
						wrk_con->converge_case(ref, wrk_y());
					// if this is only working connection and reference has previous connection
					} else if (ref_con != reference->end()) {
						ref_con->diverge_case(wrk, wrk_y());
					// if this is first connection for reference and working
					} else {
						// continue case: parse left of working
						//show_point(point_2d(wrk->seg.a, wrk_y()).center(), RGB_WHITE, 3, false);
						wrk->left_cont_case(ref, wrk_y());
					}

					wrk_con = ref;
					ref_con = wrk;
				}

				std::list<group_seg<>>::iterator wrk_bak = wrk;
				std::list<group_seg<>>::iterator ref_bak = ref;

				// if working seg ends
				if (wrk_bak->seg.b <= ref_bak->seg.b) {
					// if nothing connected to working when it ends
					if (wrk_con == working->end()) {
						wrk_bak->spawn_case(wrk_y());
					} else if (wrk_con->seg.b <= wrk_bak->seg.b) {
						wrk_bak->right_cont_case_a(wrk_con, wrk_y());
					}

					wrk_con = working->end();
					wrk++;
				}

				// if reference seg ends
				if (wrk_bak->seg.b >= ref_bak->seg.b) {
					// if nothing connectd to reference when it ends
					if (ref_con == reference->end()) {
						ref_bak->term_case(ref_y());
					//} else if (wrk_bak->seg.b <= wrk_con->seg.b) {
					} else if (ref_con->seg.b < ref_bak->seg.b) {
						// right continue case b
						ref_con->right_cont_case_b(ref_bak, wrk_y());
					}

					ref_con = reference->end();
					ref++;
				} 
			}
		}
	}
};


void group_filter(
	draw_surface* in,
	quantum_2d bot_lef = VECT_2D_INT_ORIG, quantum_2d top_rig = quantum_2d(-1, -1), 
	float threshhold = 255.0 * 0.01 * 8.0) 
{
	typedef rgb sample_type;

	if (top_rig.x < 0) {
		top_rig = in->get_top_rig();
	}

	sample_reader<sample_type> samp((rgb*)in->get_direct_ptr(), in->get_dims(), bot_lef, top_rig);
	group_seg_list<> segs(bot_lef.x, top_rig.x, 0);

	for (samp.init_y(); samp.y_cont(); samp.advance_y()) {
		segs.parse_row(samp, threshhold);
		segs.color_segs(samp, bot_lef, top_rig);
		segs.trace_segs();
	}

	segs.trace_last();
}


void delta_filter(
	draw_surface* in,
	quantum_2d bot_lef = VECT_2D_INT_ORIG, quantum_2d top_rig = quantum_2d(-1, -1), 
	float threshhold = 255.0 * 0.01 * 0.15) 
{
	typedef rgb sample_type;

	if (top_rig.x < 0) {
		top_rig = in->get_top_rig();
	}

	//show_point(point_2d(50, 50), RGB_RED, 3, true);

	sample_reader<sample_type> samp((rgb*)in->get_direct_ptr(), in->get_dims(), bot_lef, top_rig);
	worm<>* worm_y = new worm<float>[in->get_dims().x];

	for (int x = 0; x < in->get_dims().x; x++) {
		worm_y[x].init((float)bot_lef.y, 0.0f);
	}

	for (int y = 0; y < top_rig.y; ++y) {

		worm<> worm_x((float)bot_lef.x, 0.0f);
		int i = 0;

		//for (samp.init_x(); worm_x.head_ind() < top_rig.x; samp.advance_x()) {
		for (int x = 0; x < top_rig.x * 2; ++x) {
			float t = 0.5f;
			double mt = 15.0;
			int m;
			bool hit = false;

			i = worm_x.head_ind();

			int mid = (worm_x.head_ind() + worm_x.tail_ind()) / 2;

			if (worm_x.test(t) || worm_x.head_ind() >= top_rig.x) {
				worm_x.update_tail((float)samp.map(worm_x.tail_ind(), y).sum());
			} else {
				m = worm_x.slope();

				//if (!approx_zero(m)) {
					//show_text("m", m, 15);
				color_point(point_2d(mid, y), RGB_WHITE.interpolate_sign(m, RGB_GREEN, RGB_RED), 3, false);
				//}

				if (true) {
					//show_point(point_2d(i - 1, y).left(), RGB_ORANGE, 3, false);
					hit = true;
				} else if (m == -1) {
					//show_point(point_2d(i - 1, y).left(), RGB_ORANGE, 3, false);
				}

				while (worm_y[i].test(t) || worm_y[i].head_ind() >= top_rig.y) {
					worm_y[i].update_tail((float)samp.map(i, worm_y[i].tail_ind()).color_percept2().sum());
				};

				m = worm_y[i].slope();

				//if (!approx_zero(m)) {
					//show_text("m", m, 16);
					//show_point(point_2d(i, worm_y[i].head_ind() - 1), RGB_WHITE.interpolate_sign(m, RGB_GREEN, RGB_RED), 3, false);
				//}

				if (true) {
					//show_point(point_2d(i, worm_y[i].head_ind() - 1).bottom(), RGB_ORANGE, 3, false);
					hit = true;
				} else if (m == -1) {
					//show_point(point_2d(i, worm_y[i].head_ind() - 1).bottom(), RGB_ORANGE, 3, false);
				}

				if (hit) {
					//show_point(point_2d(i - 1, worm_y[i].head_ind() - 1).center(), RGB_ORANGE, 3, false);
				}

				worm_y[i].update_head((float)samp.map(i, worm_y[i].head_ind()).color_percept2().sum());
				worm_x.update_head((float)samp.map(i, y).sum());
			}
			
			/*
			show_text("m(y)", worm_y[i].slope(), 17);
			cull_geo();
			show_point(point_2d(i, worm_y[i].tail_ind()).bottom(), RGB_BLUE);
			show_point(point_2d(i, worm_y[i].head_ind()).bottom(), RGB_ORANGE);
			*/

		}
	}

	delete worm_y;

	show_point(0);
}


#endif // INC_COMPUTER_VISION2