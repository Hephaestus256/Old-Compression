#ifndef INC_RUN_LEN
#define INC_RUN_LEN

#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\bit_file.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\varint.h"
#include <deque>

// min_span >= 2 * sizeof(run_type) + 1

template <class data_type = int, data_type init_value = -1, class run_type = short int, int min_span = 5, int data_bits = 8>
class run_len2 {

	class run {
	public:
		data_type data;
		int span;

		run ()
		{
		}

		run (data_type p_data, int p_span)
		{
			data = p_data;
			span = p_span;
		}

		void init()
		{
			data = init_value;
			span = 0;
		}
	};

	class run_point {
	public:
		int x, y;
		run r;
		int pos;

		run_point ()
		{
		}

		run_point (const int p_x, const int p_y, const data_type p_data, const int p_run, int p_pos)
		{
			//cout << "def y: " << p_y << endl;
			//getchar();

			x = p_x - p_run;
			y = p_y;
			r = run(p_data, p_run);
			pos = p_pos;
		}
	};

	class run_trace {
	public:
		run curr;
		deque<run_point> spans;
		deque<run_point> skips;
		int pos;

		run_trace ()
		{
			init();
		}

		void init ()
		{
			curr.init();
			pos = 0;
		}

		bool update(const int p_x, const int p_y, const data_type p_data)
		{
			bool ret;

			if (curr.data == p_data) {
				ret = false;
				curr.span++;
			} else {
				ret = is_span();
				if (ret) cut_span(p_x, p_y);
				curr.span = 0;
			}

			curr.data = p_data;

			return ret;
		}

		bool is_span ()
		{
			return curr.span >= min_span - 1;
		}

		void cut_span (const int p_x, int const p_y)
		{
			spans.push_back(run_point(p_x, p_y, curr.data, curr.span + 1, pos));
			//skips.push_back(run_point(
			/*
			if (curr.span > 0) {
				spans.push_back(run_point(p_x, p_y, curr.data, curr.span + 1));
			} else {
				spans.push_back(run_point(p_x, p_y, curr.data, 0));
			}
			*/
		}

		int get_state ()
		{
			if (curr.span == 0) {
				return 0; // 
			} else if (curr.span == min_span - 1) {
				return 1;
			} else {
				return 2;
			}
		}
	};

	run_trace horz;
	//run_trace* vert;
	int x, y, offset;
	int width, height;
	int data_count;
	int span_index;
	deque<data_type> data;
	deque<int> meta;

	void init()
	{
		stream.clear();
		x = 0;
		y = 0;
		offset = 0;
	}

public:
	
	deque<data_type> stream;

	run_len2 (const int p_width, const int p_height)
	{
		width = p_width;
		height = p_height;
		//vert = new run_trace[width];
	}

	~run_len2 ()
	{
		//delete vert;
	}

	// encode functions

	void encode_start ()
	{
		init();
	}

	void encode_next (data_type p_data)
	{
		stream.push_back(p_data);
		horz.update(x, y, p_data);
		//vert[x].update(x, y, p_data);
		x++;

		if (x >= width) {
			horz.cut_span(x, y);
			horz.init();
			x = 0;
			y++;
		}
	}

	void encode_end ()
	{
		horz.cut_span(x, y);
		
		for (int i = 0; i < width; i++) {
			//vert[i].cut_span(i, y - 1);
		}

		for (int i = 0; i < horz.spans.size(); i++) {
			//cout << "span dump: x: " << horz.spans[i].x << " y: " << horz.spans[i].y << " span: " << horz.spans[i].r.span << endl;
			//getchar();
		}
	}

	void process_meta_start ()
	{
		x = 0;
		y = 0;
		span_index = 0;
	}

	int calc_offset(int p_x, int p_y)
	{
		return p_x + width * p_y;
	}

	int calc_x(int p_offset)
	{
		return p_offset % width;
	}

	int calc_y(int p_offset)
	{
		return p_offset / width;
	}

	bool process_meta_next (int& p_run)
	{
		if (y >= height) return false;

		// spit out spans and gaps between spans.  spans are +, gaps are -	
		if (span_index >= horz.spans.size()) {
			p_run = width - x;
			x = width;
		} else {
			run_type s;
			int is_skip = 0;
			bool rep;
			p_run = 0;

			/*
			if (horz.spans[span_index].r.span == 0 && span_index < horz.spans.size()) {
				x += horz.spans[span_index].r.span;
				span_index++;
			}
			*/
				//cout << "y position: " << horz.spans[span_index].y << " span_index " << span_index << endl;

				//if (calc_offset(x, y) >= calc_offset(horz.spans[span_index].x, horz.spans[span_index].y)) {
				/*
				if (horz.spans[span_index].r.span == 0) {
					x = 0;
					span_index++;
				}
				*/

				if (x >= horz.spans[span_index].x) {
					s = horz.spans[span_index].r.span; // run
					p_run += s;
					x += s;
					//offset += s;
					span_index++;
					
					//cout << "x: " << x << " y: " << y << " run: " << p_run << " sx: " << horz.spans[span_index].x << " span index: " << span_index << endl;
					//getchar();
					//x = calc_x(offset);
					//y = calc_y(offset);
				} else {
					s = -(horz.spans[span_index].x - x); // skip
					p_run += s;
					x = horz.spans[span_index].x;

					//y = horz.spans[span_index].y;
					//offset = calc_offset(x, y);
					
					//cout << "x: " << x << " y: " << y << " skip: " << p_run << " sx: " << horz.spans[span_index].x << " span index: " << span_index << endl;
					//getchar();
				}


			//do  {

			//} while (s == 0);
		}

		if (x >= width) {
			x = 0;
			y++;
		}

		return true;
	}

	void process_meta_end ()
	{
	}

	void process_data_start ()
	{
		x = 0;
		y = 0;
		span_index = 0;
	}

	bool process_data_next (data_type& p_data)
	{
		if (y >= height) return false;

		p_data = stream[x + width * y];

		// only represent span data in stream once and skip x to end
		if (span_index < horz.spans.size() && horz.spans[span_index].x == x && horz.spans[span_index].y == y) {
			x += horz.spans[span_index].r.span;
			span_index++;
		} else {
			x++;
		}

		if (x >= width) {
			x = 0;
			y++;
		}

		return true;
	}

	void process_data_end ()
	{
	}

	void decode_meta_start ()
	{
		x = 0;
		meta.clear();
	}

	bool decode_meta_next (int p_run)
	{
		meta.push_back(p_run);

		x += abs(p_run);

		return x < width * height;
	}

	void decode_meta_end ()
	{
	}

	int decode_data_start ()
	{
		// calculate data points
		data_count = 0;
		for (int i = 0; i < meta.size(); ++i) {
			int m = meta[i];

			if (m > 0) {
				data_count++;
			} else if (m < 0) {
				data_count += -m;
			}
		}

		data.clear();

		return data_count;
	}

	bool decode_data_next (data_type p_data)
	{
		data.push_back(p_data);
		return data.size() < data_count;
	}

	void decode_data_end ()
	{
	}

	void process_encode (bit_file& f)
	{
		data_type dat;
		int run;

		cout << "span count: " << horz.spans.size() << endl;
		getchar();

		// process metadata

		process_meta_start();

		while (process_meta_next(run)) {
			//cout << "value: " << run << endl;
			f.write((run_type)run);
		}

		process_meta_end();

		// process data 

		process_data_start();

		while (process_data_next(dat)) {
			f.write(dat, data_bits);
		};

		process_data_end();
	}

	int process_decode ()
	{
		int d = 0;

		stream.clear();

		for (int i = 0; i < meta.size(); ++i) {
			int m = meta[i];

			if (m > 0) { // run
				for (int j = 0; j < m; ++j) {
					stream.push_back(data[d]);
				}
				d++;
			} else if (m < 0) { // literal (skip)
				for (int j = 0; j < -m; ++j) {
					stream.push_back(data[d]);
					d++;
				}
			}
		}

		return (int)stream.size();
	}

	// export to container
	template <class gen>
	void export_cont(deque<gen>& p_cont)
	{
		for (int i = 0; i < stream.size(); ++i) {
			p_cont.push_back(stream[i]);
		}
	}
};

template <int data_bits = 8, int run_bits = data_bits, class data_type = int, class run_cont = int, int limit = -1, int rep_min = 2>
class run_len {

	typedef varint<run_bits, true, run_cont, 0, 1> run_type;

	std::deque<data_type> data_buffer;
	std::deque<run_type> run_buffer;
	data_type prev;
	int run;
	int span;
	bool first;
	bool rep;
	int data_count;
	bool eos; // end of stream

	int update(data_type d)
	{
		if (!first) {
			if (prev == d) {
				run++;
			} else {
				run = 0;
			}
		}

		first = false;
		prev = d;
		return run;
	}

	void encode_zero_run()
	{
		run_type run = 0;
		for (int i = 0; run.read_continue(i); ++i) {
			run_buffer.push_back(run.read_seg(i));
		};
	}

	void encode_rep_run(const int p_run)
	{
		if (p_run) {
			//cout << "rep: " << p_run << endl;
			//getchar();

			run_type run = p_run;
			for (int i = 0; run.read_continue(i); ++i) {
				run_buffer.push_back(run.read_seg(i));
			};
		}
	}

	void encode_nonrep_run(const int p_run)
	{
		if (p_run) {
			//cout << "non-rep: " << p_run << endl;
			//getchar();
		}

		encode_rep_run(-p_run);
	}

public:

	void init()
	{
		run = 0;
		first = true;
		rep = false;
		span = 0;
		data_count = 0;
		eos = false;
	}

	run_len ()
	{
		init();
	}

	/*** Encode methods ***/

	void encode_start()
	{
		init();
	}

	void encode_next(data_type p_data)
	{
		int len = update(p_data);

		if (len >= rep_min) { // repeating
			if (!rep) { // going from non-repeating to repeating
				data_count = (int)data_buffer.size() - rep_min;
				encode_nonrep_run(data_count);
				span = 0;
			}

			span++;
			rep = true;
		} else { // non-repeating
			if (rep) { // going from repeating to non-repeating
				data_count = 1;

				for (int i = 0; i < rep_min - 1; i++) {
					data_buffer.pop_front();
				}

				encode_rep_run(span + rep_min);
				span = 0;
			} else if (limit > 0 && ((int)data_buffer.size() - rep_min) >= limit) {
				data_count = (int)data_buffer.size() - rep_min;
				encode_nonrep_run(data_count);
			}

			data_buffer.push_back(p_data);

			rep = false;
		}
	}

	// returns number of bits encoded
	template <class gen>
	int encode_process(gen& p_data)
	{
		if (run_buffer.size() > 0) {
			p_data = (gen)run_buffer.front();
			run_buffer.pop_front();
			return run_bits;
		} else if (data_buffer.size() > 0 && data_count > 0) {
			p_data = (gen)data_buffer.front();
			data_buffer.pop_front();
			data_count--;
			return data_bits;
		} else if (eos) {
			p_data = run_type(0).container;
			eos = false;
			return run_bits;
		} else {
			return 0;
		}		
	}

	void encode_end()
	{
		if (rep) {
			data_count = 1;
			encode_rep_run(span + rep_min);
		} else {
			data_count = (int)data_buffer.size();
			encode_nonrep_run(data_count);
		}

		eos = true;
	}

	void encode_end(bit_file& file)
	{
		data_type dat;
		int bits;

		encode_end();

		do {
			bits = encode_process(dat);
			file.write(dat, bits);
		} while (bits > 0);
	}

	template <class gen>
	void encode_end(deque<gen>& p_cont)
	{
		data_type dat;
		int bits;

		encode_end();

		do {
			bits = encode_process(dat);
			if (bits > 0) {
				p_cont.push_back((gen)dat);
			}
		} while (bits > 0);
	}

	void encode_next(bit_file& file, data_type p_data)
	{
		data_type dat;
		int bits;

		encode_next(p_data);

		do {
			bits = encode_process(dat);
			if (bits > 0) {
				file.write(dat, bits);
				//cout << "value: " << (int)dat << " bits: " << bits << endl;
				//getchar();
			}
		} while (bits > 0);
	}

	template <class gen>
	void encode_next(deque<gen>& p_cont, const data_type p_data)
	{
		data_type dat;
		int bits;

		encode_next(p_data);

		do {
			bits = encode_process(dat);
			if (bits > 0) {
				p_cont.push_back((gen)dat);
			}
		} while (bits > 0);
	}

	template<class gen>
	void encode(const deque<gen>& p_from, deque<gen>& p_to)
	{
		encode_start();

		for (unsigned int i = 0; i < p_from.size(); i++) {
			encode_next(p_to, p_from[i]);
		}

		encode_end(p_to);
	}

	/*** Decode methods ***/

	bool is_repeating()
	{
		return span != 0 && run > 0 && span != run;
	}

	// returns 1 if in repete run
	int get_state()
	{
		if (is_repeating()) {
			return 1;
		} else {
			return 0;
		}
	}

	void decode_start()
	{
		run = -1;
		span = 0;
		data_count = 0;
	}

	// returns 0 if done, 1 if output number, -1 if not ready for this number
	int decode_next(data_type& p_data)
	{
		int ret;
		static run_type len;

		if (span == 0) {
			if (run == 0) {
				ret = 0;
			} else {
				bool ctrl = len.write_seg(p_data, data_count);
				data_count++;
				ret = -1;

				if (!ctrl) {
					run = (int)len;
					span = abs(run);
					data_count = 0;

					if (run == 0) ret = 0;
				}
			}
		} else if (run > 0) { // repeating run
			if (span == run) {
				prev = p_data;
			} else {
				p_data = prev;
			}
			ret = 1;
			span--;
		} else {
			ret = 1;
			span--;
		}

		return ret;
	}

	bool decode_next(bit_file& p_file, data_type& p_data)
	{
		if (span == 0) {
			run_type len;
			p_file.read(len);
			run = (int)len;
			span = abs(run);
		}

		if (run > 0) { // repeating run
			if (span == run) {
				p_file.read(p_data, data_bits);
				prev = p_data;
			} else {
				p_data = prev;
			}
		} else if (run < 0) { // literal run
			p_file.read(p_data, data_bits);
		}

		span--;

		return run != 0;
	}

	void decode_end()
	{
	}

	template <class gen>
	void decode(const deque<gen>& p_from, deque<gen>& p_to)
	{
		int res;
		int from_ind = 0;
		int to_ind = 0;

		decode_start();

		do {
			data_type dat;

			if (!is_repeating()) {	
				dat = p_from[from_ind];
				from_ind++;
			}

			res = decode_next(dat);
			
			if (res == 1) {
				p_to.push_back(dat);
				to_ind++;
			}
		} while (res != 0);

		decode_end();
	}
};

#endif