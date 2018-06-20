#ifndef INC_RUN_LEN_2D
#define INC_RUN_LEN_2D

#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\bit_file.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\varint.h"
#include <deque>

// min_span >= 2 * sizeof(run_type) + 1

// birds: min_run = 2: 641 kb
// birds: min_run = 3: 716 kb

template <class data_type = int, class run_type = int, class pos_type = int, data_type init_data = -1, int min_run = 4>
class run_len_2d {

	int width, height;
	int segment_index;
	pos_type pos;
	int data_ct;

	class span {
	public:
		pos_type start;
		run_type len;

		span ()
		{
		}

		span (pos_type p_start, run_type p_len)
		{
			start = p_start;
			len = p_len;
		}
	};

	class run_tracer {
		data_type last_data;
		pos_type pos;
		run_type run_len;

	public:

		deque<span> segments;

		run_tracer ()
		{
			init();
		}

		void init ()		
		{
			last_data = init_data;
			pos = 0;
			run_len = 0;
			segments.clear();
		}

		pos_type get_pos()
		{
			return pos;
		}

		bool is_span ()
		{
			return run_len + 1 >= min_run;
		}

		bool update (data_type p_data)
		{
			bool ret;

			if (p_data == last_data) {
				ret = false;
				run_len++;
			} else {
				ret = is_span();
				if (ret) cut();
				run_len = 0;
			}

			pos++;
			last_data = p_data;

			return ret;
		}

		void cut ()
		{
			segments.push_back(span(pos - run_len - 1, run_len + 1));
		}

		void cleanup_last ()
		{
			if (is_span()) {
				cut();
			}
		}
	};

	run_tracer horz;

public:

	deque<data_type> data_stream;

	run_len_2d (int p_width, int p_height)
	{
		init(p_width, p_height);
	}

	run_len_2d ()
	{
	}

	void init (int p_width, int p_height)
	{
		width = p_width;
		height = p_height;
	}

	void encode_start()
	{
		horz.init();
	}

	bool encode_next(data_type p_data)
	{
		data_stream.push_back(p_data);
		horz.update(p_data);
		return horz.get_pos() < width * height;
	}

	void encode_end()
	{
		horz.cleanup_last();

		for (int i = 0; i < horz.segments.size(); ++i) {
			//cout << "pos: " << horz.segments[i].start << " len: " << horz.segments[i].len << endl;
			//getchar();
		}
	}

	void encode_process_meta_start()
	{
		segment_index = 0;
		pos = 0;
	}

	bool encode_process_meta_next(run_type& p_meta)
	{
		if (segment_index > horz.segments.size()) {
			return false;
		} else if (segment_index == horz.segments.size()) {
			if (pos >= width * height - 1) { // last segment was span (run)
				return false;
			} else { // last segment was skip
				p_meta = -(width * height - pos);
				segment_index++;
				return true;
			} 
		} else {
			pos_type s = horz.segments[segment_index].start;
			run_type l = horz.segments[segment_index].len;

			if (pos < s) {
				p_meta = -(s - pos); // skip is negative
				pos = s;
			} else {
				p_meta = l; // run is positive
				pos += l;
				segment_index++;
			}

			return true;
		}
	}

	void encode_process_meta_end()
	{
	}

	int process_encode(bit_file& f)
	{
		run_type meta;
		data_type data;
		int r = 0;
		
		encode_process_meta_start();

		while (encode_process_meta_next(meta)) {
			//cout << "meta: " << meta << endl;
			//getchar();
			//f.write((varint<5, true>)meta); // 157 kb
			//f.write((varint<8, true>)meta); // 157 kb
			//f.write((varint<7, true>)meta); // 157 kb
			//f.write((varint<16, true>)meta); //  155 kb
			//f.write((varint<12, true>)meta); // 155 kb
			f.write((varint<8, true, int, 0, 1>)meta); // 157 kb

			//f.write(meta); //  kb
			r += abs(meta);
		}

		encode_process_meta_end();

		encode_process_data_start();

		while (encode_process_data_next(data)) {
			//cout << "pos: " << pos << " data: " << data << endl;
			//getchar();
			f.write((unsigned char)data);
		}

		return r;
	}

	void encode_process_data_start()
	{
		segment_index = 0;
		pos = 0;
	}

	bool encode_process_data_next(data_type& p_data)
	{
		pos_type s;
		run_type l;
		pos_type total_len = width * height;

		if (pos >= total_len) return false;

		if (segment_index >= horz.segments.size()) {
			s = total_len;
			l = s - pos;
		} else {
			s = horz.segments[segment_index].start;
			l = horz.segments[segment_index].len;
		}


		//cout << "pos: " << pos << " seg ind: " << segment_index << " start: " << s << " len: " << l << endl;
		//getchar();

		p_data = data_stream[pos];

		if (pos < s) {
			pos++;
		} else {
			pos += l;
			segment_index++;
		}

		return true;
	}

	void encode_process_data_end()
	{
	}

	void decode_meta_start(int p_width, int p_height)
	{
		init(p_width, p_height);

		segment_index = 0;
		pos = 0;
		horz.segments.clear();
	}

	bool decode_meta_next(const run_type p_meta)
	{
		run_type len;

		if (p_meta > 0) {
			//cout << "span: " << p_meta << " pos: " << pos << endl; getchar();
			len = (pos_type)p_meta;
			horz.segments.push_back(span(pos, len));
			data_ct++;
			segment_index++;
		} else {
			//cout << "run: " << p_meta << " pos: " << pos << endl; getchar();
			len = -(pos_type)p_meta;
			data_ct += len;
		}

		pos += (pos_type)len;

		return pos < width * height;
	}

	int decode_meta_end()
	{
		for (int i = 0; i < horz.segments.size(); ++i) {
			//cout << "meta decode: " << horz.segments[i] << " i: " << i << endl;
		}

		return data_ct;
	}

	void decode_data_start()
	{
		segment_index = 0;
		pos = 0;
		data_stream.clear();
	}

	bool decode_data_next(const data_type p_data)
	{
		//cout << "segs: " << horz.segments.size() << endl; getchar();
		pos_type s;
		run_type l;

		if (segment_index < horz.segments.size()) {
			s = horz.segments[segment_index].start;
			l = horz.segments[segment_index].len;
		} else {
			s = width * height;
			l = s - pos;
		}

		if (pos < s) {
			//cout << "skip: " << endl; getchar();
			data_stream.push_back(p_data);
			pos++;
		} else {
			//cout << "span: primary: " << p_data << " len: " << l << endl; getchar();
			for (int i = 0; i < l; ++i) {
				data_stream.push_back(p_data);
			}
			segment_index++;
			pos += l;
		}

		return pos < width * height;
	}

	void decode_data_end()
	{
	}

	void decode_process(bit_file& f, int p_width, int p_height)
	{
		varint<8, true, int, 0, 1> meta;
		//int meta;
		unsigned char data;

		decode_meta_start(p_width, p_height);

		do {
			f.read(meta);
			//cout << "meta: " << (int)meta << endl; getchar();
		} while (decode_meta_next((run_type)meta));

		decode_meta_end();

		decode_data_start();

		do {
			f.read(data);
		} while (decode_data_next((data_type)data));

		decode_data_end();
	}
};

#endif // INC_RUN_LEN_2D