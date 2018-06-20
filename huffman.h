	/*

	Example of writing huffman to a file:

	bit_file file("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\huffman.abc", true);
	huffman<0, 255, 8> huf;

	huf.encode_start();

	huf.encode_next(21);
	huf.encode_next(22);
	huf.encode_next(23);
	huf.encode_next(24);
	huf.encode_next(25);
	for (int i = 0; i < 100; i++)
		huf.encode_next(8);
	huf.encode_next(5);
	huf.encode_next(125);
	for (int i = 0; i < 100; i++)
		huf.encode_next(15);
	huf.encode_next(225);
	huf.encode_next(225);

	huf.encode_end();

	huf.write_trie(file);
	huf.write_data(file);

	file.close();

	getchar();
	return 0;


	Example of reading huffman:

	bit_file file_in("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\huffman.abc");
	huffman<0, 511, 9> huf_in;

	huf_in.read_trie(file_in);
	int len_in = huf_in.decode_start(file_in);
	for (int i = 0; i < len_in; i++) {
		cout << "data: " << (int)huf_in.decode_next(file_in) << endl;
	}

	huf_in.decode_end(file_in);

	getchar();
	return 0;

	*/

#ifndef INC_HUFFMAN
#define INC_HUFFMAN

#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\static_array.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\bit_file.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\run_len.h"
#include <deque>
#include <vector>

template <
	int value_min = 0, int value_max = 255, int value_bits = 8, 
	class value_type = int, class freq_type = unsigned int, bool rle = true>
class huffman {

	class trie {
	public:
		trie* left;
		trie* right;
		freq_type freq;
		value_type value;
		std::vector<bool> code;

		trie ()
		{
		}

		trie (value_type p_value)
		{
			value = p_value;
			left = NULL;
		}

		trie (freq_type p_freq, value_type p_value)
		{
			freq = p_freq;
			value = p_value;
			left = NULL;
		}

		trie (trie* p_left, trie* p_right)
		{
			freq = p_left->freq + p_right->freq;
			value = 0;
			left = p_left;
			right = p_right;
		}

		bool is_node()
		{
			return left == NULL;
		}
	};

	static_array<value_max - value_min + 1, freq_type, value_type> frequencies;
	trie* trie_ptrs[value_max - value_min + 1]; // array of pointers to trie, keyed on value
	trie* trie_root; // beginning of huffman trie
	trie* curr_node;
	std::deque<value_type> values;
	size_t value_count;
	run_len<value_bits, value_bits> run_len_enc;
	
public:

	huffman ()
	{
		trie_root = NULL;
	}

	void encode_start()
	{
		if (rle) {
			run_len_enc.encode_start();
		}

		for (int i = 0; i < value_max - value_min + 1; i++) {
			frequencies.push_back(0, i);
		}
	}

	void encode_next(value_type p_value)
	{
		value_type dat;
		int bits;

		p_value -= value_min;

		if (rle) {
			run_len_enc.encode_next(p_value);
			do {
				bits = run_len_enc.encode_process(dat);

				if (bits > 0) {
					values.push_back(dat);
					frequencies.array[(int)dat].key++;
				}
			} while (bits > 0);
		}
	}

	trie* create_lesser(static_array<value_max - value_min + 1, freq_type, value_type>& p_nodes, std::deque<trie*>& p_queue)
	{
		trie* next;

		if (p_queue.empty()) {
			next = new trie(p_nodes.front().key, p_nodes.front().data);
			p_nodes.pop_front();
		} else {
			next = p_queue.front();

			if (p_nodes.empty() || next->freq < p_nodes.front().key) {
				p_queue.pop_front();
			} else {
				next = new trie(p_nodes.front().key, p_nodes.front().data);
				p_nodes.pop_front();
			}
		}

		return next;
	}

	void assign_binary(trie* p_node)
	{
		std::vector<bool> p_code;
		assign_binary(p_node, p_code);
	}

	void assign_binary(trie* p_node, std::vector<bool>& p_code)
	{
		if (p_node->is_node()) { // is a freq
			p_node->code = p_code;
			trie_ptrs[p_node->value] = p_node;
		} else { // is a sum
			std::vector<bool> lef = p_code;
			std::vector<bool> rig = p_code;
			lef.push_back(false); // left turn is 0
			rig.push_back(true);  // right turn is 1
			assign_binary(p_node->left, lef);
			assign_binary(p_node->right, rig);
		}
	}

	void encode_end()
	{
		trie* a = NULL;
		trie* b = NULL;
		std::deque<trie*> queue;
		value_type dat;
		int bits;

		// encode the last few items since it's been rle
		
		run_len_enc.encode_end();

		do {
			bits = run_len_enc.encode_process(dat);
			values.push_back(dat);
			frequencies.array[(int)dat].key++;
		} while (bits > 0);

		// sort histogram of unique values by frequency
		frequencies.sort_array();
		
		// skip forward past all of the zero frequencies
		while (!frequencies.empty() && frequencies.front().key == 0) {
			frequencies.pop_front();
		}

		// build huffman tree
		while (!frequencies.empty() || queue.size() > 1) {
			a = create_lesser(frequencies, queue);
			b = create_lesser(frequencies, queue);
			queue.push_back(new trie(a, b));
		}

		trie_root = queue.front();
		
		assign_binary(trie_root);
	}

	void clear_trie(const trie* p_start)
	{
		if (p_start->left != NULL) {
			clear_trie(p_start->left);
			clear_trie(p_start->right);
		}

		delete p_start;
	}

	void clear_trie()
	{
		if (trie_root != NULL) {
			clear_trie(trie_root);
			trie_root = NULL;
		}
	}

	void send()
	{
		for (unsigned i = 0; i < values.size(); i++) {
			std::vector<bool> n = trie_ptrs[values[i]]->code;
			for (unsigned j = 0; j < n.size(); j++) {
				cout << n[j];
			}
			cout << ",";
		}
	}

	size_t coded_length()
	{
		size_t total = 0;

		for (unsigned i = 0; i < values.size(); i++) {
			std::vector<bool> n = trie_ptrs[values[i]]->code;
			total += n.size();
		}

		return total;
	}

	float calc_ratio()
	{
		return float(values.size() * 8) / (float)coded_length();
	}

	void write_trie(bit_file& p_file)
	{
		write_trie(p_file, trie_root);
	}

	void write_trie(bit_file& p_file, trie* p_node)
	{
		if (p_node->is_node()) { // is an actual value (so no left or right branches)
			p_file.write_bit(true);	
			p_file.write(p_node->value, value_bits);
		} else { // is a connector so explore left and right branches
			p_file.write_bit(false);
			write_trie(p_file, p_node->left);
			write_trie(p_file, p_node->right);
		}
	}

	void read_trie(bit_file& p_file)
	{
		clear_trie();
		trie_root = read_trie_node(p_file);
	}

	trie* read_trie_node(bit_file& p_file)
	{
		trie* ret;

		bool is_value = p_file.read_bit();
		if (is_value) {
			value_type value;
			p_file.read(value, value_bits);
			ret = new trie(value);
		} else {
			trie* lef = read_trie_node(p_file);
			trie* rig = read_trie_node(p_file);
			ret = new trie(lef, rig);
		}

		return ret;
	}

	void write_data(bit_file& p_file)
	{
		p_file.write(values.size());

		for (size_t i = 0; i < values.size(); ++i) {
			trie* t = trie_ptrs[values[i]];
			for (size_t j = 0; j < t->code.size(); ++j) {
				bool b = t->code[j];
				p_file.write_bit(b);
				//cout << b;
			}
			//cout << endl;
		}
	}

	// initializes huffman decoding and returns number of encoded values
	size_t decode_start(bit_file& p_file)
	{
		p_file.read(value_count);
		curr_node = trie_root;
		
		if (rle) {
			run_len_enc.decode_start();
		}

		return value_count;
	}

	void decode_end(bit_file& p_file)
	{
		if (rle) {
			run_len_enc.decode_end();
		}
	}

	bool decode_next(bit_file& p_file, value_type& p_data)
	{
		bool ret;

		if (rle) {
			int rl;

			do {
				if (run_len_enc.get_state() != 1) {	
					while (!curr_node->is_node()) {
						if (p_file.read_bit()) {
							curr_node = curr_node->right; // right is 1
						} else {
							curr_node = curr_node->left; // left is 0
						}
					};
				}

				p_data = curr_node->value;
				curr_node = trie_root;

				rl = run_len_enc.decode_next(p_data);
			} while (rl == -1);

			ret = (rl != 0);
		} else {
			// todo: huffman without RLE
		}

		return ret;
	}
};

#endif