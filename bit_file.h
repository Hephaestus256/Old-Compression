#ifndef INC_BIT_FILE
#define INC_BIT_FILE

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\varint.h"

using namespace std;

class bit_file {
	
	unsigned int bits;

	void write_file()
	{
		file.write((const char*)&buffer, sizeof(buffer));
		buffer = 0;
	}

public:
	std::fstream file;
	unsigned char buffer;

	void open (const string p_file, bool p_truc = false)
	{
		int mode = ios::in | ios::out | ios::binary;
		if (p_truc) mode |= ios::trunc;
		file.open(p_file, mode);
	}

	void close()
	{
		write_file();
		file.close();
	}

	bit_file (const string p_file, bool p_trunc = false)
	{
		bits = 0;
		buffer = 0;
		open(p_file, p_trunc);
	}

	~bit_file()
	{
		close();
	}

	/*
	void write_file(unsigned int p_data, unsigned int p_bits)
	{
		int bytes = p_bits / 8;

		file.write((const char*)&p_data, bytes);
		for (unsigned int i = bytes * 8 + 1; i <= p_bits; i++) {
			write((p_data & (1 << i)) != 0);
		}
	}
	*/

	bool read_bit()
	{
		int buf_bits = sizeof(buffer) * 8;
		int buf_mask = 1 << ((buf_bits - 1) - (bits % buf_bits));

		if (bits % buf_bits == 0) {
			file.read((char*)&buffer, sizeof(buffer));
		}
		bits++;

		return (buffer & buf_mask) != 0;
	}

	template <int seg_len, bool signd, class cont_type, int zero, int xtra_n>
	void read(varint<seg_len, signd, cont_type, zero, xtra_n>& p_data)
	{
		cont_type data;
		int i = 0;
		bool cont;

		p_data.clear();

		do {
			read(data, seg_len);
			cont = p_data.write_seg(data, i);
			i++;
		} while (cont);
	}

	template <class gen>
	void read(gen& p_data, unsigned int p_bits = sizeof(gen) * 8)
	{
		p_data = 0;
		
		for (int i = p_bits - 1; i >= 0; i--) {
			bool d = read_bit();
			p_data = p_data | ((gen)d << (gen)i);
		}
	}

	void read(string& p_data)
	{
		char c;

		while (file.read(&c, 1)) {
			p_data += c;
		}
	}

	void read_unary(int& n, int base = 0)
	{
		n = base;

		while (read_bit()) {
			n++;
		}
	}

	void write_bit(bool p_data)
	{
		unsigned int buf_bits = sizeof(buffer) * 8;

		if (bits > 0 && bits % buf_bits == 0) {
			write_file();
		}

		buffer |= (p_data << ((buf_bits - 1) - (bits & (buf_bits - 1)))); 
		bits++;
	}

	template <int seg_len, bool signd, class cont_type, int zero, int xtra_n>
	void write(varint<seg_len, signd, cont_type, zero, xtra_n> p_data)
	{
		int i = 0;

		do {
			write(p_data.read_seg(i), seg_len);
		} while (p_data.read_seg_ctrl(i++));
	}

	template <class gen>
	void write(gen p_data, int bit_len = sizeof(gen) * 8)
	{
		for (int i = bit_len - 1; i >= 0; i--) {
			bool d = bool((p_data & (1 << i)) != 0);
			write_bit(d);
		}
	}

	void write_unary(int n, int base = 0)
	{
		n -= base;

		for (int i = 0; i < n - 1; ++i) {
			write_bit(true);
		}
		write_bit(false);
	}
};

#endif