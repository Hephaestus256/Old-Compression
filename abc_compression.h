#ifndef INC_ABC_COMPRESSION
#define INC_ABC_COMPRESSION

#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_1d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_2d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\huffman.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\bitmap.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\fixd.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\bit_file.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\varint.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\run_len_2d.h"

#include <string>

class abc_compression {
	int width, height;
public:
	abc_compression()
	{
	}

	~abc_compression()
	{
	}

	void write_header(bit_file& p_file, int p_width = -1, int p_height = -1)
	{
		if (p_width != -1) {
			width = p_width;
			height = p_height;
		}

		p_file.write(width);
		p_file.write(height);
	}

	void read_header(bit_file& p_file)
	{
		p_file.read(width);
		p_file.read(height);
	}

	void compress_from_bmp(string p_bmp_file, string p_abc_file)
	{
		bitmap bmp(p_bmp_file);
		bit_file abc_out(p_abc_file, true);
		//huffman<0, 255, 8, int, unsigned int, true> huff;
		run_len<8, 8> rle;

		write_header(abc_out, bmp.dims.x, bmp.dims.y);		
		
		// encode abc data

		//huff.encode_start();
		rle.encode_start();
		
		for (int i = 0; i < 3; i++) {
			for (int y = bmp.dims.y - 1; y >= 0; y--) {
				for (int x = bmp.dims.x - 1; x >= 0; x--) {
					rgb a, b, c, d;

					if (x > 0) {
						a = bmp.read_rgb(x - 1, y);
						if (y > 0) {
							b = bmp.read_rgb(x - 1, y - 1);
						} else {
							b = 0;
						}
					} else {
						a = 0;
						b = 0;
					}

					if (y > 0) {
						c = bmp.read_rgb(x, y - 1);
						//if (x < 
						//rgb d = bmp.read_rgb(x + 1, y - 1);
					} else {
						c = 0;
					}
					d = 0;
					pixel<int> mean = (a + b + c + d) * 0.3333333333333333333;

					pixel<int> curr = bmp.read_rgb(x, y).cv_int() - mean;
					rle.encode_next(abc_out, (int)curr.channel[i]);
					//huff.encode_next((int)c.channel[i]);
				}
			}
		}

		rle.encode_end(abc_out);
		//huff.encode_end();
		//huff.write_trie(abc_out);
		//huff.write_data(abc_out);

		abc_out.close();
	}
		
	void decompress_to_bmp(const string p_abc_file, const string p_bmp_file)
	{
		cout << "abc: " << p_abc_file << " bmp: " << p_bmp_file << endl;
		getchar();

		run_len<> rle;
		//huffman<0, 255, 8, int, unsigned int, true> huff;

		bit_file abc_in(p_abc_file);
		int dat;
		int x, y;

		read_header(abc_in);
		bitmap bmp(width, height);

		cout << "width: " << width << endl;
		cout << "height: " << height << endl;
		getchar();

		rle.decode_start();
		//huff.read_trie(abc_in);
		//huff.decode_start(abc_in);

		for (int i = 0; i < 3; i++) {
			x = 0;
			y = 0;
			while (y < height) { // && rle.decode_next(abc_in, dat)) {
				rle.decode_next(abc_in, dat);
				//huff.decode_next(abc_in, dat);
				bmp.write_primary(x, y, i, dat);
				x++;

				if (x == width) {
					x = 0;
					y++;
					//cout << "y: " << y << endl;
					//getchar();
				}
			};
		}

		rle.decode_end();
		//huff.decode_end(abc_in);
		abc_in.close();

		bmp.save_to_file(p_bmp_file);
	}

	void test_compression(const string dir, const string file_in, const string file_out)
	{
		abc_compression abc;
		abc_compression abc2;

		string a = dir + "\\" + file_in + ".bmp";
		string b = dir + "\\" + file_in + ".abc";
		string c = dir + "\\" + file_out + ".bmp";

		abc.compress_from_bmp(a, b);
		abc2.decompress_to_bmp(b, c);

		bitmap a1(a);
		bitmap a2(c);

		a1.compare(a2);
		getchar();
	}
};

#endif