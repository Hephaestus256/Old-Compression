#ifndef INC_BITMAP
#define INC_BITMAP

#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_1d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\math_2d.h"
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\color.h"

class bitmap {

public:

	rgb* data;
	quantum_2d dims;
	int width_bytes;

	struct {
		unsigned short int type;             /* Magic identifier            */
		unsigned int size;                       /* File size in bytes          */
		unsigned short int reserved1, reserved2;
		unsigned int offset;                   /* Offset to image data, bytes */
	} bmp_header;

	struct {
		unsigned int size;               /* Header size in bytes      */
		int width, height;                /* Width and height of image */
		unsigned short int planes;       /* Number of colour planes   */
		unsigned short int bits;         /* Bits per pixel            */
		unsigned int compression;        /* Compression type          */
		unsigned int imagesize;          /* Image size in bytes       */
		int xresolution, yresolution;     /* Pixels per meter          */
		unsigned int ncolors;           /* Number of colours         */
		unsigned int importantcolors;   /* Important colours         */
	} bmp_info;

	bitmap (const string p_filename)
	{
		load_from_file(p_filename);
	}

	bitmap (int x_res, int y_res, rgb p_color = RGB_BLACK)
	{
		reserve(x_res, y_res);
		paint_color(p_color);
	}

	int calc_width_pad()
	{
		return (4 - ((dims.x * 3) % 4)) % 4;
	}

	int calc_padded_width()
	{
		return dims.x * 3 + calc_width_pad();
	}

	int calc_image_size()
	{
		return width_bytes * dims.y;
	}

	void reserve(int x_res, int y_res)
	{
		dims.x = x_res;
		dims.y = y_res;
		width_bytes = calc_padded_width();
		int size = calc_image_size();
		data = new rgb[size]; // allocate 3 bytes per pixel
	}

	void paint_color(rgb p_color = RGB_BLACK)
	{
		for (int y = 0; y < dims.y; y++) {
			for (int x = 0; x < dims.x; x++) {
				write_rgb(x, y, p_color);
			}
		}
	}

	bool load_from_file(const string p_filename)
	{
		FILE* f;
		errno_t res = fopen_s(&f, p_filename.c_str(), "rb");
		fread((unsigned char*)&bmp_header.type, 2, 1, f);

		if (bmp_header.type != 19778) return false;

		fread((unsigned char*)&bmp_header.size, 4, 1, f);
		fread((unsigned char*)&bmp_header.reserved1, 2, 1, f);
		fread((unsigned char*)&bmp_header.reserved2, 2, 1, f);
		fread((unsigned char*)&bmp_header.offset, 4, 1, f);

		fread((unsigned char*)&bmp_info, 40, 1, f);

		// extract image height and width from header
		reserve(bmp_info.width, bmp_info.height);
		fread(data, sizeof(rgb), calc_image_size(), f); // read the rest of the data at once
		fclose(f);

		return true;
	}

	void save_to_file(const string p_filename)
	{
		bmp_header.type = 19778;
		bmp_header.size = 54 + dims.x * dims.y;
		bmp_header.reserved1 = 0;
		bmp_header.reserved2 = 0;
		bmp_header.offset = 54;

		bmp_info.size = 40;
		bmp_info.width = dims.x;
		bmp_info.height = dims.y;
		bmp_info.planes = 1;
		bmp_info.bits = 24;
		bmp_info.compression = 0;
		bmp_info.imagesize = calc_image_size();
		bmp_info.xresolution = 0;
		bmp_info.yresolution = 0;
		bmp_info.ncolors = 0;
		bmp_info.importantcolors = 0;

		std::fstream file;

		int mode = ios::in | ios::out | ios::binary;
		mode |= ios::trunc;
		file.open(p_filename, mode);

		file.write((const char*)&bmp_header.type, 2);
		file.write((const char*)&bmp_header.size, 4);
		file.write((const char*)&bmp_header.reserved1, 2);
		file.write((const char*)&bmp_header.reserved2, 2);
		file.write((const char*)&bmp_header.offset, 4);

		file.write((const char*)&bmp_info, sizeof(bmp_info));
		file.write((const char*)data, calc_image_size());

		file.close();
	}

	template <class gen>
	void load_from_deques (
		deque<gen> p_data_r, 
		deque<gen> p_data_g, 
		deque<gen> p_data_b,
		int p_width, int p_height)
	{
		int x, y;
		int len = p_width * p_height;
		reserve(p_width, p_height);

		x = 0;
		y = 0;
		for (int i = len - 1; i >= 0 ; --i) {
			rgb c = RGB_BLACK;
			if ((int)p_data_r.size() - 1 >= i) c.rgb_r = p_data_r[i];
			if ((int)p_data_g.size() - 1 >= i) c.rgb_g = p_data_g[i];
			if ((int)p_data_b.size() - 1 >= i) c.rgb_b = p_data_b[i];

			write_rgb(x, y, c);

			x++;
			if (x >= p_width) {
				x = 0;
				y++;
			}
		}
	}

	rgb read_rgb (int p_x, int p_y)
	{
		char* p = ((char*)data + p_x * 3 + p_y * width_bytes);
		return ((rgb*)p)->reverse();
	}

	void write_rgb (int p_x, int p_y, rgb p_color)
	{
		char* p = (char*)data + (dims.x - 1 - p_x) * 3 + (dims.y - 1 - p_y) * width_bytes;
		*((rgb*)p) = p_color.reverse();
	}

	/*
		this does NOT work, IDK why
	void write_rgb_channel (int p_x, int p_y, int p_value, int p_channel)
	{
		char* p = (char*)data + (dims.x - 1 - p_x) * 3 + (dims.y - 1 - p_y) * width_bytes;
		*(p + (2 - p_channel)) = (char)p_value;
	}
	*/

	void write_primary (int p_x, int p_y, int p_chan, int p_prim)
	{
		char* p = ((char*)data + ((dims.x - 1 - p_x) * 3 + (dims.y - 1 - p_y) * width_bytes));
		((rgb*)p)->channel[2 - p_chan] = p_prim;
	}

	void compare(bitmap& p_bmp)
	{
		bool ret = true;

		quantum_2d r_fault = -1;
		quantum_2d g_fault = -1;
		quantum_2d b_fault = -1;

		if (dims.x != p_bmp.dims.x) {
			cout << "X dim different" << endl;
			ret = false;
		}

		if (dims.y != p_bmp.dims.y) {
			cout << "Y dim different" << endl;
			ret = false;
		}

		if (bmp_info.bits != p_bmp.bmp_info.bits) {
			cout << "Bits per pixel different" << endl;
			ret = false;
		}

		if (ret) {
			pixel<int> diff = 0;
			pixel<int> d;

			for (int y = 0; y < dims.y; y++) {
				for (int x = 0; x < dims.x; x++) {
					d = (read_rgb(x, y) - p_bmp.read_rgb(x, y));
					diff += d;

					if (r_fault.x < 0 && d.rgb_r != 0) r_fault = quantum_2d(x, y);
					if (g_fault.x < 0 && d.rgb_g != 0) g_fault = quantum_2d(x, y);
					if (b_fault.x < 0 && d.rgb_b != 0) b_fault = quantum_2d(x, y);
				}
			}

			cout << "Delta R: " << diff.rgb_r << "(" << (double)diff.rgb_r / (dims.x * dims.y) << " per pixel) fault: " << r_fault.report() << endl;
			cout << "Delta G: " << diff.rgb_g << "(" << (double)diff.rgb_g / (dims.x * dims.y) << " per pixel) fault: " << g_fault.report() << endl;
			cout << "Delta B: " << diff.rgb_b << "(" << (double)diff.rgb_b / (dims.x * dims.y) << " per pixel) fault: " << b_fault.report() << endl;
		}
	}
};

#endif INC_BITMAP