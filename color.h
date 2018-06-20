#ifndef INC_COLOR
#define INC_COLOR

#include <math.h>
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\fixd.h"

char color_h_report_text[50];

#define RGB_BLACK rgb(0, 0, 0)
#define RGB_WHITE rgb(255, 255, 255)
#define RGB_GREY rgb(127, 127, 127)
#define RGB_RED rgb(255, 0, 0)
#define RGB_GREEN rgb(0, 255, 0)
#define RGB_BLUE rgb(0, 0, 255)
#define RGB_YELLOW rgb(255, 255, 0)
#define RGB_AQUA rgb(0, 255, 255)
#define RGB_PURPLE rgb(255, 0, 255)
#define RGB_ORANGE rgb(255, 165, 0)
#define COLOR_PERCEPT pixel<float>(0.2989, 0.5870, 0.1140)
#define RGB_MAX pixel<float>(255.0, 255.0, 255.0)

template <
	class data_type1 = int, class data_type2 = data_type1, class data_type3 = data_type2,
	class array_type = data_type1>
class pixel {
public:
	union {
		struct {
			data_type1 rgb_r;
			data_type2 rgb_g;
			data_type3 rgb_b;
		};
		struct {
			data_type1 xyz_x;
			data_type2 xyz_y;
			data_type3 xyz_z;
		};
		struct {
			data_type1 lab_l;
			data_type2 lab_a;
			data_type3 lab_b;
		};
		struct {
			data_type1 lch_l;
			data_type2 lch_c;
			data_type3 lch_h;
		};
		struct {
			data_type1 ycbcr_y;
			data_type2 ycbcr_cb;
			data_type3 ycbcr_cr;
		};
		struct {
			data_type1 ycocg_y;
			data_type2 ycocg_co;
			data_type3 ycocg_cg;
		};
		struct {
			array_type channel[3];
		};
	};

	pixel ()
	{
	}

	template <class gen>
	pixel (gen lightness)
	{
		rgb_r = (data_type1)lightness;
		rgb_g = (data_type2)lightness;
		rgb_b = (data_type3)lightness;
	}

	template <class gen1, class gen2, class gen3, class gen4>
	pixel (pixel<gen1, gen2, gen3, gen4> p)
	{
		rgb_r = (data_type1)p.rgb_r;
		rgb_g = (data_type2)p.rgb_g;
		rgb_b = (data_type3)p.rgb_b;
	}

	template <class gen1, class gen2, class gen3>
	pixel (gen1 p_r, gen2 p_g, gen3 p_b)
	{
		rgb_r = (data_type1)p_r;
		rgb_g = (data_type2)p_g;
		rgb_b = (data_type3)p_b;
	}

	bool operator == (pixel c)
	{
		return rgb_r == c.rgb_r && rgb_g == c.rgb_g && rgb_b == c.rgb_b;
	}

	bool operator != (pixel c)
	{
		return rgb_r != c.rgb_r || rgb_g != c.rgb_g || rgb_b != c.rgb_b;
	}

	bool operator <= (pixel c)
	{
		return rgb_r <= c.rgb_r && rgb_g <= c.rgb_g && rgb_b <= c.rgb_b;
	}

	bool operator < (pixel c)
	{
		return rgb_r < c.rgb_r && rgb_g < c.rgb_g && rgb_b < c.rgb_b;
	}

	bool operator >= (pixel c)
	{
		return rgb_r >= c.rgb_r && rgb_g >= c.rgb_g && rgb_b >= c.rgb_b;
	}

	bool operator > (pixel c)
	{
		return rgb_r > c.rgb_r && rgb_g > c.rgb_g && rgb_b > c.rgb_b;
	}

	pixel<int> cv_int()
	{
		return pixel<int>(
			int(rgb_r), 
			int(rgb_g), 
			int(rgb_b)
		);
	}

	pixel<int> cv_int_round()
	{
		return pixel<int>(
			int(rgb_r + data_type1(0.499)), 
			int(rgb_g + data_type2(0.499)), 
			int(rgb_b + data_type3(0.499))
		);
	}

	pixel<unsigned char> cv_byte()
	{
		return pixel<unsigned char>(
			unsigned char(cap_range((unsigned char)rgb_r, (unsigned char)0, (unsigned char)255)), 
			unsigned char(cap_range((unsigned char)rgb_g, (unsigned char)0, (unsigned char)255)), 
			unsigned char(cap_range((unsigned char)rgb_b, (unsigned char)0, (unsigned char)255))
		);
	}

	pixel<float> cv_float()
	{
		return pixel<float>(float(rgb_r), float(rgb_g), float(rgb_b));
	}

	pixel operator * (pixel<float> p)
	{
		return pixel(data_type1(rgb_r * p.rgb_r), data_type2(rgb_g * p.g), data_type3(rgb_b * p.b));
	}
	
	pixel reverse()
	{
		return pixel(rgb_b, rgb_g, rgb_r);
	}

	pixel<float> denormalize(float p_max0 = 255, float p_max1 = 255, float p_max2 = 255)
	{
		return pixel<float>(
			float(channel[0]) * float(p_max0), 
			float(channel[1]) * float(p_max1), 
			float(channel[2]) * float(p_max2)
		); 
	}

	pixel<float> normalize(float p_max0 = 255, float p_max1 = 255, float p_max2 = 255)
	{
		//return denormalize(1.0f / p_max1, 1.0f / p_max2, 1.0f / p_max3);
		return pixel<float>(channel[0] / p_max0, channel[1] / p_max1, channel[2] / p_max2);
	}

	template <class gen>
	pixel gamma(gen gamma)
	{
		pixel<float> norm = normalize();
		pixel<float> p(
			(float)pow((float)norm.rgb_r, (float)gamma),
			(float)pow((float)norm.rgb_g, (float)gamma),
			(float)pow((float)norm.rgb_b, (float)gamma)
		);
		
		return pixel(p.denormalize());
	}

	template <class gen>
	pixel<float> gamma_norm(gen gamma)
	{
		pixel<float> norm = normalize();
		pixel<float> p(
			(float)pow((float)norm.rgb_r, (float)gamma),
			(float)pow((float)norm.rgb_g, (float)gamma),
			(float)pow((float)norm.rgb_b, (float)gamma)
		);
		
		return p;
	}

	int sum()
	{
		return rgb_r + rgb_g + rgb_b;
	}

	int abs_sum()
	{
		return abs(rgb_r) + abs(rgb_g) + abs(rgb_b);
	}

	// 2PI >= a >= 0.0
	static pixel show_angle(double a)
	{
		double f = (a / (2.0 * PI));

		if (is_nan(a)) {
			return RGB_BLACK;
		} else if (f < 0.25) {
			return (rgb)interpolate(f, 0.0, 0.25, RGB_RED.cv_float(), RGB_YELLOW.cv_float());
		} else if (f < 0.50) {
			return (rgb)interpolate(f - 0.25, 0.0, 0.25, RGB_YELLOW.cv_float(), RGB_GREEN.cv_float());
		} else if (f < 0.75) {
			return (rgb)interpolate(f - 0.50, 0.0, 0.25, RGB_GREEN.cv_float(), RGB_BLUE.cv_float());
		} else {
			return (rgb)interpolate(f - 0.75, 0.0, 0.25, RGB_BLUE.cv_float(), RGB_RED.cv_float());
		}
	}

	// euclidean distance from color space origin
	float dist_euc()
	{
		return
			sqrt(
				float(rgb_r) * float(rgb_r) 
				+ float(rgb_g) * float(rgb_g) 
				+ float(rgb_b) * float(rgb_b)
			);
	}

	float dist_euc_sqrd()
	{
		return
			float(rgb_r) * float(rgb_r) 
			+ float(rgb_g) * float(rgb_g) 
			+ float(rgb_b) * float(rgb_b)
		;
	}

	// euclidean distance from 
	float dist_euc(pixel c)
	{
		return (cv_int() - c.cv_int()).dist_euc();
	}

	float dist_euc_sqrd(pixel c)
	{
		return (*this - c).dist_euc_sqrd();
	}

	bool euc_dist_thresh(pixel c, float thresh)
	{
		return dist_euc_sqrd(c) <= thresh * thresh;
	}

	// color perception 1: rgb weighted
	pixel color_percept1()
	{
		return pixel(*this * COLOR_PERCEPT);
	}

	// color perception 2: rgb weighted and gamma adjusted
	pixel color_percept2()
	{
		return color_percept1().gamma(1.0 / 2.4);
	}

	pixel color_percept3()
	{
		return rgb_to_xyz();
	}

	pixel color_percept4()
	{
		return rgb_to_lab();
	}

	pixel color_percept5()
	{
		return rgb_to_lch();
	}

	pixel color_percept6()
	{
		return rgb_to_lab() * pixel<float>(0.5, 1.0, 1.0);
	}

	pixel false_color()
	{
		return pixel(
			(rgb_r * 321 + rgb_g * 789 + rgb_b * 555) & 255,
			(rgb_r * 999 + rgb_g * 876 + rgb_b * 444) & 255,
			(rgb_r * 333 + rgb_g * 432 + rgb_b * 888) & 255
		);
	}

	template <class gen>
	pixel operator * (gen n)
	{
		return pixel(data_type1(rgb_r * n), data_type2(rgb_g * n), data_type3(rgb_b * n));
	}

	pixel operator + (pixel c)
	{
		return pixel(rgb_r + c.rgb_r, rgb_g + c.rgb_g, rgb_b + c.rgb_b);
	}

	pixel operator - (pixel c)
	{
		return pixel(rgb_r - c.rgb_r, rgb_g - c.rgb_g, rgb_b - c.rgb_b);
	}

	pixel operator += (pixel c)
	{
		return pixel(rgb_r += c.rgb_r, rgb_g += c.rgb_g, rgb_b += c.rgb_b);
	}

	pixel operator -= (pixel c)
	{
		return pixel(rgb_r -= c.rgb_r, rgb_g -= c.rgb_g, rgb_b -= c.rgb_b);
	}

	pixel operator - ()
	{
		return pixel(-rgb_r, -rgb_g, -rgb_b);
	}

	pixel min_value(pixel c)
	{
		return pixel(min(r, c.rgb_r), min(rgb_g, c.rgb_g), min(rgb_b, c.rgb_b));
	}

	pixel max_value(pixel c)
	{
		return pixel(max(r, c.rgb_r), max(g, c.rgb_g), max(b, c.rgb_b));
	}

	pixel alpha_blend(pixel c, double alpha = 0.5)
	{
		return pixel(interpolate(alpha, rgb_r, c.rgb_r), interpolate(alpha, rgb_g, c.rgb_g), interpolate(alpha, rgb_b, c.rgb_b));
	}

	pixel jiggle(int n)
	{
		// 0, -1, +1
		data_type1 j1 = n % 3 - 1;
		data_type2 j2 = (n / 3) % 3 - 1;
		data_type3 j3 = (n / 9) % 3 - 1;

		return pixel(rgb_r + j3, rgb_g + j2, rgb_b + j1);
	}

	pixel<float> rgb_to_xyz()
	{
		double var_r = (rgb_r / 255.0);
		double var_g = (rgb_g / 255.0);
		double var_b = (rgb_b / 255.0);

		if (var_r > 0.04045) {
			var_r = pow(((var_r + 0.055) / 1.055), 2.4);
		} else {                   
			var_r = var_r / 12.92;
		}

		if (var_g > 0.04045) {
			var_g = pow(((var_g + 0.055) / 1.055), 2.4);
		} else {
			var_g = var_g / 12.92;
		}

		if (var_b > 0.04045) {
			var_b = pow(((var_b + 0.055) / 1.055), 2.4);
		} else {
			var_b = var_b / 12.92;
		}

		var_r = var_r * 100.0;
		var_g = var_g * 100.0;
		var_b = var_b * 100.0;

		//Observer. = 2°, Illuminant = D65
		return pixel<float> (
			var_r * 0.4124 + var_g * 0.3576 + var_b * 0.1805,
			var_r * 0.2126 + var_g * 0.7152 + var_b * 0.0722,
			var_r * 0.0193 + var_g * 0.1192 + var_b * 0.9505
		);
	}

	pixel<float> rgb_to_abc()
	{
		double var_r = (r / 255.0);
		double var_g = (g / 255.0);
		double var_b = (b / 255.0);

		if (var_r > 0.04045 ) {
			var_r = pow(((var_r + 0.055) / 1.055), 0.4);
		} else {                   
			var_r = var_r / 12.92;
		}

		if (var_g > 0.04045) {
			var_g = pow(((var_g + 0.055) / 1.055), 0.4);
		} else {
			var_g = var_g / 12.92;
		}

		if (var_b > 0.04045 ) {
			var_b = pow(((var_b + 0.055) / 1.055), 0.4);
		} else {
			var_b = var_b / 12.92;
		}

		var_r = var_r * 100.0;
		var_g = var_g * 100.0;
		var_b = var_b * 100.0;

		//Observer. = 2°, Illuminant = D65
		return pixel<float> (
			var_r * 0.4124 + var_g * 0.3576 + var_b * 0.1805,
			var_r * 0.2126 + var_g * 0.7152 + var_b * 0.0722,
			var_r * 0.0193 + var_g * 0.1192 + var_b * 0.9505
		);
	}

	pixel<float> xyz_to_lab()
	{
		double var_x = xyz_x / 95.047; // Observer= 2°, Illuminant= D65
		double var_y = xyz_y / 100.0;          
		double var_z = xyz_z / 108.883;

		if (var_y > 0.008856) {
			var_x = pow(var_x, 0.3333333333);
		} else {
			var_x = (7.787 * var_x) + (16.0 / 116.0);
		}

		if (var_y > 0.008856) {
			var_y = pow(var_y, 0.3333333333);
		} else {
			var_y = (7.787 * var_y) + (16.0 / 116.0);
		}

		if (var_z > 0.008856) {
			var_z = pow(var_z, 0.3333333333);
		} else {
			var_z = (7.787 * var_z) + (16.0 / 116.0);
		}

		return pixel<float> (
			(116.0 * var_y) - 16.0,     // L
			500.0 * (var_x - var_y), // A
			200.0 * (var_y - var_z)  // B
		);
	}

	pixel<float> lab_to_lch()
	{
		float var_h = atan2(lab_b, lab_a);

		if (var_h > 0) {
			var_h = (var_h / PI ) * 180.0;
		} else {
			var_h = 360.0 - (fabs(var_h) / PI) * 180.0;
		}

		return pixel<float> (
			lab_l,
			(float)sqrt(lab_a * lab_a + lab_b * lab_b),
			var_h
		);
	}

	pixel<float> rgb_to_lab()
	{
		return rgb_to_xyz().xyz_to_lab();
	}

	pixel<float> rgb_to_lch()
	{
		return rgb_to_lab().lab_to_lch();
	}

	pixel<int> rgb_to_ycocg()
	{
		int co = rgb_r - rgb_b; // r & b
		int t = rgb_b + ((co + 1) >> 1); // b & co
		int cg = rgb_g - t; // g & t
		int y = t + ((cg + 1) >> 1); // t & cg
		return pixel<int>(y, co, cg);
	}

	pixel<int> ycocg_to_rgb()
	{
		int t = ycocg_y - ((ycocg_cg + 1) >> 1); // y & cg
		int g = ycocg_cg + t; // t & cg
		int b = t - ((ycocg_co + 1) >> 1); // t & co
		int r = b + ycocg_co; // b & co
		return pixel<int>(r, g, b);
	}

	pixel rgb_to_ycbcr(double gamma = 1.0 / 2.2)
	{
		pixel<float> gam = this->gamma(gamma);

		return pixel<float>(
			data_type1(000.0f) + gam.rgb_r * +0.29900f +  gam.rgb_g * +0.587000f + gam.rgb_b * +0.114000f,  // Y
			data_type2(128.0f) + gam.rgb_r * -0.168736f + gam.rgb_g * -0.331264f +  gam.rgb_b * +0.500000f, // Cb
			data_type3(128.0f) + gam.rgb_r * +0.500000f + gam.rgb_g * -0.418688f + gam.rgb_b * -0.081312f   // Cr
		);
	}

	pixel ycbcr_to_rgb(double gamma = 2.2)
	{
		return pixel(
			ycbcr_y +                                                     (ycbcr_cr - 128.0f) * 1.402f ,  // R
			ycbcr_y - (ycbcr_cb - 128.0f) * 0.34414f - (ycbcr_cr - 128.0f) * 0.71414f, // G
			ycbcr_y + (ycbcr_cb - 128.0f)	* 1.77200f   												  // B
		).gamma(gamma);
	}

	pixel<data_type1, data_type2, data_type3> rgb_to_ycbcr_i()
	{
		return pixel<data_type1, data_type2, data_type3>(
			data_type1(000.0f + float(rgb_r) * +0.29900f +  float(rgb_g) * +0.587000f + float(rgb_b) * +0.114000f), // Y
			data_type2(128.0f + float(rgb_r) * -0.168736f + float(rgb_g) * -0.331264f +  float(rgb_b) * +0.500000f), // Cb
			data_type3(128.0f + float(rgb_r) * +0.500000f + float(rgb_g) * -0.418688f + float(rgb_b) * -0.081312f)   // Cr
		);
	}

	pixel<data_type1, data_type2, data_type3, array_type> ycbcr_to_rgb_i()
	{
		return pixel<data_type1, data_type2, data_type3, array_type>(
			data_type1(float(ycbcr_y) +                                                               (float(ycbcr_cr) - 128.0f) * 1.40200f), // R
			data_type2(float(ycbcr_y) - (float(ycbcr_cb) - 128.0f) * 0.34414f - (float(ycbcr_cr) - 128.0f) * 0.71414f), // G
			data_type3(float(ycbcr_y) + (float(ycbcr_cb) - 128.0f) * 1.77200f       								                        )  // B
		);
	}

	pixel<float> lab_to_rgb()
	{
		return lab_to_xyz().xyz_to_rgb();
	}

	// tested
	pixel<int> xyz_to_rgb()
	{
		float X = xyz_x / 100.0f;        //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
		float Y = xyz_y / 100.0f;        //Y from 0 to 100.000
		float Z = xyz_z / 100.0f;        //Z from 0 to 108.883

		float R = X * 3.2406f + Y * -1.5372f + Z * -0.4986f;
		float G = X * -0.9689f + Y *  1.8758f + Z *  0.0415f;
		float B = X *  0.0557f + Y * -0.2040f + Z *  1.0570f;

		if (R > 0.0031308f) R = 1.055f * pow(R, 1.0f / 2.4f) - 0.055f;
		else R = 12.92f * R;

		if (G > 0.0031308f) G = 1.055f * pow(G, 1.0f / 2.4f) - 0.055f;
		else G = 12.92f * G;

		if (B > 0.0031308f) B = 1.055f * pow(B, 1.0f / 2.4f) - 0.055f;
		else B = 12.92f * B;

		return pixel<int>(
			(int)cap_max(cap_min(round(R * 255.0f)), 255),
			(int)cap_max(cap_min(round(G * 255.0f)), 255), 
			(int)cap_max(cap_min(round(B * 255.0f)), 255)
		);
	}

	pixel<float> lab_to_xyz()
	{
		float Y = (lab_l + 16.0f) / 116.0f;
		float X = lab_a / 500.0f + Y;
		float Z = Y - lab_b / 200.0f;

		if (Y * Y * Y > 0.008856 ) Y = Y * Y * Y;
		else Y = (Y - 16.0f / 116.0f) / 7.787f;

		if (X * X * X > 0.008856 ) X = X * X * X;
		else X = (X - 16.0f / 116.0f) / 7.787f;

		if (Z * Z * Z > 0.008856 ) Z = Z * Z * Z;
		else Z = (Z - 16.0f / 116.0f) / 7.787f;

		 // Observer= 2°, Illuminant= D65
		return(pixel<float>(95.047f * X, 100.0f * Y, 108.883f * Z)); 
	}

	// if positive, interpolate to c1 [0, +1]
	// if negative, interopolate to c2 [-1, 0]
	// if zero then return this pixel
	pixel interpolate_sign(double n, pixel c1, pixel c2)
	{
		if (n > 0.0) {
			return alpha_blend(c1, min(1.0, n));
		} else {
			return alpha_blend(c2, max(-1.0, -n));
		}
	}

	pixel debug_byte()
	{
		int c = pixel_case();

		if (c == -1) {
			return pixel(RGB_BLUE.r, RGB_BLUE.g, RGB_BLUE.b);
		} else if (c == 1) {
			return pixel(RGB_ORANGE.r, RGB_ORANGE.g, RGB_ORANGE.b);
		} else {
			return *this;
		}
	}

	int pixel_case()
	{
		if (rgb_r < 0 || rgb_g < 0 || rgb_b < 0) {
			return -1;
		} else if (rgb_r > 255 || rgb_g > 255 || rgb_b > 255) {
			return 1;
		} else {
			return 0;
		}
	}

	static pixel false_color(int seed)
	{
		return pixel(
			(seed * 54321) & 255,
			(seed * 99999) & 255,
			(seed * 33333) & 255
		);
	}

	pixel static median (pixel a, pixel b, pixel c)
	{
		return pixel(
			median_value(a.r, b.r, c.r),
			median_value(a.g, b.g, c.g),
			median_value(a.b, b.b, c.b)
		);
	}

	pixel static median_min (pixel a, pixel b, pixel c, pixel d)
	{
		return pixel(
			median_min_value(a.r, b.r, c.r, d.r),
			median_min_value(a.g, b.g, c.g, d.g),
			median_min_value(a.b, b.b, c.b, d.b)
		);
	}

	pixel static median_max (pixel a, pixel b, pixel c, pixel d)
	{
		return pixel(
			median_max_value(a.r, b.r, c.r, d.r),
			median_max_value(a.g, b.g, c.g, d.g),
			median_max_value(a.b, b.b, c.b, d.b)
		);
	}

	static pixel display(float x, pixel p_base = RGB_WHITE, pixel p_max_out = RGB_RED)
	{
		if (x > 1.0) {
			return p_max_out;
		} else {
			return p_base * x;
		}
	}

	pixel view_channel(
		int p_channel, 
		array_type in_min = (array_type)0, array_type in_max = (array_type)255,
		float p_gamma = 2.2, pixel<array_type> p_color = RGB_WHITE)
	{
		//float chan = 8;//pow((float)scale((float)channel[p_channel], (float)in_min, (float)in_max), p_gamma);
		float chan = pow(rescale(channel[p_channel - 1], in_min, in_max), p_gamma);
		//float chan = channel[p_channel - 1] / 255.0f;
		return p_color * chan;
	}

	inline char* report()
	{
		pixel<float> tmp = *this; // calls conversion constructor
		sprintf_s(color_h_report_text, "(%f, %f, %f)", (float)tmp.rgb_r, (float)tmp.rgb_g, (float)tmp.rgb_b);
		return color_h_report_text;
	}

	inline char* report_int()
	{
		sprintf_s(color_h_report_text, "(%i, %i, %i)", (int)r, (int)g, (int)b);
		return color_h_report_text;
	}
};

typedef pixel<unsigned char> rgb;

#endif // INC_COLOR

