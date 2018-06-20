#ifndef INC_VARINT
#define INC_VARINT

template <int seg_len = 8, bool signd = false, class cont_type = unsigned int, int base = 0, int xtra_n = 1>
class varint {

	inline cont_type segment_mask (int n)
	{
		return ((cont_type)1 << cont_type(seg_len * n)) - (cont_type)1;
	}

	inline void set_extend_bit (int n, bool b = true)
	{
		const int pos = (cont_type(seg_len * n) - 1);

		if (b) {
			container |= (cont_type)1 << pos;
		} else {
			container &= !((cont_type)1 << pos);
		}
	}

public:

	cont_type container;

	varint ()
	{
	}

	inline varint (int n)
	{
		n -= base;
		int seg_ct = 1;

		if (signd) {
			bool s;
			
			if (xtra_n == 0) { // no extra number for pos or neg
				s = (n < 0);
				if (n < 0) n = -n;
			} else if (xtra_n == 1) { // an extra number for pos
				s = (n <= 0);
				if (n > 0) { // -1: 0011, 0: 0001, 1: 0000, 2: 0010, 3: 0100
					n--;
				} else {
					n = -n;
				}
			} else { // an extra number for neg
				s = (n < 0);
				if (n < 0) n = -1 - n;
			}

			n <<= 1;
			n += (cont_type)s;
		}

		cont_type range = (cont_type)(1 << ((seg_len - 1) * seg_ct));
		container = (cont_type)0;

		while ((int)range <= n) {
			set_extend_bit(seg_ct);
			seg_ct++;
			range += cont_type(1 << ((seg_len - 1) * seg_ct));
		};
		range -= cont_type(1 << ((seg_len - 1) * seg_ct));

		for (int i = 0; i < seg_ct; ++i) {
			cont_type mask = (((1 << (seg_len - 1)) - 1) << i * (seg_len - 1));
			container |= (((n - range) & mask) << i);
		}
	}

	operator int()
	{
		cont_type cont = 0;
		cont_type mask = (1 << seg_len) - 1;
		int n = 0;
		int ret;

		do {
			cont += (container & mask) >> n;
			mask = mask << seg_len;
			n++;
		} while (container & ((1 << (n * seg_len)) - 1));

		ret = base + int(cont);

		/* todo: make casting to int reflect the extra N */
		if (signd) {
			int s = (ret & 1);
			ret >>= 1;
			
			if (xtra_n == 0) {
				if (s) ret = -ret;
			} else if (xtra_n == 1) { // -1: 0011, 0: 0001, 1: 0000, 2: 0010, 3: 0100
				if (s) {
					ret = -ret;
				} else {
					ret++;
				}
			} else { // xtra is -1
				if (s) {
					ret = -ret;
					ret--;
				}
			}
		} 

		return ret;
	}

	inline void clear()
	{
		container = 0;
	}

	inline size_t size()
	{
		int ct = 1;
		int mask = 1 << (seg_len - 1);

		while (container & mask) {
			mask = mask << seg_len;
		}

		return ct * seg_len;
	}

	// reads the nth seg value, including control bit
	inline cont_type read_seg(int seg)
	{
		return (container >> (seg * seg_len)) & ((1 << seg_len) - 1);
	}

	inline cont_type read_seg_value(int seg)
	{
		return (container >> (seg * seg_len)) & ((1 << (seg_len - 1)) - 1);
	}

	inline bool read_seg_ctrl(int seg)
	{
		const int pos = (seg + 1) * seg_len - 1;
		return (container & (1 << pos)) > 0;
	}
	
	inline bool read_continue(int seg)
	{
		if (seg == 0) {
			return true;
		} else {
			return read_seg_ctrl(seg - 1);
		}
	}

	// writes nth segment including control bit
	// returs control bit
	inline bool write_seg(cont_type p_value, int p_seg)
	{
		if (p_seg == 0) container = 0;

		container |= (p_value << (p_seg * seg_len));
		return (p_value & (1 << ((p_seg + 1) * seg_len - 1))) != 0;
	}
};

typedef varint<8, false> varint_unsigned;
typedef varint<8, true> varint_signed;

#endif