#ifndef INC_FIXD
#define INC_FIXD

template <int frac_bits = 16, class space_type = int>
class fixd {
public:
	space_type space;

	fixd() 
	{
	}

	template <class gen>
	fixd(gen n)
	{
		space = space_type(n * gen((space_type)1 << frac_bits) + gen(0.5));
	}

	int round(int p_decimal = 0)
	{
		fixd ret = *this;
		ret.space += (1 << (frac_bits - 1));
		ret.space &= (space_type(-1) << (frac_bits - p_decimal));
		return ret;
	}

	// truncate to bit place (decimal = 0 = 1s place, + is frac, - is whole >= 2)
	fixd truncate(int p_decimal = 0)
	{
		fixd ret = *this;
		ret.space &= space_type(-1) << (frac_bits - p_decimal);
		return ret;
	}

	operator unsigned char()
	{
		return (unsigned char)space >> (unsigned char)frac_bits;
	}

	operator char()
	{
		return (char)space >> (char)frac_bits;
	}

	operator unsigned int()
	{
		return (unsigned int)space >> (unsigned int)frac_bits;
	}

	operator int()
	{
		return (int)space >> (int)frac_bits;
	}

	operator float()
	{
		return float(space) * (1.0f / float((space_type)1 << (space_type)frac_bits));
	}

	operator double()
	{
		return double(space) * (1.0 / double((space_type)1 << (space_type)frac_bits));
	}

	fixd operator + (fixd n)
	{
		fixd ret;
		ret.space = space + n.space;
		return ret;
	}

	fixd operator - (fixd n)
	{
		fixd ret;
		ret.space = space - n.space;
		return ret;
	}

	fixd operator * (fixd n)
	{
		fixd ret;
		ret.space = ((__int64)space * (__int64)n.space) >> frac_bits;
		return ret;
	}

	fixd operator / (fixd n)
	{
		fixd ret;
		ret.space = ((__int64)space << frac_bits) / n.space;
		return ret;
	}

	bool operator == (fixd n)
	{
		return space == n.space;
	}

	bool operator != (fixd n)
	{
		return space != n.space;
	}
};

template <int bits, class gen>
int round (fixd<bits, gen> n)
{
	return n.round();
}

template <int bits, class gen>
int floor (fixd<bits, gen> n)
{
	return n.floor();
}

#endif