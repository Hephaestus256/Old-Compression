#ifndef INC_STATIC_ARRAY
#define INC_STATIC_ARRAY

template <int size, class key_type = int, class data_type = int>
class static_array {
	int first, last;

public:
	class node_type {
	public:
		key_type key;
		data_type data;

		node_type ()
		{
		}

		node_type (key_type p_key, data_type p_data)
		{
			key = p_key;
			data = p_data;
		}
	};

	node_type* array;

	static_array ()
	{
		first = 0;
		last = -1;
		array = new node_type[size];
	}

	~static_array ()
	{
		delete array;
	}

	int get_start()
	{
		return first;
	}

	int get_end()
	{
		return last;
	}

	int get_length()
	{
		return last - first + 1;
	}

	bool empty()
	{
		return last < first;
	}

	node_type front()
	{
		return array[first];
	}

	node_type back()
	{
		return array[last];
	}

	node_type pop_front ()
	{
		node_type ret = front();
		first++;
		return ret;
	}

	void push_back (const key_type p_key, const data_type p_data)
	{
		array[++last] = node_type(p_key, p_data);
	}

	node_type pop_back ()
	{
		node_type ret = back();
		last--;
		return ret;
	}

	// assumes that the array is sorted
	// executes in log(n)
	inline bool find_key_range (const key_type p_key, key_type& p_min, key_type& p_max, key_type p_guess)
	{
		do { 
			int curr_key = array[p_guess].key;

			if (curr_key > p_key) {
				p_max = p_guess;
			} else if (curr_key < p_key) {
				p_min = p_guess;
			} else {
				return true;
			}

			p_guess = (p_min + p_max) / 2;
		} while (p_min + 1 < p_max);

		return false;
	}

	void merge(node_type* tmp, int start, int mid, int end)
	{
		int i = start, j = mid + 1, n = start;

		while (i <= mid && j <= end) {
			if (array[i].key < array[j].key) {
				tmp[n++] = array[i++];
			} else {
				tmp[n++] = array[j++];
			}
		};

		while (i <= mid) {
			tmp[n++] = array[i++];
		};

		while (j <= end) {
			tmp[n++] = array[j++];
		};

		for (n = start; n <= end; n++) {
			array[n] = tmp[n];
		}
	}

	void sort_array(node_type* tmp, int start = 0, int end = size - 1)
	{
		if (start >= end) return; 

		int mid = (start + end) / 2;

		sort_array(tmp, start, mid);
		sort_array(tmp, mid + 1, end);
		merge(tmp, start, mid, end);
	}

	// performs merge sort on array
	void sort_array(int start = 0, int end = size - 1)
	{
		node_type tmp[size];
		sort_array(tmp, start, end);
	}
};

#endif