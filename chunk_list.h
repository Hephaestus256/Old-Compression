#ifndef INC_CHUNK_LIST
#define INC_CHUNK_LIST


template <class chunk_type = unsigned int, unsigned int chunk_len = 32 / sizeof(chunk_type)>
class chunk_list {

	class chunk {
	public:
		chunk_type data[chunk_len];
		chunk* next;

		chunk()
		{
			next = NULL;
		}
	};

	chunk* first;
	chunk* last;
	unsigned int pos;

public:

	chunk_list ()
	{
		first = last = NULL;
		pos = chunk_len * sizeof(chunk_type) * 8;
	}

	template <class gen>
	void push_back(const gen& p_data, unsigned int p_bits = sizeof(gen) * 8)
	{
		if (pos + p_bits > chunk_len * sizeof(chunk_type) * 8) {
			chunk* c = new chunk;
			
			if (first == NULL) {
				first = c;
			} else {
				last->next = c;
			}

			last = c;

			if (pos >= chunk_len * sizeof(chunk_type) * 8) {
				pos = 0;
			}
		}

		cout << "pos: " << pos << " chunk len: " << chunk_len << endl;
		getchar();

		memcpy((char*)last->data + pos / 8, &p_data, p_bits / 8);
		pos += p_bits;
	}

	void report_cout()
	{
		for (chunk* i = first; i != NULL; i = i->next) {
			if (i->next == NULL) {
				for (unsigned int j = 0; j < pos / (sizeof(chunk_type) * 8); j++) {
					cout << i->data[j] << ", ";
				}
			} else {
				for (unsigned int j = 0; j < chunk_len; j++) {
					cout << i->data[j] << ", ";
				}
			}
		}
		cout << endl;
	}
};

#endif