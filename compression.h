#ifndef COMPRESSION
#define COMPRESSION

#include <vector>
#include <deque>
#include <unordered_map>
#include <string>

using namespace std;

class compression {
	
	class rep {
	public:
		int start;
		string word;

		rep (const int p_start, const string p_word)
		{
			start = p_start;
			word = p_word;
		}
	};

	class span_set {
		vector<int> span;
		int tail, head;
		int prev_sign;
		int end_of_span;
		int size;
		int word_len;

	public:

		span_set ()
		{
		}

		int get_iterator ()
		{
			return head;
		}

		void init (int p_word_len, int p_size) 
		{
			size = p_size;
			word_len = p_word_len;
			span.reserve(p_size);
			span[0] = p_size;
		}

		void skip_nulls ()
		{
			if (head + word_len - 1 >= end_of_span) {
				while (span[head] < 0) { // skip past negative spans
					head += - span[head];
				}

				end_of_span = head + span[head];
			}
		}

		void start ()
		{
			tail = head = 0;
			prev_sign = -1;
			end_of_span = 0;

			skip_nulls();
		}

		bool next (bool p_cond)
		{
			int sign = p_cond ? 1: -1;

			if (sign != prev_sign) {
				span[tail] = (head - tail) * prev_sign;
				tail = head;
			}

			prev_sign = sign;

			return navigate();
		}

		bool navigate ()
		{
			skip_nulls();
			head++;
			return head + word_len - 1 < size;
		}

		void end ()
		{
			span[tail] = (head - tail) * prev_sign;
		}
	};

public:

	string data;
	unordered_map<string, int> dict;
	deque<rep> reps;
	span_set spans;
	string curr_string;
	unordered_map<string, int>::iterator curr_find;

	compression(const string p_data)
	{
		data = p_data;
		spans.init(p_data.length());
	}

	bool find (int i, int p_len)
	{
		curr_string = data.substr(i, p_len);
		curr_find = dict.find(curr_string);
		return curr_find != dict.end();
	}

	void init_spans()
	{
		spans.start();
	}

	bool count_freq (const int p_len)
	{
		if (find(spans.get_iterator(), p_len)) { // string found
			cout << "Inc: " << curr_string << endl;
			curr_find->second++;
		} else { // string NOT found
			cout << "Inserting: " << curr_string << endl;
			dict.insert(pair<string, int>(curr_string, 1));
		}
		return spans.navigate();
	}

	bool proc_freq (const int p_len)
	{
		bool b = find(spans.get_iterator(), p_len) && curr_find->second > 1;
		return spans.next(b);
	}

	void process(int p_len)
	{
		int n = 0;

		// determine frequencies
		while (count_freq(p_len));

		// process repeats
		n = 0;

		for (int i = 0; i < (int)data.length() - p_len + 1; ++i) {
			int l = skip[n];

			//cout << "l: " << l << endl;

			if (l < 0) {
				i += -l;
				n += -l;
			} else {
				proc_freq(i, l, p_len);
				i += l - p_len;
				n += l;
			} // ! if (l < 0)
		} // for (int i = 0
		
		cout << "span start: " << span_start << endl;
		cout << "skip start: " << skip_start << endl;
		getchar();

		if (cont) {
			//skip[span_start] = (int)data.length() - span_start;
		} else {
			//skip[skip_start] = -((int)data.length() - skip_start);
		}

		cout << endl;
		for (int i = 0; i < data.size(); ++i) {
			//cout << "skip " << i << ": " << skip[i] << endl;
		}
	}

	void encode()
	{
		cout << "processing: " << data << endl;

		process(1);
		process(2);
	}
};

#ifdef abc123
// compression.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <deque>
#include <unordered_map>
#include <string>
#include <iostream>
#include "C:\Users\Nathan\Documents\Visual Studio 2012\Libraries\compression.h"

using namespace std;

class part {
public:
	int start;
	string word;

	part (int p_start, string p_word)
	{
		start = p_start;
		word = p_word;
	}
};

class word_stat {
public:
	int count;
	int last;

	word_stat (int p_count, int p_last)
	{
		count = p_count;
		last = p_last;
	}
};

void remove_oneoffs (
	const vector<part>& p_curr_data, unordered_map<string, word_stat>& p_freqs, deque<part>& p_dict, const int p_word_len)
{
	int pos;

	if (p_curr_data.empty()) {
		return;
	}

	/*
	cout << "Parsing (" << p_word_len << "): ";

	for (int i = 0; i < p_curr_data.size(); i++) {
		cout << p_curr_data[i] << ", ";
	}

	cout << endl;
	*/

	pos = 0;
	for (unsigned int i = 0; i < p_curr_data.size(); i++) {
		//cout << "start a: " << p_curr_data[i].start << endl;
		for (unsigned int j = 0; j <= p_curr_data[i].word.length() - p_word_len; j++) {
			string s = p_curr_data[i].word.substr(j, p_word_len);
			unordered_map<string, word_stat>::iterator f;

			f = p_freqs.find(s);

			if (f == p_freqs.end()) {
				p_freqs.insert(pair<string, word_stat>(s, word_stat(1, pos)));
			} else if (f->second.last + p_word_len - 1 < pos) {
				f->second.count++;
				f->second.last = pos;
				//cout << "pos a: str: " << s << " pos: " << pos << endl;
			}
			pos++;
		} // for (unsigned int j = 0
		pos += p_word_len - 1;
	} // for (unsigned int i = 0

	for (unordered_map<string, word_stat>::iterator i = p_freqs.begin(); i != p_freqs.end(); ++i) {
		i->second.last = -p_word_len - 1;
	}

	string curr;
	int start;
	unsigned int j;

	vector<part> next_data;
	pos = 0;

	for (unsigned int i = 0; i < p_curr_data.size(); i++) {
		start = 0;

		for (j = 0; j <= p_curr_data[i].word.length() - p_word_len; j++) {
			if (true) {
				string s = p_curr_data[i].word.substr(j, p_word_len);
				unordered_map<string, word_stat>::iterator f;

				f = p_freqs.find(s);

				if (f->second.count > 1 && f->second.last + p_word_len - 1 < pos) {
					//cout << "found: '" << s << "' len: " << p_word_len << " ct: " << f->second.count << " start: " << start <<endl;
					//cout << "found: " << s << " pos: " << pos << " len: " << p_word_len << endl;
					//p_lens[pos] = p_word_len;
					p_dict.push_back(part(p_curr_data[i].start + pos, s));
					f->second.last = pos;
					//cout << "str: " << s << " pos: " << pos << endl;
				} else {
					if (j - start + p_word_len - 1 > 0) {
						curr = p_curr_data[i].word.substr(start, j - start + p_word_len - 1);
			
						if (p_word_len + 1 <= curr.length()) {
							next_data.push_back(part(p_curr_data[i].start, curr));
							start = j + p_word_len - 1;
						} else {
							start = j + 1;
						}
					} else {
						start = j + 1;
					}
				}
			} // if (j >= start 

			pos++;
		} // for (j = 0

		if (j - start + p_word_len - 1 > 0) {
			curr = p_curr_data[i].word.substr(start, j - start + p_word_len - 1);

			if (p_word_len + 1 <= curr.length()) {
				next_data.push_back(part(p_curr_data[i].start, curr));
			}
		}

		pos += p_word_len - 1;
	} // for (unsigned int i = 0

	remove_oneoffs(next_data, p_freqs, p_dict, p_word_len + 1);
}

void compress (const string p_str_in, const int p_word_len = 1)
{
	vector<part> curr_data;
	unordered_map<string, word_stat> freqs;
	deque<part> dict;

	curr_data.push_back(part(0, p_str_in));
	remove_oneoffs (curr_data, freqs, dict, p_word_len);

	cout << "defining: " << p_str_in << endl;

	for (int i = 0; i < dict.size(); ++i) {
		//cout << "part: start: " << dict[i].start << " word: " << dict[i].word << endl;
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	compression comp("abc123abc456abc789");

	comp.encode();

	getchar();
	return 0;

	//compress("abcdefgflyingmonkeysgfedcbaflyingmonkeysoqrsmtabuv");
	//compress("xyzabcdefghhhabcdefguvwh");
	//compress("abcxxxxabcyyyy");
	//compress("abbabba123123123"); // X
	//compress("zzzzzzzzzz"); // ? need run len detect
	//compress("abcd123ab456cd");
	//compress("lsjglirrbklwellijb;kdf;;lkwwl,dffkjaww;ldsdlbsdfl;jsdddhdfblalhsfblknsdd;lfbjdkmawd;lmcgnljnwkmbknwjnbl,aljnfl,wkb");

	/*
	for (unsigned int i = 0; i < out.size(); ++i) {
		for (unsigned int j = 0; j < out[i]->size(); ++j) {
			cout << "out " << i << ": " << (*out[i])[j] << ", ";
		}
	}
	*/

	getchar();
	return 0;
}

#endif

#endif