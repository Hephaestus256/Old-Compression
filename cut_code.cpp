	// Read ABC file and write BMP out

	/*
	if (true) {
		int w, h;
		deque<int> bmp_dat[3];
		bit_file in_file("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\" + fn + ".abc");

		in_file.read(w, 16);
		in_file.read(h, 16);
		bitmap bmp_out(w, h);

		for (int channel = 0; channel < 3; channel++) {
			if (channels[channel]) {
				cout << "********** starting channel: " << channel << endl;

				run_len2<int> r_load(w, h);

				r_load.decode_meta_start();

				short int run;
				do {
					in_file.read(run);
				} while (r_load.decode_meta_next(run));
				
				r_load.decode_meta_end();

				r_load.decode_data_start();

				int dat;
				do {
					in_file.read(dat, 8);
				} while (r_load.decode_data_next(dat));

				r_load.decode_data_end();
		
				cout << "stream count: " << r_load.process_decode() << " expected: " << w * h << endl;
				getchar();

				r_load.export_cont(bmp_dat[channel]);
			}
		}

		bmp_out.load_from_deques(bmp_dat[0], bmp_dat[1], bmp_dat[2], w, h);

		bmp_out.save_to_file("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\" + fn + "_out.bmp");

		bitmap orig("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\" + fn + ".bmp");
		bmp_out.compare(orig);
	}
	*/

	/*
	r.encode_end();

	r.process_meta_start();

	cout << "starting meta" << endl;
	getchar();

	while (r.process_meta_next(run)) {
		cout << "meta: " << run << endl;
		getchar();
	}

	r.process_meta_end();

	r.process_data_start();
	
	while (r.process_data_next(dat)) {
		cout << "data: " << dat << endl;
		getchar();
	};

	r.process_data_end();

	return 0;
	*/

	if (false) {
		run_len<> rle;
		deque<int> orig;
		deque<int> q_data;
		deque<int> last;

		/* load values */

		for (int i = 0; i < 8; i++)
			orig.push_back(10);
	
		for (int i = 0; i < 6; i++)
			orig.push_back(77);

		for (int i = 0; i < 4; i++)
			orig.push_back(25);

		for (int i = 0; i < 5; i++)
			orig.push_back(254);

		/* encode test */
		rle.encode(orig, q_data);

		for (unsigned int i = 0; i < q_data.size(); i++) {
			cout << "enc: " << q_data[i] << ", ";
		}
		getchar();

		/* decode test */
		rle.decode(q_data, last);

		for (unsigned int i = 0; i < last.size(); i++) {
			//cout << "dec: " << last[i] << ", ";
		}

		getchar();
		return 0;

		/*
		run_len<> rle2;

		rle2.decode_start();
		rle2.decode_next(
		*/

		return 0;
	}

	abc_compression abc_out;

	abc_out.test_compression(
		"C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files",
		"lena", "test25"
	);

	return 0;

	abc_out.compress_from_bmp(
		"C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\birds.bmp",
		"C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\birds.abc"
	);

	return 0;

	/*
	bitmap bin("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\rectangle.bmp");

	cout << "type: " << bin.bmp_header.type << endl;
	cout << "size: " << bin.bmp_header.size << endl;
	cout << "res1: " << bin.bmp_header.reserved1 << endl;
	cout << "res2: " << bin.bmp_header.reserved2 << endl;
	cout << "offset: " << bin.bmp_header.offset << endl;

	cout << "size: " << bin.bmp_info.size << endl;
	cout << "width: " << bin.bmp_info.width << endl;
	cout << "height: " << bin.bmp_info.height << endl;
	cout << "planes: " << bin.bmp_info.planes << endl;
	cout << "bits: " << bin.bmp_info.bits << endl;
	cout << "comp: " << bin.bmp_info.compression << endl;
	cout << "image size: " << bin.bmp_info.imagesize << endl;
	cout << "x pixels per meter: " << bin.bmp_info.xresolution << endl;
	cout << "y pixels per meter: " << bin.bmp_info.yresolution << endl;
	cout << "colors: " << bin.bmp_info.ncolors << endl;
	cout << "important colors: " << bin.bmp_info.importantcolors << endl;

	getchar();
	return 0;
	*/

	bitmap bmp_a("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\lena.bmp");
	bitmap bmp_b("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\pepper.bmp");

	for (int y = 0; y < bmp_a.dims.y; y++) {
		for (int x = 0; x < bmp_a.dims.x; x++) {
			//bmp_a.write_rgb(x, y, bmp_b.read_rgb(x, y));
		}
	}

	bmp_a.compare(bmp_b);
	getchar();

	return 0;

	/*
	bit_file file_in("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\huffman.abc");
	huffman<0, 511, 9> huf_in;

	huf_in.read_trie(file_in);
	int len_in = huf_in.decode_start(file_in);
	for (int i = 0; i < len_in; i++) {
		cout << "data: " << (int)huf_in.decode_next(file_in) << endl;
	}

	getchar();
	return 0;
	*/

	/*
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
	*/

	bit_file file2("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\huffman.abc");
	huffman<0, 255, 8> huf2;

	huf2.read_trie(file2);
	huf2.decode_start(file2);

	int data;

	while (huf2.decode_next(file2, data)) {
		cout << "final data: " << data << endl;
		getchar();
	};

	/*
	for (size_t i = 0; i < len; ++i) {
		cout << "value: " << huf2.decode_next(file2) << endl;
	}
	*/

	huf2.decode_end(file2);

	file2.close();

	getchar();
	return 0;

	if (false) {
		bitmap bmp("C:\\Users\\Nathan\\Documents\\Visual Studio 2012\\Test Files\\birds.bmp");

		cout << "width: " << bmp.dims.x << " height: " << bmp.dims.y << endl << endl;
		getchar();

		huffman<-512 * 8, 512 * 8, 13, int> h;

		h.encode_start();
	
		typedef fixd<2> enc_type;
		enc_type predict = 0;
		for (int chan = 0; chan < 3; chan++) {
			for (int y = 0; y < bmp.dims.y; y++) {
				enc_type prev = 0;
				for (int x = 0; x < bmp.dims.x; x++) {
					pixel<enc_type> p = bmp.read_rgb(x, y);
					pixel<enc_type> cpy = p;
					enc_type curr = p.rgb_to_ycbcr_i().channel[chan];
					cpy.channel[chan] = cpy.channel[chan].truncate(2);
				
					if (cpy.ycbcr_to_rgb_i() == p.ycbcr_to_rgb_i()) {
						curr = cpy.rgb_to_ycbcr_i().channel[chan];
					}

					if (y >= 1) {
						pixel<enc_type> pred = bmp.read_rgb(x, y);
						predict = ((pred.rgb_to_ycbcr_i().channel[chan] + fixd<2>(prev)) / fixd<2>(2));
					} else {
						predict = prev;
					}
				
					enc_type delta = curr - predict;
					enc_type delta_cpy = delta;
					delta = delta.round(2);

					h.encode_next(delta.space);
					prev = curr;
				}
			}
		}

		cout << "tally made" << endl;
		getchar();

		h.encode_end();
		//h.send();
		cout << "ratio: " << h.calc_ratio();
		getchar();

		return 0;

		/* end of main */
	}
