#include "scantop.h"

int main(int argc, const char* argv[]) {
	if(argc != 7) {
		cerr << "Usage:\n";
		cerr << "scantop <ace_file> <motif_number> <sequence_file> <background_file> <site_count> <output_file>\n";
		exit(0);
	}

	string ace_file(argv[1]);
	string motnum(argv[2]);
	string seq_file(argv[3]);
	string bg_seq_file(argv[4]);
	unsigned int count = atoi(argv[5]);
	string out_file(argv[6]);

	// Read background sequences and create background model
	cerr << "Reading background sequences ... ";
	vector<string> bgseqs, bgnameset;
	get_fasta_fast(bg_seq_file.c_str(), bgseqs, bgnameset, 0);
	Seqset bg_sseq(bgseqs);
	BGModel bgm(bg_sseq);
	cerr << "done.\n";

	// Read sequences to scan
	cerr << "Reading sequences to scan ... ";
	vector<string> seqs, nameset;
	get_fasta_fast(seq_file.c_str(), seqs, nameset, 0);
	Seqset scan_sseq(seqs);
	cerr << "done.\n";
	// cerr << scan_sseq.num_seqs() << " sequences read\n";

	// Read motif from ACE file
	cerr << "Reading motif... ";
	vector<string> hits;
	string line;
	ifstream acein(ace_file.c_str());
	bool found = false;
	while(acein.good()) {
		getline(acein, line);
		istringstream iss(line);
		vector<string> fields;
		copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter<vector<string> >(fields));
		if(found) {
			if(line[0] != '*') {
				hits.push_back(fields[0]);
			} else {
				break;
			}
		} else {
			if(fields.size() > 0 
				&& (fields[0].compare("Motif") == 0)
				&& (fields[1].compare(motnum) == 0)) {
				found = true;
			}
		}
	}
	if(! found) {
		cerr << "Motif number " << motnum << " not found in ACE file '" << ace_file << "'\n";
		exit(1);
	}
	Seqset ace_sseq(hits);
	Motif mot(ace_sseq, 12, bgm.gcgenome());
	mot.clear_cols();
	for(unsigned int i = 0; i < line.length(); i++) {
		if(line[i] == '*') mot.add_col(i);
	}
	for(unsigned int i = 0; i < hits.size(); i++) {
		mot.add_site(i, 0, 1);
	}
	acein.close();
	cerr << "done.\n";

	// Scan sequences
	cerr << "Scanning sequences...\n";
	vector<list <double> > topsites(scan_sseq.num_seqs());
	double* sm = new double[mot.ncols() * 4];
	double* am = new double[mot.ncols()];
	mot.calc_score_matrix(sm, am);
	int width = mot.get_width();
	double scw, scc, msc, bsc;
	for(int c = 0; c < scan_sseq.num_seqs(); c++) {
		if(scan_sseq.len_seq(c) < mot.get_width()) {
			for(unsigned int i = 0; i < count; i++)
				topsites[c].push_back(0.0);
		}
		for(int p = 0; p < scan_sseq.len_seq(c) - width + 1; p++) {
			msc = mot.score_site(scan_sseq, sm, am, c, p, 1);
			bsc = bgm.score_site(mot, c, p, 1);
			scw = msc - bsc;
			//cerr << c << '\t' << p << '\t' << seqs[c].substr(p, width) << '\t' << 1 << '\t' << setprecision(4) << msc << '\t' << bsc << '\t' << scw << '\n';
			if(topsites[c].size() < count || topsites[c].front() < scw) {
				if(topsites[c].size() == count) {
					topsites[c].pop_front();
				}
				topsites[c].push_front(scw);
				topsites[c].sort();
			}
			msc = mot.score_site(scan_sseq, sm, am, c, p, 0);
			bsc = bgm.score_site(mot, c, p, 0);
			scc = msc - bsc;
			//cerr << c << '\t' << p << '\t' << seqs[c].substr(p, width) << '\t' << 0 << '\t' << setprecision(4) << msc << '\t' << bsc << '\t' << scc << '\n';
			if(topsites[c].size() < count || topsites[c].front() < scc) {
				if(topsites[c].size() == count) {
					topsites[c].pop_front();
				}
				topsites[c].push_front(scc);
				topsites[c].sort();
			}
		}
	}
	delete [] sm;
	delete [] am;

	// Print out the sites
	ofstream scout(out_file.c_str());
	list<double>::iterator it;
	for(int c = 0; c < scan_sseq.num_seqs(); c++) {
		scout << nameset[c];
		// scout << topsites[c].posit() << '\t';
		// scout << topsites[c].posit() + width - 1 << '\t';
		// scout << topsites[c].strand() << "\t";
		for(it = topsites[c].begin(); it != topsites[c].end(); ++it)
			scout << '\t' << *it;
		scout << '\n';
	}
	scout.close();
}
