#include "scanmot.h"

int main(int argc, const char* argv[]) {
	if(argc != 6) {
		cerr << "Usage:\n";
		cerr << "scanmot <ace_file> <motif_number> <sequence_file> <output_file> <max_sites>\n";
		exit(0);
	}

	string ace_file(argv[1]);
	string motnum(argv[2]);
	string seq_file(argv[3]);
	string out_file(argv[4]);
	unsigned int ms = atoi(argv[5]);

	// Read sequences and create background model
	vector<string> seqs, nameset;
	get_fasta_fast(seq_file.c_str(), seqs, nameset, 12);
	Seqset scan_sseq(seqs);
	BGModel bgm(scan_sseq);

	// Read motif from ACE file
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

	// Scan sequences
	SiteHeap sh(ms);
	double* sm = new double[mot.ncols() * 4];
	double* am = new double[mot.ncols()];
	mot.calc_score_matrix(sm, am);
	int width = mot.get_width();
	int n = 0;
	double msc, bsc, scw, scc;
	for(int c = 0; c < scan_sseq.num_seqs(); c++) {
		for(int p = 0; p < scan_sseq.len_seq(c) - width; p++) {
			msc = mot.score_site(scan_sseq, sm, am, c, p, 1);
			bsc = bgm.score_site(mot, c, p, 1);
			scw = msc - bsc;
			// cerr << c << '\t' << p << '\t' << seqs[c].substr(p, width) << '\t' << 1 << '\t' << setprecision(4) << msc << '\t' << bsc << '\t' << scw << '\n'; 
			ScoredSite ssw(scw, c, p, 1);
			sh.add(ssw);
			msc = mot.score_site(scan_sseq, sm, am, c, p, 0);
			bsc = bgm.score_site(mot, c, p, 0);
			scc = msc - bsc;
			// cerr << c << '\t' << p << '\t' << seqs[c].substr(p, width) << '\t' << 0 << '\t' << setprecision(4) << msc << '\t' << bsc << '\t' << scc << '\n';
			ScoredSite ssc(scc, c, p, 0);
			sh.add(ssc);
			n += 2;
		}
	}
	delete [] sm;

	// Print out the sites
	ofstream scout(out_file.c_str());
	for(int i = 0; i < sh.sitecount(); i++) {
		ScoredSite bestss = sh.get_site(i);
		scout << nameset[bestss.site.chrom()] << '\t';
		scout << bestss.site.posit() << '\t';
		scout << bestss.site.posit() + width - 1 << '\t';
		scout << bestss.site.strand() << "\t";
		if(bestss.site.strand())
			scout << seqs[bestss.site.chrom()].substr(bestss.site.posit(), width) << '\t';
		else
			scout << reverse_comp(seqs[bestss.site.chrom()].substr(bestss.site.posit(), width)) << '\t';
		scout << bestss.score << '\n';
	}
	scout.close();
}
