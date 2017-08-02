#include "bgmodel.h"

BGModel::BGModel(const Seqset& s) :
seqset(s),
total_seq_len(0),
gc_genome(0),
model(4),
wbgscores(seqset.num_seqs()),
cbgscores(seqset.num_seqs()) {
	const vector<vector <int> > ss_seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	for(int i = 0; i < ss_num_seqs; i++) {
		for(int j = 0; j < seqset.len_seq(i); j++) {
			if(ss_seq[i][j] == 1 || ss_seq[i][j] == 2)
				gc_genome++;
			if(ss_seq[i][j] != -1)
				total_seq_len++;
		}
	}
	gc_genome /= total_seq_len;

	train_background();
	calc_bg_scores();
}

double BGModel::score_site(const Motif& motif, const int c, const int p, const bool s) const {
	double L = 0.0;
	int width = motif.get_width();
	int matpos;
	vector<int>::const_iterator col_iter = motif.first_column();
	vector<int>::const_iterator last_col = motif.last_column();
	if(s) {
		matpos = 0;
		for(; col_iter != last_col; ++col_iter) {
			assert(p + *col_iter >= 0);
			assert(p + *col_iter < seqset.len_seq(c));
			L += wbgscores[c][p + *col_iter];
			matpos += 4;
		}
	} else {
		matpos = 0;
		for(; col_iter != last_col; ++col_iter) {
			assert(p + width - 1 - *col_iter >= 0);
			assert(p + width - 1 - *col_iter < seqset.len_seq(c));
			L += cbgscores[c][p + width - 1 - *col_iter];
			matpos += 4;
		}
	}
	return L;
}

void BGModel::train_background() {
	for(int i = 0; i < 4; i++) {
		model[i] = 0;
	}
	
	// Just use genome-wide GC content
	model[0] = (1 - gc_genome)/2;
	model[1] = gc_genome/2;
	model[2] = model[1];
	model[3] = model[0];
}

void BGModel::calc_bg_scores() {
	const vector<vector <int> >& ss_seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = ss_seq[i].size();
		wbgscores[i].clear();
		cbgscores[i].clear();
		wbgscores[i].reserve(len);
		cbgscores[i].reserve(len);

		float bgavg = 0.0;
		bgavg += 2 * (1.0 - gc_genome) * log2(0.5 - gc_genome/2.0);
		bgavg += 2 * gc_genome * log2(gc_genome/2.0);
		
		// Use 0th order model for all bases
		for(int j = 0; j < len; j++) {
			if(ss_seq[i][j] != -1) {
				wbgscores[i].push_back(log2(model[ss_seq[i][j]]));
			} else {
				wbgscores[i].push_back(bgavg);
			}
		}
		for(int j = 0; j < len; j++) {
			if(ss_seq[i][j] != -1) {
				cbgscores[i].push_back(log2(model[(3 - ss_seq[i][j])]));
			} else {
				cbgscores[i].push_back(bgavg);
			}
		}
		
		assert(wbgscores[i].size() == (unsigned int) len);
		assert(cbgscores[i].size() == (unsigned int) len);
	}
}

void BGModel::print_model(ofstream& out) {
	map<int,char> nt;
	nt[0] = 'A';
	nt[1] = 'C';
	nt[2] = 'G';
	nt[3] = 'T';

	out << "Model order:\n";
	for(int i = 0; i < 4; i++) {
		out << nt[i % 4] << '\t' << model[i] << '\n';
	}
	out << "\n\n";
}
