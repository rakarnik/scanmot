#ifndef _bgmodel
#define _bgmodel

#include "standard.h"
#include "seqset.h"
#include "motif.h"

class BGModel {
	const Seqset& seqset;
	int total_seq_len;
	int order;
	float gc_genome;
	vector<float> model;
	vector<vector <float> > wbgscores;
	vector<vector <float> > cbgscores;

	void train_background();                                     // Train 0th order background model
	void calc_bg_scores();                                       // Calculate scores using 0th order background model

public:
	BGModel(const Seqset& s);
	BGModel(const Seqset& s, const float gc);
	float tot_seq_len() const { return total_seq_len; }           // Return total length of all sequences
	float gcgenome() const { return gc_genome; }                  // Return overall GC content
	double score_site(const Motif& motif, const int c, const int p, const bool s) const;
	void print_model(ofstream& out);
};

#endif
