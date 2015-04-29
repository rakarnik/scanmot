#ifndef _scoredsites
#define _scoredsites

#include "site.h"

class ScoredSite {
 public:
	double score;
	Site site;
	ScoredSite(const double sc, const Site& st);
	ScoredSite(const double sc, const int c, const int p, const bool s);
};

class CompareScoredSite {
 public:
	bool operator()(const ScoredSite& ss1, const ScoredSite& ss2);
};

#endif
