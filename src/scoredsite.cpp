#include "scoredsite.h"

ScoredSite::ScoredSite(const double sc, const Site& st):
score(sc),
site(st) {
}

ScoredSite::ScoredSite(const double sc, const int c, const int p, const bool s):
score(sc),
site(c, p, s) {
}

bool CompareScoredSite::operator()(const ScoredSite& ss1, const ScoredSite& ss2) {
	return ss1.score > ss2.score;
}
