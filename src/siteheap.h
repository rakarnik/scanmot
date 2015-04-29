#include "scoredsite.h"

class SiteHeap {
	unsigned int maxsize;
	double minscore;
	bool dirty;
	vector<ScoredSite> store;
	CompareScoredSite css;

 public:
	SiteHeap(const unsigned int ms);
	int sitecount() const;
	void add(const ScoredSite& ss);
	void cleanup();
	ScoredSite& get_site(const int num);
};
