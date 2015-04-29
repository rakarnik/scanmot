#include "siteheap.h"

SiteHeap::SiteHeap(const unsigned int ms):
maxsize(ms),
minscore(DBL_MAX),
dirty(false) {
}

int SiteHeap::sitecount() const {
	return store.size();
}

void SiteHeap::add(const ScoredSite& ss) {
	if(ss.score > minscore || store.size() < maxsize) {
		store.insert(store.end(), ss);
		if(ss.score < minscore)
			minscore = ss.score;
		dirty = true;
	}

	if(store.size() > 5 * maxsize)
		cleanup();
}

void SiteHeap::cleanup() {
	sort(store.begin(), store.end(), css);
	vector<ScoredSite> tmp;
	int newsize = maxsize < store.size()?  maxsize : store.size();
	tmp.assign(store.begin(), store.begin() + newsize);
	store.clear();
	store.assign(tmp.begin(), tmp.end());
	minscore = store.back().score;
	dirty = false;
}

ScoredSite& SiteHeap::get_site(const int num) {
	if(dirty)
		cleanup();
	return store[num];
}
