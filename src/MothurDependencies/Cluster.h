//
// Created by Gregory Johnson on 6/14/24.
//

#ifndef CLUSTER_H
#define CLUSTER_H



#include "SparseDistanceMatrix.h"
#include "OptiMatrix.h"
#include "PDistCell.h"
#include "RAbundVector.h"

class ListVector;

// typedef unsigned long ull;
class Cluster  {

public:
	Cluster(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, std::string, float);
    Cluster() = default;
    virtual ~Cluster() = default;
    virtual bool update(double&);
	virtual std::string getTag() = 0;
	// virtual void setMapWanted(bool m);
	virtual std::map<std::string, int> getSeqtoBin()  {  return seq2Bin;	}
	virtual bool updateDistance(PDistCell& colCell, PDistCell& rowCell) = 0;
protected:
	virtual bool clusterBins();
	virtual bool clusterNames();
	virtual bool updateMap();


	RAbundVector* rabund{};
	ListVector* list{};
	SparseDistanceMatrix* dMatrix{};

	unsigned long smallRow{}, smallCol{}, nRowCells{}, nColCells{};
	float cutoff{};
	std::string method;
	float smallDist{}, adjust{};
	bool mapWanted{};
	std::map<std::string, int> seq2Bin;


    Utils util;
};


#endif //CLUSTER_H
