//
// Created by gregj on 7/9/2024.
//

#ifndef CLUSTERTESTFIXTURE_H
#define CLUSTERTESTFIXTURE_H
#include "TestFixture.h"
#include "../MothurDependencies/Cluster.h"


class ClusterTestFixture : public TestFixture , public Cluster {
public:
    ClusterTestFixture() = default;
    bool TestGetTagGetter(Cluster *clust, const std::string &expectedResult);
    bool TestGetSeqToBin(Cluster* clust, const std::map<std::string, int>& expectedMap);
    bool TestClusterBins(RAbundVector *rAbdund, SparseDistanceMatrix *distanceMatrix, bool expectedResult);
    bool TestClusterNames(ListVector *listVector, SparseDistanceMatrix *distanceMatrix, bool expectedResult);
    bool TestUpdateMap(ListVector *listVector, SparseDistanceMatrix *distanceMatrix, const bool expectedResult);
    ~ClusterTestFixture() override = default;

private:
    void Setup() override;
    void TearDown() override;
    std::string getTag() override {return "";}
    bool updateDistance(PDistCell &colCell, PDistCell &rowCell) override {return false;}

private:
    Cluster* cluster = nullptr;
};



#endif //CLUSTERTESTFIXTURE_H
