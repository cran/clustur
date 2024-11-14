//
// Created by Gregory Johnson on 8/1/24.
//

#ifndef OTUABUNDANCEPAIR_H
#define OTUABUNDANCEPAIR_H
#include <string>
#include <utility>
struct OTUAbundancePair {
    OTUAbundancePair(const int abundance, std::string otu, std::string group, std::string label)
        : abundance(abundance),
          otu(std::move(otu)),
          group(std::move(group)),
          label(std::move(label)){
    }
    OTUAbundancePair() = default;

    int abundance = 0;
    std::string otu;
    std::string group;
    std::string label;
};

struct SharedAbundance {
    SharedAbundance(std::string group, std::string otu,
        std::string label, const double group_abundance)
        : group(std::move(group)),
          OTU(std::move(otu)),
          label(std::move(label)),
          groupAbundance(group_abundance) {
    }

    std::string group;
    std::string OTU;
    std::string label;
    double groupAbundance;
};
#endif //OTUABUNDANCEPAIR_H
