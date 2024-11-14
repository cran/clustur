//
// Created by Gregory Johnson on 6/28/24.
//
#include "Adapters/MatrixAdapter.h"
#include <map>
#include <utility>

MatrixAdapter::MatrixAdapter(const std::vector<int> &iIndexes, const std::vector<int> &jIndexes,
                             const std::vector<double> &dataValues, const double cutOff, const bool isSimularity,
                             CountTableAdapter table): cutoff(cutOff), countTable(std::move(table)),xPosition(iIndexes), yPosition(jIndexes),
                            data(dataValues) {
}

SparseDistanceMatrix MatrixAdapter::CreateSparseMatrix() {
    // The indexes are +1, i need to push them back so that 1 -> 0, 2-> 1, etc (name map maybe?)
    std::set<std::string> names;
    SparseDistanceMatrix sparseMatrix;
    const int nSeqs = static_cast<int>(data.size());
    if(nSeqs <= 0)
        return {};
    std::unordered_map<int, int> positionsOfIndexs;
    std::unordered_map<int, std::string> positionsToNames;
    auto samples = countTable.GetSamples();
    names.insert(samples.begin(), samples.end());
    const int nameSize = static_cast<int>(names.size());
    sparseMatrix.resize(nameSize);
    // Not size, but the largest index inside of xPostions
    const int maxXValue = *std::max_element(xPosition.begin(), xPosition.end());
    const int maxYValue = *std::max_element(yPosition.begin(), yPosition.end());
    if(std::max(maxXValue, maxYValue) > nameSize) { // There are values that should exist
        std::set<std::string> unknownNames;
        for(int i = nameSize; i < static_cast<int>(xPosition.size()); i++) {
            if(i >= nameSize + 2)
                break;
            unknownNames.insert(std::to_string(xPosition[i]));
        }
        const Utils util;
        util.CheckForDistanceFileError(unknownNames);
    }
    matrixNames = std::vector<std::string>(nameSize);
    for (int i = 0; i < nameSize; i++) {
        positionsToNames[i] = countTable.GetNameByIndex(i);// Not going to work, I need a way to link my names to the sparse matix indices
        matrixNames[i] = positionsToNames[i];
    }

    for (int i = 0; i < nSeqs;  i++) {

        double currentDist = data[i];
        if(currentDist > cutoff) continue;
        if(currentDist < 0) {
            currentDist = 0;
        }
        const int xIndex = xPosition[i]; // Coming from r -> c++, indeces start at 1 in r
        const int yIndex = yPosition[i];
      
        if(xIndex > yIndex)
            sparseMatrix.addCell(yIndex, PDistCell(xIndex, static_cast<float>(currentDist)));
        else
            sparseMatrix.addCell(xIndex, PDistCell(yIndex, static_cast<float>(currentDist)));
    }
    return sparseMatrix;
}

ListVector MatrixAdapter::CreateListVector() const {
    ListVector vector;
    vector.resize(static_cast<int>(matrixNames.size()));
    int count = 0;
    for(const auto& name : matrixNames) {
        vector.set(count++, name);
    }
    return vector;
}
