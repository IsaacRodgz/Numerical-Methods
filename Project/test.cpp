#include<bits/stdc++.h>
#include "Spectral.hpp"

using namespace std;

void test_Spectral(){

    Spectral spectral = Spectral(1.0);
    int eigenSize = 3;
    int numClusters = 4;

    spectral.readData(1, "data.txt");

    //spectral.printData();

    spectral.buildGraph();

    //spectral.printGraph();

    spectral.buildLaplacian();

    //spectral.printLaplacian();

    spectral.computeEigen(eigenSize, 10000, 0.000001);

    //spectral.printEigen();

    //spectral.printLaplacian();

    spectral.transformData(eigenSize);
    //spectral.printEigen();

    spectral.kMeansEig(500, numClusters);

    //spectral.kMeans(500, numClusters);

    cout << "\nClustering result: \n" << endl;
    spectral.printClusters();
}

int main(int argc, char const *argv[]) {

    test_Spectral();
}
