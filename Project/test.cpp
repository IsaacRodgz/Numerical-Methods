#include<bits/stdc++.h>
#include "Spectral.hpp"

using namespace std;

void test_Spectral(){

    Spectral spectral = Spectral(1.0);
    int eigenSize = 3;
    int numClusters = 4;

    cout << "\n Reading data...\n" << endl;
    spectral.readData(1, "data.txt");

    //spectral.printData();

    cout << "\n Building graph...\n" << endl;
    spectral.buildGraph();

    //spectral.printGraph();

    cout << "\n Building laplacian...\n" << endl;
    spectral.buildLaplacian();

    //spectral.printLaplacian();

    cout << "\n Solving eigensystem...\n" << endl;
    spectral.computeEigen(eigenSize, 10000, 0.000001);

    //spectral.printEigen();

    //spectral.printLaplacian();

    cout << "\n Transforming data...\n" << endl;
    spectral.transformData(eigenSize);
    //spectral.printEigen();

    cout << "\n Clustering points...\n" << endl;
    spectral.kMeansEig(500, numClusters);

    cout << "\nClustering result: \n" << endl;
    spectral.printClusters();
    
    //spectral.kMeans(500, numClusters);
    //spectral.plot();

}

int main(int argc, char const *argv[]) {

    test_Spectral();
}
