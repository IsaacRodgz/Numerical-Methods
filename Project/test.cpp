#include<bits/stdc++.h>
#include "Spectral.hpp"

using namespace std;

void test_Spectral(){

    Spectral spectral = Spectral(1.0);

    spectral.readData(1, "data_small.txt");

    //spectral.printData();

    spectral.buildGraph();

    //spectral.printGraph();

    spectral.buildLaplacian();

    spectral.computeEigen(4, 10000, 0.000001);

}

int main(int argc, char const *argv[]) {

    test_Spectral();
}
