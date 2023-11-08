#include <stdio.h>

#include "dfUEqn.H"
#include "dfMatrixDataBase.H"

void makeTestData(){
    num_cells = 0;

    double arr[3] = {1.1, 1.2, 1.3};
    volume = arr;
    printf("%f\n",volume);
    printf("%f\n",volume[1]);

    rdelta_t = 0;
}

int main(int argc, char const *argv[]){
    //makeTestData();
    //process();
    double arr[3] = {1.1, 1.2, 1.3};
    volume = arr;
    printf("%f\n",volume[0]);
    test0(volume);
    //test0(&volume);
    printf("%f\n",volume[0]);
}

void test0(double *volume){
    volume[0] = 99999999.;
}