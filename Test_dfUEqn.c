#include <stdio.h>

#include "dfUEqn.H"
#include "dfMatrixDataBase.H"

void makeTestData(){
    num_cells = 0;

    double arr[3] = {1.1, 1.2, 1.3};
    volume = arr;
    //printf("%f\n",volume);
    //printf("%f\n",volume[1]);

    rdelta_t = 0;
}

double test0(double volume,
            double i){
    return((double)(volume + i));
}

int main(int argc, char const *argv[]){
    //makeTestData();
    //process();
    //double arr[3] = {1.1, 1.2, 1.3};
    //volume = arr;
    //printf("%f\n",volume[0]);
    //test0(volume);
    //test0(&volume);
    //printf("%f\n",volume[0]);
    for (int i = 0; i< 1; i++){
        int num = 3;
    }
    for (int i = 0; i< 1; i++){
        int num = 4;
    }
    double re = test0(1. , 2.);
    printf("%f\n",re);
}
