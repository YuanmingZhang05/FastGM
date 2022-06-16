#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <random>
#include <map>
#include <set>
using namespace std;

vector<default_random_engine> generator_vector;
double v_sum = 0;

struct NextCUST{
    double x_i;
    //int z_i;
    int c;
};

struct GBM_org_var{
    double y_j;
    int s_j;
};

vector<NextCUST> CUSTVec;
vector<unsigned int> PermutationVec;
vector<GBM_org_var> GBM_sketch;

void Initial(const vector<double> data, const uint32_t k){
    int n = data.size();
    v_sum = 0;
    NextCUST CUST;
    GBM_org_var GBM_VAR;
    CUSTVec.resize(n);
    PermutationVec.resize(n*k);
    GBM_sketch.resize(k);

    for(uint32_t i=0;i<n;i++) {
        v_sum = v_sum + data[i];
        CUST.x_i = 0;
        //CUST.z_i = 0;
        CUSTVec[i] = CUST;
        for(int j=0;j<k;j++) PermutationVec[i*k+j] = j;
    }

    for(uint32_t i=0;i<k;i++) {
        GBM_VAR.y_j = -1;
        GBM_VAR.s_j = -1;
        GBM_sketch[i] = GBM_VAR;
    }
}

void AscEXP_STM(int i, int k, int z){
    default_random_engine generator = generator_vector[i];
    uniform_int_distribution<int> distribution1(1, RAND_MAX-1);
    int r = distribution1(generator);
    double u = (r) / double(RAND_MAX);
    double x=0;
    int j=0;

    //z = CUSTVec[i].z_i+1;
    j = (r % (k-z+1))+z;
    x = -(log(u))/(k+1-z);


    int temp = PermutationVec[i * k + z-1];
    PermutationVec[i * k + z-1] = PermutationVec[i * k + j];
    PermutationVec[i * k + j] = temp;

    CUSTVec[i].x_i += x;
    //CUSTVec[i].z_i = z;
    CUSTVec[i].c = PermutationVec[i*k + z-1];
    generator_vector[i] = generator;
}

int argmax(vector<GBM_org_var> GBM_sketch){
    int j_star=0;
    double max = GBM_sketch[0].y_j;
    for (int i = 0; i < GBM_sketch.size(); i++) {
        if (max < GBM_sketch[i].y_j){
            max = GBM_sketch[i].y_j;
            j_star = i;
        }
    }   
    return j_star;
}

vector<GBM_org_var> FastGM(const vector<double> data, const uint32_t k, int seed){
    uint32_t delta = 0, k_star = k;   
    int j_star = 0, Prune_Module_Flag = false;
    for (int i=0;i<data.size();i++){
        double v_i = data[i];
        default_random_engine generator((i << 16) + seed);
        generator_vector.push_back(generator);
        if(v_i == 0) continue;
        for (int l = 1; l < k+1; l++){
            AscEXP_STM(i, k, l);
            double b_i = CUSTVec[i].x_i /(v_i);
            int c = CUSTVec[i].c;
            if (Prune_Module_Flag == false){
                if (GBM_sketch[c].y_j < 0 ){
                        GBM_sketch[c].y_j = b_i;
                        k_star = k_star - 1;
                        GBM_sketch[c].s_j = i;
                        if (k_star == 0){
                            Prune_Module_Flag = true;
                            j_star = argmax(GBM_sketch);
                        }
                    }
                    else if (b_i < GBM_sketch[c].y_j){
                        GBM_sketch[c].y_j = b_i;
                        GBM_sketch[c].s_j = i;
                    }
            }
            else{
                if (GBM_sketch[j_star].y_j < b_i) break;
                else if (b_i < GBM_sketch[c].y_j){
                    GBM_sketch[c].y_j = b_i;
                    GBM_sketch[c].s_j = i;
                    if(c == j_star){
                        j_star = argmax(GBM_sketch);
                    }
                }
            }
        }
    }
    generator_vector.clear();
    return GBM_sketch;
}

int main(){
    const uint32_t n = 10000;
    const uint32_t k = 400;
    const uint32_t round = 1000;
    srand(clock());
    vector<double> TestData;
    for (int i = 0; i<n;i++){
        TestData.emplace_back(rand()/double(RAND_MAX));
    }
    time_t start = clock();

    for (int i = 0; i < round; i++) {
        int seed = i+1;
        Initial(TestData, k);
        vector<GBM_org_var> SketchA = FastGM(TestData, k, seed);
    }
    time_t end = clock();
    cout << "Average Time: " << double(end - start) / CLOCKS_PER_SEC / round << endl;
    return 0;
}