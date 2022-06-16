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
    int z_i;
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
        CUST.z_i = 0;
        CUSTVec[i] = CUST;
        for(int j=0;j<k;j++) PermutationVec[i*k+j] = j;
    }

    for(uint32_t i=0;i<k;i++) {
        GBM_VAR.y_j = -1;
        GBM_VAR.s_j = -1;
        GBM_sketch[i] = GBM_VAR;
    }
}

void AscEXP(int i, int k, int seed){
    double x=0;
    int j=0, z=0;
    z = CUSTVec[i].z_i+1;
    //default_random_engine generator(i << 12 + seed);
    //generator_vector.push_back(generator);
    default_random_engine generator = generator_vector[i];
    uniform_int_distribution<unsigned> distribution1(1, RAND_MAX-1);
    unsigned r = distribution1(generator);
    double u = (r) / double(RAND_MAX);
    j = (r % (k-z+1))+z;
    x = -(log(u))/(k+1-z);


    int temp = PermutationVec[i * k + z-1];
    PermutationVec[i * k + z-1] = PermutationVec[i * k + j];
    PermutationVec[i * k + j] = temp;

    CUSTVec[i].x_i += x;
    CUSTVec[i].z_i = z;
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
    int j_star=0;
    while (k_star != 0) {
        delta = delta + k;
        for (int i=0;i<data.size();i++){
            double v_i = data[i];
            if(CUSTVec[i].z_i == 0) {
                default_random_engine generator((i << 9) + seed);
                generator_vector.push_back(generator);
            }
            if(v_i == 0) continue;
            uint64_t r_i = delta * v_i / v_sum;

            while(CUSTVec[i].z_i < r_i){
                AscEXP(i, k, seed);
                double b_i = CUSTVec[i].x_i /(v_i);
                int c = CUSTVec[i].c;
                if (GBM_sketch[c].y_j < 0 ){
                    GBM_sketch[c].y_j = b_i;
                    k_star = k_star - 1;
                    GBM_sketch[c].s_j = i;
                }
                else if (b_i < GBM_sketch[c].y_j){
                    GBM_sketch[c].y_j = b_i;
                    GBM_sketch[c].s_j = i;
                }
            }
        }
    }
//FastPrune
    j_star = argmax(GBM_sketch);
    for(int i=0; i!=data.size();i++){
        double v_i = data[i];

        while(CUSTVec[i].z_i <= k){
            AscEXP(i, k, seed);
            int c = CUSTVec[i].c;

            double b_i = CUSTVec[i].x_i /(v_i);
            if (GBM_sketch[j_star].y_j < b_i or CUSTVec[i].z_i == k) break;
            else if (b_i < GBM_sketch[c].y_j){
                GBM_sketch[c].y_j = b_i;
                GBM_sketch[c].s_j = i;
                if(c == j_star){
                    j_star = argmax(GBM_sketch);
                }
            }
        }
    }
    generator_vector.clear();
    return GBM_sketch;
}


vector<GBM_org_var> BBM_Basic(const vector<double> data, const uint32_t k, int seed){
    //vector<double> yVec(k, -1);
    //vector<int> Sketch(k, -1);
    for(int i=0; i<data.size(); i++){
        srand((i<<8) + seed);
        for (int j = 0; j < k; j++) {
            double u = rand() / double(RAND_MAX);
            double exp = -log(u) / data[i];
            if(GBM_sketch[j].y_j == -1 or exp < GBM_sketch[j].y_j) {
                GBM_sketch[j].y_j = exp;
                GBM_sketch[j].s_j = i;
            }
        }
    }
    return GBM_sketch;
}

int main(){
    const uint32_t n = 10000;
    const uint32_t k = 400;
    const uint32_t round = 100;
    srand(clock());
    vector<double> TestData;
    double ture_Lambda = 0;
    double bias = 0;
    for (int i = 0; i<n;i++){
        double w_i = rand()/double(RAND_MAX);
        TestData.emplace_back(w_i);
        ture_Lambda = ture_Lambda + w_i;
    }
    cout << "Ture weighted cardinality: " << ture_Lambda << endl;
    time_t start = clock();

    for (int i = 0; i < round; i++) {
        int seed = i+1;
        double hat_Lambda = 0;
        double G_k =0;
        Initial(TestData, k);
        vector<GBM_org_var> SketchA = FastGM(TestData, k, seed);
//        vector<GBM_org_var> SketchA = BBM_Basic(TestData, k, seed);
        for(int m =0; m<k; m++){
            G_k = G_k + SketchA[m].y_j;
        }
        hat_Lambda = (k)/G_k;
        //cout << "Estimation weighted cardinality: " << hat_Lambda << endl;
        bias = bias + (ture_Lambda - hat_Lambda);

    }
    time_t end = clock();
    cout << "Average Time: " << double(end - start) / CLOCKS_PER_SEC / round << endl;
    cout << "Average Bias: " << bias/ round << endl;
    return 0;
}