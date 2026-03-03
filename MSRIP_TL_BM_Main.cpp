// Mating System and Reproductive Isolation in Parapatry
// Two-Loci C++ Models
// BDMi mutations � MAIN
// Lucas Marie-Orleach

#include <iostream>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>
#include <vector>
#include "MSRIP_TL_MyFunctions.h"

using namespace std;

// Define parameters
unsigned long long int threshold(0);
unsigned long long int N_iter(0);
unsigned int N_1(0);
unsigned int N_2(0);

double mu_Aa_1(.0);
double mu_aA_1(.0);
double mu_Bb_1(.0);
double mu_bB_1(.0);
double mu_Aa_2(.0);
double mu_aA_2(.0);
double mu_Bb_2(.0);
double mu_bB_2(.0);

double self_r_1(.0);
double self_r_2(.0);

double m_h_1(.0);
double m_h_2(.0);
double m_d_1(.0);
double m_d_2(.0);

double Fitness_1[10] = { 1,1,1,1,1,1,1,1,1,1 };
double Fitness_2[10] = { 1,1,1,1,1,1,1,1,1,1 };
double ha_1(.0);
double sa_1(.0);
double hb_1(.0); 
double sb_1(.0);
double ha_2(.0);
double sa_2(.0);
double hb_2(.0); 
double sb_2(.0);

double rec_1(.0);
double rec_2(.0);

double s_B_1(0);
double h_B_1(0);
double k_B_1(0);
double s_B_2(0);
double h_B_2(0);
double k_B_2(0);


int span (0);
int interval (0);

int main(int, char* argv[]) {

  // const gsl_rng_type* T;
  // gsl_rng* r;

  // create a generator chosen by the environment variable GSL_RNG_TYPE
  // gsl_rng_env_setup();
  // gsl_rng_default_seed = (unsigned long)time(0);
  // T = gsl_rng_default;
  // r = gsl_rng_alloc(T);
 
  // use argument passed in command line
  threshold=atoll(argv[1]);
  N_iter=atoll(argv[2]);
  span=atoi(argv[3]);
  interval=atoi(argv[4]);  
  
  N_1=atoi(argv[5]);
  N_2=atoi(argv[6]);

  self_r_1=atof(argv[7]);
  self_r_2=atof(argv[8]);
  
  mu_Aa_1=mu_Bb_1=atof(argv[9]);
  mu_Aa_2=mu_Bb_2=atof(argv[10]);
  
  ha_1=atof(argv[11]);
  sa_1=atof(argv[12]);
  ha_2=atof(argv[13]);
  sa_2=atof(argv[14]);
  
  rec_1=atof(argv[15]);
  rec_2=atof(argv[16]);
  
  h_B_1=k_B_1=atof(argv[17]);
  s_B_1=atof(argv[18])*pow(10,atof(argv[19]));
  h_B_2=k_B_2=atof(argv[20]);
  s_B_2=atof(argv[21])*pow(10,atof(argv[22]));

  m_h_1=atof(argv[23]);
  m_h_2=atof(argv[24]);
  m_d_1=atof(argv[25]);
  m_d_2=atof(argv[26]);
            
  // Compute the Meiosis_Mutation Matrix for population 1 
  double Me_Matrix_1[10][4];
  double Mu_Matrix_1[4][4];
  double Me_Mu_Matrix_1[10][4];
  
  Me_MATRIX_COMP(rec_1, Me_Matrix_1);
  Mu_MATRIX_COMP(mu_Aa_1, mu_aA_1, mu_Bb_1, mu_bB_1, Mu_Matrix_1);
  Me_Mu_MATRIX_COMP(Me_Matrix_1, Mu_Matrix_1, Me_Mu_Matrix_1);

  // Compute the Meiosis_Mutation Matrix for population 2 
  double Me_Matrix_2[10][4];
  double Mu_Matrix_2[4][4];
  double Me_Mu_Matrix_2[10][4];
  
  Me_MATRIX_COMP(rec_2, Me_Matrix_2);
  Mu_MATRIX_COMP(mu_Aa_2, mu_aA_2, mu_Bb_2, mu_bB_2, Mu_Matrix_2);
  Me_Mu_MATRIX_COMP(Me_Matrix_2, Mu_Matrix_2, Me_Mu_Matrix_2);
  		
  //Compute Fitness landsacpe
  FITNESS_LANDSCAPE_BM(sa_1, ha_1, sb_1, hb_1, s_B_1, h_B_1, k_B_1, Fitness_1);
  FITNESS_LANDSCAPE_BM(sa_2, ha_2, sb_2, hb_2, s_B_2, h_B_2, k_B_2, Fitness_2);
  	
  for (int k(0); k < (int)N_iter; ++k) {
     
    //condition initialisation
    double dip_FREQ_1[10] = {};
    double dip_FREQ_2[10] = {};
    double final_1[10] = {};
    double final_2[10] = {};
    unsigned int dip_IND_1[10] = { N_1 };
    unsigned int dip_IND_2[10] = { N_2 };
    double al_FREQ_1[4] = {};
    double al_FREQ_2[4] = {};
      
    unsigned long long int gen(0);
    bool isfin(0);
    
    vector<vector<double>> gen_FREQ_1(0,vector<double>(10));
    vector<vector<double>> gen_FREQ_2(0,vector<double>(10));
   
    //life cycle
    while (isfin == 0) {
      // Remise à zéro des fréquences temporaires
      for(int i=0; i<10; ++i){ dip_FREQ_1[i]=0; dip_FREQ_2[i]=0; }

      //record genotypic frequencies on the first gen
      if(span!=0 && gen==0) {
        gen_FREQ_1.push_back(vector<double>(10));
        gen_FREQ_2.push_back(vector<double>(10));
        for (int i(0); i<=9; ++i) {
          gen_FREQ_1[ceil(gen/interval)][i]=(double)dip_IND_1[i]/N_1;
          gen_FREQ_2[ceil(gen/interval)][i]=(double)dip_IND_2[i]/N_2;
        }
      }
   
      // Reproduction Pop 1 (donne dip_FREQ_1)
      REPRODUCTION(self_r_1, dip_IND_1, Me_Mu_Matrix_1, dip_FREQ_1, m_h_1, dip_IND_2, Me_Mu_Matrix_2);
      // Reproduction Pop 2 (donne dip_FREQ_2)
      REPRODUCTION(self_r_2, dip_IND_2, Me_Mu_Matrix_2, dip_FREQ_2, m_h_2, dip_IND_1, Me_Mu_Matrix_1);

      SEED_MIGRATION(m_d_1, m_d_2, Fitness_1, Fitness_2, dip_FREQ_1, dip_FREQ_2, final_1, final_2);

      for(int i=0; i<10; ++i) {
        dip_IND_1[i] = (unsigned int)(final_1[i] * N_1);
        dip_IND_2[i] = (unsigned int)(final_2[i] * N_2);
      }

      //Stop conditions						
      if (dip_IND_1[7]+dip_IND_1[8]+dip_IND_1[9]==N_1 || dip_IND_1[4]+dip_IND_1[6]+dip_IND_1[9]==N_1 ||
          dip_IND_2[7]+dip_IND_2[8]+dip_IND_2[9]==N_2 || dip_IND_2[4]+dip_IND_2[6]+dip_IND_2[9]==N_2 ||
          gen>=threshold) { isfin = 1; }

      //record genotypic frequencies every interval
      if((gen!=0 && span!=0 && gen%interval==0) || isfin == 1) {
        gen_FREQ_1.push_back(vector<double>(10));
        gen_FREQ_2.push_back(vector<double>(10));
        for (int i(0); i<=9; ++i) {
          gen_FREQ_1[(int)gen_FREQ_1.size()-1][i]=(double)dip_IND_1[i]/N_1;
          gen_FREQ_2[(int)gen_FREQ_2.size()-1][i]=(double)dip_IND_2[i]/N_2;
        }
      }
      
      //add generation
      ++gen; 
    }
   
    //compute final allele frequencies
    ALLELE_FREQ_COMP(final_1, al_FREQ_1);
    ALLELE_FREQ_COMP(final_2, al_FREQ_2);
 
    //Fill in output file
    std::ofstream outfile;
    outfile.open("Output_TL_BM.csv", std::ios_base::app);
    outfile << threshold << "," << N_iter << "," << span << "," << interval << "," << N_1 << "," << N_2 << "," << self_r_1 << "," << self_r_2 << ",";
    outfile << mu_Aa_1 << "," << mu_Aa_2 << ","<< mu_aA_1 << ","<< mu_aA_2 << ","<< mu_Bb_1 << ","<< mu_bB_2 << ","<< ha_1 << "," << sa_1 << "," << ha_2 << "," << sa_2 << ","<< rec_1 << "," << rec_2 << ","<< h_B_1 << "," ;
    outfile << k_B_1 << "," << s_B_1 << ","<< h_B_2 << "," << k_B_2 << "," << s_B_2 << ","<< m_h_1 << "," << m_h_2 << "," << m_d_1 << "," << m_d_2 << ","<< gen;
            
    if(gen != threshold+1) {
      outfile << "," << al_FREQ_1[0] << "," << al_FREQ_1[1] << "," << al_FREQ_1[2] << "," << al_FREQ_1[3];
      outfile << "," << al_FREQ_2[0] << "," << al_FREQ_2[1] << "," << al_FREQ_2[2] << "," << al_FREQ_2[3];
      for (int i=(int)gen_FREQ_1.size()-1; (i>=0) && (((int)gen_FREQ_1.size()-1)-i<=span); --i) {
        for (int j(0); j<=9; ++j) {
          outfile << "," << gen_FREQ_1[i][j] << "," << gen_FREQ_2[i][j];
        }
      }
    } outfile << std::endl;						
  }
}