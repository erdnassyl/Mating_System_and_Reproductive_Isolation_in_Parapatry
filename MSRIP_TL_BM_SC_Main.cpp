// Mating System and Reproductive Isolation in Parapatry
// Two-Loci C++ Models with BDMi mutations
// Lyssandre Marchand M2 intership, codes inspired from Lucas Marie-Orleach et al. 2022. 

// Deterministic model with infinite population 
// Secondary Contact Case

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

// Mutation parameters
double mu_Aa_1(.0);
double mu_aA_1(.0);
double mu_Bb_1(.0);
double mu_bB_1(.0);
double mu_Aa_2(.0);
double mu_aA_2(.0);
double mu_Bb_2(.0);
double mu_bB_2(.0);

// Selfing rates parameters
double self_r_1(.0);
double self_r_2(.0);

// Haploid and diploid migration rates 
double m_h_1(.0);
double m_h_2(.0);
double m_d_1(.0);
double m_d_2(.0);

// Fitness initialization, selection and dominance rates 
double Fitness_1[10] = { 1,1,1,1,1,1,1,1,1,1 };
double Fitness_2[10] = { 1,1,1,1,1,1,1,1,1,1 };
double alpha_1(.0);
double beta_1(.0);
double gamma_1(.0);
double alpha_2(.0);
double beta_2(.0);
double gamma_2(.0); 

// Recombination rates
double rec_1(.0);
double rec_2(.0);

// Selection and dominance for BDMi mutations
double s_B_1(0);
double h_B_1(0);
double k_B_1(0);
double s_B_2(0);
double h_B_2(0);
double k_B_2(0);

int span (0);
int interval (0);

int main(int, char* argv[]) {

  //////// PARAMETERS

  // use argument passed in command line
  threshold=atoll(argv[1]);
  N_iter=atoll(argv[2]);
  span=atoi(argv[3]);
  interval=atoi(argv[4]);  

  self_r_1=atof(argv[5]);
  self_r_2=atof(argv[6]);
  
  mu_Aa_1=mu_Bb_1=atof(argv[7]);
  mu_Aa_2=mu_Bb_2=atof(argv[8]);
  
  alpha_1=atof(argv[9]);
  beta_1=atof(argv[10]); 
  gamma_1=atof(argv[11]);
  alpha_2=atof(argv[12]);
  beta_2=atof(argv[13]);
  gamma_2=atof(argv[14]);
  
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

  // Fixation threshold 
  double epsilon=1e-3;

  //////// MEIOSIS_MUTATION MATRICES AND FITNESS
  
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
  		
  // Compute Fitness landsacpe
  FITNESS_LANDSCAPE(alpha_1, beta_1, gamma_1, Fitness_1);
  FITNESS_LANDSCAPE(alpha_2, beta_2, gamma_2, Fitness_2);
  	
  for (int k(0); k < (int)N_iter; ++k) {
     
    // Condition initialisation
    double dip_FREQ_1[10] ={1.0,0,0,0,0,0,0,0,0,0}; // AABB dans pop 1
    double dip_FREQ_2[10] = {0,0,0,0,0,0,0,0,0,1.0}; // aabb dans pop 2 

    double after_repro_1[10] = {};
    double after_repro_2[10] = {};
    double final_1[10] = {};
    double final_2[10] = {};
    
    double pre_al_FREQ_1[4] = {}; // Stores the previous allelic frequences at the beginning of the cycle
    double pre_al_FREQ_2[4] = {};
    double al_FREQ_1[4] = {}; // Stores the new allelic frequences at the end of the cycle
    double al_FREQ_2[4] = {};
    
    unsigned long long int gen(0);
    bool isfin(0);
    
    vector<vector<double>> gen_FREQ_1(0,vector<double>(10));
    vector<vector<double>> gen_FREQ_2(0,vector<double>(10));
   
    //////// LIFE CYCLE 
    
    while (isfin == 0) {

      // Compute the allelic frequencies before the life cycle 
      // Will be used to compute the stop conditions
      ALLELE_FREQ_COMP(dip_FREQ_1, pre_al_FREQ_1);
      ALLELE_FREQ_COMP(dip_FREQ_2, pre_al_FREQ_2);
      
      // Remise à zéro des fréquences temporaires
      for(int i=0; i<10; ++i){ 
        after_repro_1[i]=0; 
        after_repro_2[i]=0; 
        final_1[i]=0; 
        final_2[i]=0; 
      }

      // Record genotypic frequencies on the first gen
      if(span!=0 && gen==0) {
        gen_FREQ_1.push_back(vector<double>(10));
        gen_FREQ_2.push_back(vector<double>(10));
        for (int i(0); i<=9; ++i) {
          gen_FREQ_1[ceil(gen/interval)][i]=dip_FREQ_1[i];
          gen_FREQ_2[ceil(gen/interval)][i]=dip_FREQ_2[i];
        }
      }
   
      // Reproduction 
      REPRODUCTION_POP1(self_r_1, dip_FREQ_1, Me_Mu_Matrix_1, after_repro_1, m_h_1, dip_FREQ_2, Me_Mu_Matrix_2);
      REPRODUCTION_POP2(self_r_2, dip_FREQ_2, Me_Mu_Matrix_2, after_repro_2, m_h_2, dip_FREQ_1, Me_Mu_Matrix_1);

      // Migration
      SEED_MIGRATION(m_d_1, m_d_2, Fitness_1, Fitness_2, after_repro_1, after_repro_2, final_1, final_2);

      // Storing the new allelic frequences
      for(int i=0; i<10; ++i) {
        dip_FREQ_1[i] = final_1[i];
        dip_FREQ_2[i] = final_2[i];
      }

      // Compute the new allelic frequences and the delta 
      ALLELE_FREQ_COMP(dip_FREQ_1, al_FREQ_1);
      ALLELE_FREQ_COMP(dip_FREQ_2, al_FREQ_2);

      double delta = 0.0;
      for (int i = 0; i < 3; ++i){
        double d1 = std::abs(al_FREQ_1[i] - pre_al_FREQ_1[i]); // deltas computations for pop 1
        double d2 = std::abs(al_FREQ_2[i] - pre_al_FREQ_2[i]); // deltas computations for pop 2

        // We take the highest value
        if (d1 > delta) delta = d1;
        if (d2 > delta) delta = d2;
      }

      //Stop conditions
        if (delta < epsilon || gen >= threshold){
          isfin = 1;
        }

      //record genotypic frequencies every interval
      if((gen!=0 && span!=0 && gen%interval==0) || isfin == 1) {
        gen_FREQ_1.push_back(vector<double>(10));
        gen_FREQ_2.push_back(vector<double>(10));
        for (int i(0); i<=9; ++i) {
          gen_FREQ_1[(int)gen_FREQ_1.size()-1][i]=dip_FREQ_1[i];
          gen_FREQ_2[(int)gen_FREQ_2.size()-1][i]=dip_FREQ_2[i];
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
    outfile.open("Output_TL_BM_SC.csv", std::ios_base::app);
    outfile << threshold << "," << N_iter << "," << span << "," << interval << "," << self_r_1 << "," << self_r_2 << ",";
    outfile << mu_Aa_1 << "," << mu_aA_1 << "," << mu_Bb_1 << ","<< mu_bB_1 <<","<< mu_Aa_2 << "," << mu_aA_2 << "," << mu_Bb_2 << ","<< mu_bB_2 << ","<< alpha_1 << "," << beta_1 << "," << gamma_1 << "," << alpha_2 << "," << beta_2 << "," << gamma_2 << "," << rec_1 << "," << rec_2 << ",";
    outfile << ","<< m_h_1 << "," << m_h_2 << "," << m_d_1 << "," << m_d_2 << ","<< gen;
            
    
    outfile << "," << al_FREQ_1[0] << "," << al_FREQ_1[1] << "," << al_FREQ_1[2] << "," << al_FREQ_1[3];
    outfile << "," << al_FREQ_2[0] << "," << al_FREQ_2[1] << "," << al_FREQ_2[2] << "," << al_FREQ_2[3];
      for (int i=(int)gen_FREQ_1.size()-1; (i>=0) && (((int)gen_FREQ_1.size()-1)-i<=span); --i) {
        for (int j(0); j<=9; ++j) {
          outfile << "," << gen_FREQ_1[i][j] << "," << gen_FREQ_2[i][j];
        }
      }
    outfile << std::endl;						
  }
}
