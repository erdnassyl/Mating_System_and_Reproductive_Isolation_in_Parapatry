// Mating System and Reproductive Isolation in Parapatry
// Two-Loci C++ Models with BDMi mutations
// Lyssandre Marchand M2 intership, codes inspired from Lucas Marie-Orleach et al. 2022. 

// Stochastic model
// Secondary Contact Case in Contient-Island Model

#include <iostream>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <algorithm>
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
unsigned N1(0);
unsigned N2(0);

// Mutation parameters (will all be set at 0 because we are fixed)
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
double ha_1(.0);
double sa_1(.0);
double hb_1(.0);
double sb_1(.0);
double epsilon_1_1(.0);
double epsilon_2_1(.0);
double epsilon_3_1(.0);
double epsilon_4_1(.0);
double ha_2(.0);
double sa_2(.0);
double hb_2(.0);
double sb_2(.0);
double epsilon_1_2(.0);
double epsilon_2_2(.0);
double epsilon_3_2(.0);
double epsilon_4_2(.0);

// Recombination rates
double rec_1(.0);
double rec_2(.0);

int span (0);
int interval (0);

int main(int, char* argv[]) {
    
  const gsl_rng_type* T;
  gsl_rng* r;

  // create a generator chosen by the environment variable GSL_RNG_TYPE
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned long)time(0);
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  //////// PARAMETERS

  // use argument passed in command line
  threshold=atoll(argv[1]);
  N_iter=atoll(argv[2]);
  span=atoi(argv[3]);
  interval=atoi(argv[4]); 
  
  N1=atoi(argv[5]);
  N2=atoi(argv[6]);

  self_r_1=atof(argv[7]);
  self_r_2=atof(argv[8]);
  
  mu_Aa_1=mu_Bb_1=atof(argv[9]);
  mu_Aa_2=mu_Bb_2=atof(argv[10]);
  
  ha_1=atof(argv[11]);
  sa_1=atof(argv[12]); 
  hb_1=atof(argv[13]);
  sb_1=atof(argv[14]);
  epsilon_1_1=epsilon_2_1=atof(argv[15]);
  epsilon_3_1=atof(argv[16]);
  epsilon_4_1=atof(argv[17]);
  ha_2=atof(argv[18]);
  sa_2=atof(argv[19]); 
  hb_2=atof(argv[20]);
  sb_2=atof(argv[21]);
  epsilon_1_2=epsilon_2_2=atof(argv[22]);
  epsilon_3_2=atof(argv[23]);
  epsilon_4_2=atof(argv[24]);
  
  rec_1=atof(argv[25]);
  rec_2=atof(argv[26]);
  
  m_h_1=atof(argv[27]);
  m_h_2=atof(argv[28]);
  m_d_1=atof(argv[29]);
  m_d_2=atof(argv[30]);

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
  FITNESS_LANDSCAPE_BM(ha_1, sa_1, hb_1, sb_1, epsilon_1_1, epsilon_2_1, epsilon_3_1, epsilon_4_1, Fitness_1); // Population 1
  FITNESS_LANDSCAPE_BM(ha_2, sa_2, hb_2, sb_2, epsilon_1_2, epsilon_2_2, epsilon_3_2, epsilon_4_2, Fitness_2); // Population 2 
  	
  for (int k(0); k < (int)N_iter; ++k) {
     
    // Condition initialisation
    double dip_FREQ_1[10] ={};
    double dip_FREQ_2[10] = {}; 
    unsigned int dip_IND1[10] = { N1 };
    unsigned int dip_IND2[10] = { N2 };
    double after_repro_1[10] = {};
    double after_repro_2[10] = {};
    double final_1[10] = {};
    double final_2[10] = {};
    double al_FREQ_1[4] = {}; // Stores the new allelic frequences at the end of the cycle
    double al_FREQ_2[4] = {};
    
    unsigned long long int gen(0);
    bool isfin(0);
    
    vector<vector<double>> gen_FREQ_1(0,vector<double>(10));
    vector<vector<double>> gen_FREQ_2(0,vector<double>(10));
   
    //////// LIFE CYCLE 
    
    while (isfin == 0) {
      
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
      REPRODUCTION_POP1_STO(self_r_1, dip_IND1, Me_Mu_Matrix_1, after_repro_1, m_h_1, dip_IND2, Me_Mu_Matrix_2);
      REPRODUCTION_POP2_STO(self_r_2, dip_IND2, Me_Mu_Matrix_2, after_repro_2, m_h_2, dip_IND1, Me_Mu_Matrix_1);
		
	 // Selection
      SELECTION(after_repro_1,Fitness_1,selected_1); // Continent 
      SELECTION(after_repro_2,Fitness_2,selected_2); // Island

      // Migration
      SEED_MIGRATION_POP1(m_d_2, selected_1, selected_2, final_1);
      SEED_MIGRATION_POP2(m_d_1, selected_1, selected_2, final_2);

      // Drift 
      DRIFT(r, N1, final_1, dip_IND1);
      DRIFT(r, N2, final_2, dip_IND2);

      //Stop conditions
        if (gen >= threshold){
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
    ALLELE_FREQ_COMP_STO(dip_IND1, N1, al_FREQ_1);
    ALLELE_FREQ_COMP_STO(dip_IND2, N2, al_FREQ_2);
 
    //Fill in output file
    std::ofstream outfile;
    outfile.open("Output_TL_BM_SC_STO.csv", std::ios_base::app);
    outfile << threshold << "," << N_iter << "," << span << "," << interval << "," << N1 << "," << N2 << "," << self_r_1 << "," << self_r_2 << ",";
    outfile << ha_2 << "," << sa_2 << "," << hb_2 << "," << sb_2 << "," << epsilon_1_2 << "," << epsilon_3_2 << "," << epsilon_4_2 << ","<< rec_2 << ",";
    outfile << m_h_1 << "," << m_h_2 << "," << m_d_1 << "," << m_d_2 << "," << gen;
            
    outfile << "," << al_FREQ_1[0] << "," << al_FREQ_1[1] << "," << al_FREQ_1[2] << "," << al_FREQ_1[3];
    outfile << "," << al_FREQ_2[0] << "," << al_FREQ_2[1] << "," << al_FREQ_2[2] << "," << al_FREQ_2[3];
      for (int i=0; i<(int)gen_FREQ_1.size() && i<=span; ++i) {
        for (int j(0); j<=9; ++j) {
          outfile << "," << gen_FREQ_1[i][j] << "," << gen_FREQ_2[i][j];
        }
      }
    outfile << std::endl;						
  }
}
