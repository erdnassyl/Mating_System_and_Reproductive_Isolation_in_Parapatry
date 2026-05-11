// Mating System and Reproductive Isolation in Parapatry
// Two-Loci C++ Models with BDMi mutations
// Lyssandre Marchand M2 intership, codes inspired from Lucas Marie-Orleach et al. 2022. 

// Functions 

#pragma once
#include <iostream>

// Compute Meiosis Matrix given recombination rate
void Me_MATRIX_COMP(const double rec, double Me_Matrix[][4])
{
	double MAT[10][4] = {
		{1,  0,  0,  0},
		{.5, .5, 0,  0},
		{.5, 0,  .5, 0},
		{.5 - (rec / 2),rec / 2,rec / 2,.5 - (rec / 2)},
		{0,  1,  0,  0},
		{rec / 2,.5 - (rec / 2),.5 - (rec / 2),rec / 2},
		{0,  .5, 0,  .5},
		{0,  0,  1,  0},
		{0,  0,  .5,.5},
		{0,  0,  0,  1}
	};

	for (int i(0); i < 10; ++i) { 
		for (int j(0); j < 4; ++j) { 
			Me_Matrix[i][j] = MAT[i][j]; 
		} 
	}
}

// Compute Mutation Matrix given mutation rates
void Mu_MATRIX_COMP(const double Mu_Aa, const double mu_aA, const double mu_Bb, const double mu_bB, double Mu_Matrix[][4])
{
	double MAT[4][4] = {
		 {(1 - Mu_Aa) * (1 - mu_Bb),(1 - Mu_Aa) * mu_Bb,      Mu_Aa * (1 - mu_Bb),	     Mu_Aa * mu_Bb},
		 {(1 - Mu_Aa) * mu_bB,	     (1 - Mu_Aa) * (1 - mu_bB),Mu_Aa * mu_bB,			 Mu_Aa * (1 - mu_bB)},
		 {mu_aA * (1 - mu_Bb),	     mu_aA * mu_Bb, 		   (1 - mu_aA) * (1 - mu_Bb),(1 - mu_aA) * mu_Bb},
		 {mu_aA * mu_bB,			 mu_aA * (1 - mu_bB),	   (1 - mu_aA) * mu_bB,	     (1 - mu_aA) * (1 - mu_bB)}
	};

	for (int i(0); i < 4; ++i) {
		for (int j(0); j < 4; ++j) {
			Mu_Matrix[i][j] = MAT[i][j];
		}
	}
}

// Compute Meiose Matrix Matrix given Meiose Matrix and Mutation Matrix
void Me_Mu_MATRIX_COMP(const double Me_Matrix[][4], const double Mu_Matrix[][4], double Me_Mu_Matrix[][4])
{
	double MAT[10][4] = {};
	for (int i(0); i < 10; ++i) {
		for (int j(0); j < 4; ++j) {
			for (int k(0); k < 4; ++k) {
				//cout << Me_Matrix[i][k] << ", " << Mu_Matrix[k][j] << ", " << Me_Matrix[i][k] * Mu_Matrix[k][j] << endl;
				MAT[i][j] += Me_Matrix[i][k] * Mu_Matrix[k][j];
			}
		}
	}

	for (int i(0); i < 10; ++i) {
		for (int j(0); j < 4; ++j) {
			Me_Mu_Matrix[i][j] = MAT[i][j];
		}
	}
}

// Compute Fitness Landscape given coefficients of dominance and strenght of selection
void FITNESS_LANDSCAPE_BM(const double ha, const double sa, const double hb, const double sb, const double epsilon_1, const double epsilon_2, const double epsilon_3, const double epsilon_4, double* Fitness)
{
	Fitness[0] = 1;
	Fitness[1] = (1 + (hb * sb));
	Fitness[2] = (1 + (ha * sa));
	Fitness[3] = (1 + (ha * sa)) * (1 + (hb * sb)) * (1 + epsilon_1);
	Fitness[4] = (1 + sb);
	Fitness[5] = (1 + (ha * sa)) * (1 + (hb * sb)) * (1 + epsilon_2);
	Fitness[6] = (1 + (ha * sa)) * (1 + sb) * (1 + epsilon_1);
	Fitness[7] = (1 + sa);
	Fitness[8] = (1 + sa) * (1 + (hb * sb)) * (1 + epsilon_3);
	Fitness[9] = (1 + sa) * (1 + sb) * (1 + epsilon_4);
}

// Compute gamete haplotypes given adult genotypes and the Meiose Mutation Matrix
void GAMETE_PROD(const double dip_IND[10], const double Me_Mu_Matrix[][4], double* hap_FREQ)
{
	for (int i(0); i < 10; ++i) {
		for (int j(0); j < 4; ++j) {
			hap_FREQ[j] += dip_IND[i] * Me_Mu_Matrix[i][j];
		}
	}
}
void GAMETE_PROD_STO(const unsigned int dip_IND[10], const double Me_Mu_Matrix[][4], double* hap_FREQ)
{
	for (int i(0); i < 10; ++i) {
		for (int j(0); j < 4; ++j) {
			hap_FREQ[j] += dip_IND[i] * Me_Mu_Matrix[i][j];
		}
	}
}
// Define Reproduction without pollen migration function given selfing rate, genotypes of adult individuals, the Meiosis_Mutation matrix, and the Seeds's genotypic frequency 
void REPRODUCTION_POP1(const double self_r_1, const double dip_IND_1[10], const double Me_Mu_Matrix_1[][4], double* dip_FREQ_1,const double m_h_1, const double dip_IND_2[10], const double Me_Mu_Matrix_2[][4])	
{	// Seed genotypes produced through selfing
	double self_dip[10] = {};
	for (int i(0); i < 10; ++i) {
		self_dip[0] += dip_IND_1[i] * Me_Mu_Matrix_1[i][0] * Me_Mu_Matrix_1[i][0];
		self_dip[1] += dip_IND_1[i] * Me_Mu_Matrix_1[i][0] * Me_Mu_Matrix_1[i][1] * 2;
		self_dip[2] += dip_IND_1[i] * Me_Mu_Matrix_1[i][0] * Me_Mu_Matrix_1[i][2] * 2;
		self_dip[3] += dip_IND_1[i] * Me_Mu_Matrix_1[i][0] * Me_Mu_Matrix_1[i][3] * 2;
		self_dip[4] += dip_IND_1[i] * Me_Mu_Matrix_1[i][1] * Me_Mu_Matrix_1[i][1];
		self_dip[5] += dip_IND_1[i] * Me_Mu_Matrix_1[i][1] * Me_Mu_Matrix_1[i][2] * 2;
		self_dip[6] += dip_IND_1[i] * Me_Mu_Matrix_1[i][1] * Me_Mu_Matrix_1[i][3] * 2;
		self_dip[7] += dip_IND_1[i] * Me_Mu_Matrix_1[i][2] * Me_Mu_Matrix_1[i][2];
		self_dip[8] += dip_IND_1[i] * Me_Mu_Matrix_1[i][2] * Me_Mu_Matrix_1[i][3] * 2;
		self_dip[9] += dip_IND_1[i] * Me_Mu_Matrix_1[i][3] * Me_Mu_Matrix_1[i][3];
	}

	double self_dip_SUM(0.0);
	for (int i(0); i < 10; ++i) { self_dip_SUM += self_dip[i]; }
	double self_dip_FREQ[10];
	for (int i(0); i < 10; ++i) { self_dip_FREQ[i] = self_dip[i] / self_dip_SUM; }

	// Seed genotypes produced through outcrossing local
	double local_ovule[4] = {};
	double local_pollen[4] = {};
	GAMETE_PROD(dip_IND_1, Me_Mu_Matrix_1, local_ovule);
	GAMETE_PROD(dip_IND_1, Me_Mu_Matrix_1, local_pollen);

	double out_local_dip_1[10] = {};
	out_local_dip_1[0] = local_ovule[0] * local_pollen[0];
	out_local_dip_1[1] = local_ovule[0] * local_pollen[1] + local_ovule[1] * local_pollen[0]; 
	out_local_dip_1[2] = local_ovule[0] * local_pollen[2] + local_ovule[2] * local_pollen[0];
	out_local_dip_1[3] = local_ovule[0] * local_pollen[3] + local_ovule[3] * local_pollen[0];
	out_local_dip_1[4] = local_ovule[1] * local_pollen[1];
	out_local_dip_1[5] = local_ovule[1] * local_pollen[2] + local_ovule[2] * local_pollen[1];
	out_local_dip_1[6] = local_ovule[1] * local_pollen[3] + local_ovule[3] * local_pollen[1];
	out_local_dip_1[7] = local_ovule[2] * local_pollen[2];
	out_local_dip_1[8] = local_ovule[2] * local_pollen[3] + local_ovule[3] * local_pollen[2];
	out_local_dip_1[9] = local_ovule[3] * local_pollen[3];

	double out_local_SUM_1(0.0);
	for (int i(0); i < 10; ++i) { out_local_SUM_1 += out_local_dip_1[i]; }
	double out_local_FREQ_1[10];
	for (int i(0); i < 10; ++i) { out_local_FREQ_1[i] = out_local_dip_1[i] / out_local_SUM_1; }

	// Seed genotypes produced through outcrossing immigrants
	double immigrant_pollen[4] = {};
	GAMETE_PROD(dip_IND_2, Me_Mu_Matrix_2,immigrant_pollen);
	
	double out_immigrant_dip_1[10] = {};
	out_immigrant_dip_1[0] = local_ovule[0] * immigrant_pollen[0];
	out_immigrant_dip_1[1] = local_ovule[0] * immigrant_pollen[1] + local_ovule[1] * immigrant_pollen[0] ;
	out_immigrant_dip_1[2] = local_ovule[0] * immigrant_pollen[2] + local_ovule[2] * immigrant_pollen[0];
	out_immigrant_dip_1[3] = local_ovule[0] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[0];
	out_immigrant_dip_1[4] = local_ovule[1] * immigrant_pollen[1];
	out_immigrant_dip_1[5] = local_ovule[1] * immigrant_pollen[2] + local_ovule[2] * immigrant_pollen[1];
	out_immigrant_dip_1[6] = local_ovule[1] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[1];
	out_immigrant_dip_1[7] = local_ovule[2] * immigrant_pollen[2];
	out_immigrant_dip_1[8] = local_ovule[2] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[2];
	out_immigrant_dip_1[9] = local_ovule[3] * immigrant_pollen[3];

	double out_im_SUM_1(0.0);
	for (int i(0); i < 10; ++i) { out_im_SUM_1 += out_immigrant_dip_1[i]; }
	double out_im_FREQ_1[10];
	for (int i(0); i < 10; ++i) { out_im_FREQ_1[i] = out_immigrant_dip_1[i] / out_im_SUM_1; }

	// Total seed Genotypes
	for (int i(0); i < 10; ++i) {
		dip_FREQ_1[i] =((self_r_1 * self_dip_FREQ[i]) + ((1 - self_r_1) * (1 - m_h_1) * out_local_FREQ_1[i]) + 
		((1 - self_r_1) * m_h_1 * out_im_FREQ_1[i]));
	}
}
void REPRODUCTION_POP2(const double self_r_2, const double dip_IND_2[10], const double Me_Mu_Matrix_2[][4], double* dip_FREQ_2,const double m_h_2, const double dip_IND_1[10], const double Me_Mu_Matrix_1[][4])
{	// Seed genotypes produced through selfing
	double self_dip[10] = {};
	for (int i(0); i < 10; ++i) {
		self_dip[0] += dip_IND_2[i] * Me_Mu_Matrix_2[i][0] * Me_Mu_Matrix_2[i][0];
		self_dip[1] += dip_IND_2[i] * Me_Mu_Matrix_2[i][0] * Me_Mu_Matrix_2[i][1] * 2;
		self_dip[2] += dip_IND_2[i] * Me_Mu_Matrix_2[i][0] * Me_Mu_Matrix_2[i][2] * 2;
		self_dip[3] += dip_IND_2[i] * Me_Mu_Matrix_2[i][0] * Me_Mu_Matrix_2[i][3] * 2;
		self_dip[4] += dip_IND_2[i] * Me_Mu_Matrix_2[i][1] * Me_Mu_Matrix_2[i][1];
		self_dip[5] += dip_IND_2[i] * Me_Mu_Matrix_2[i][1] * Me_Mu_Matrix_2[i][2] * 2;
		self_dip[6] += dip_IND_2[i] * Me_Mu_Matrix_2[i][1] * Me_Mu_Matrix_2[i][3] * 2;
		self_dip[7] += dip_IND_2[i] * Me_Mu_Matrix_2[i][2] * Me_Mu_Matrix_2[i][2];
		self_dip[8] += dip_IND_2[i] * Me_Mu_Matrix_2[i][2] * Me_Mu_Matrix_2[i][3] * 2;
		self_dip[9] += dip_IND_2[i] * Me_Mu_Matrix_2[i][3] * Me_Mu_Matrix_2[i][3];
	}

	double self_dip_SUM(0.0);
	for (int i(0); i < 10; ++i) { self_dip_SUM += self_dip[i]; }
	double self_dip_FREQ[10];
	for (int i(0); i < 10; ++i) { self_dip_FREQ[i] = self_dip[i] / self_dip_SUM; }

	// Seed genotypes produced through outcrossing local
	double local_ovule[4] = {};
	double local_pollen[4] = {};
	GAMETE_PROD(dip_IND_2, Me_Mu_Matrix_2, local_ovule);
	GAMETE_PROD(dip_IND_2, Me_Mu_Matrix_2, local_pollen);

	double out_local_dip_2[10] = {};
	out_local_dip_2[0] = local_ovule[0] * local_pollen[0];
	out_local_dip_2[1] = local_ovule[0] * local_pollen[1] + local_ovule[1] * local_pollen[0]; 
	out_local_dip_2[2] = local_ovule[0] * local_pollen[2] + local_ovule[2] * local_pollen[0];
	out_local_dip_2[3] = local_ovule[0] * local_pollen[3] + local_ovule[3] * local_pollen[0];
	out_local_dip_2[4] = local_ovule[1] * local_pollen[1];
	out_local_dip_2[5] = local_ovule[1] * local_pollen[2] + local_ovule[2] * local_pollen[1];
	out_local_dip_2[6] = local_ovule[1] * local_pollen[3] + local_ovule[3] * local_pollen[1];
	out_local_dip_2[7] = local_ovule[2] * local_pollen[2];
	out_local_dip_2[8] = local_ovule[2] * local_pollen[3] + local_ovule[3] * local_pollen[2];
	out_local_dip_2[9] = local_ovule[3] * local_pollen[3];

	double out_local_SUM_2(0.0);
	for (int i(0); i < 10; ++i) { out_local_SUM_2 += out_local_dip_2[i]; }
	double out_local_FREQ_2[10];
	for (int i(0); i < 10; ++i) { out_local_FREQ_2[i] = out_local_dip_2[i] / out_local_SUM_2; }

	// Seed genotypes produced through outcrossing immigrants
	double immigrant_pollen[4] = {};
	GAMETE_PROD(dip_IND_1, Me_Mu_Matrix_1,immigrant_pollen);
	
	double out_immigrant_dip_2[10] = {};
	out_immigrant_dip_2[0] = local_ovule[0] * immigrant_pollen[0];
	out_immigrant_dip_2[1] = local_ovule[0] * immigrant_pollen[1] + local_ovule[1] * immigrant_pollen[0] ;
	out_immigrant_dip_2[2] = local_ovule[0] * immigrant_pollen[2] + local_ovule[2] * immigrant_pollen[0];
	out_immigrant_dip_2[3] = local_ovule[0] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[0];
	out_immigrant_dip_2[4] = local_ovule[1] * immigrant_pollen[1];
	out_immigrant_dip_2[5] = local_ovule[1] * immigrant_pollen[2] + local_ovule[2] * immigrant_pollen[1];
	out_immigrant_dip_2[6] = local_ovule[1] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[1];
	out_immigrant_dip_2[7] = local_ovule[2] * immigrant_pollen[2];
	out_immigrant_dip_2[8] = local_ovule[2] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[2];
	out_immigrant_dip_2[9] = local_ovule[3] * immigrant_pollen[3];

	double out_im_SUM_2(0.0);
	for (int i(0); i < 10; ++i) { out_im_SUM_2 += out_immigrant_dip_2[i]; }
	double out_im_FREQ_2[10];
	for (int i(0); i < 10; ++i) { out_im_FREQ_2[i] = out_immigrant_dip_2[i] / out_im_SUM_2; }

	// Total seed Genotypes
	for (int i(0); i < 10; ++i) {
		dip_FREQ_2[i] =((self_r_2 * self_dip_FREQ[i]) + ((1 - self_r_2) * (1 - m_h_2) * out_local_FREQ_2[i]) + 
		((1 - self_r_2) * m_h_2 * out_im_FREQ_2[i]));
	}
}
void REPRODUCTION_POP1_STO(const double self_r_1, const unsigned int dip_IND_1[10], const double Me_Mu_Matrix_1[][4], double* dip_FREQ_1,const double m_h_1, const unsigned int dip_IND_2[10], const double Me_Mu_Matrix_2[][4])
{	// Seed genotypes produced through selfing
	double self_dip[10] = {};
	for (int i(0); i < 10; ++i) {
		self_dip[0] += dip_IND_1[i] * Me_Mu_Matrix_1[i][0] * Me_Mu_Matrix_1[i][0];
		self_dip[1] += dip_IND_1[i] * Me_Mu_Matrix_1[i][0] * Me_Mu_Matrix_1[i][1] * 2;
		self_dip[2] += dip_IND_1[i] * Me_Mu_Matrix_1[i][0] * Me_Mu_Matrix_1[i][2] * 2;
		self_dip[3] += dip_IND_1[i] * Me_Mu_Matrix_1[i][0] * Me_Mu_Matrix_1[i][3] * 2;
		self_dip[4] += dip_IND_1[i] * Me_Mu_Matrix_1[i][1] * Me_Mu_Matrix_1[i][1];
		self_dip[5] += dip_IND_1[i] * Me_Mu_Matrix_1[i][1] * Me_Mu_Matrix_1[i][2] * 2;
		self_dip[6] += dip_IND_1[i] * Me_Mu_Matrix_1[i][1] * Me_Mu_Matrix_1[i][3] * 2;
		self_dip[7] += dip_IND_1[i] * Me_Mu_Matrix_1[i][2] * Me_Mu_Matrix_1[i][2];
		self_dip[8] += dip_IND_1[i] * Me_Mu_Matrix_1[i][2] * Me_Mu_Matrix_1[i][3] * 2;
		self_dip[9] += dip_IND_1[i] * Me_Mu_Matrix_1[i][3] * Me_Mu_Matrix_1[i][3];
	}

	double self_dip_SUM(0.0);
	for (int i(0); i < 10; ++i) { self_dip_SUM += self_dip[i]; }
	double self_dip_FREQ[10];
	for (int i(0); i < 10; ++i) { self_dip_FREQ[i] = self_dip[i] / self_dip_SUM; }

	// Seed genotypes produced through outcrossing local
	double local_ovule[4] = {};
	double local_pollen[4] = {};
	GAMETE_PROD_STO(dip_IND_1, Me_Mu_Matrix_1, local_ovule);
	GAMETE_PROD_STO(dip_IND_1, Me_Mu_Matrix_1, local_pollen);

	double out_local_dip_1[10] = {};
	out_local_dip_1[0] = local_ovule[0] * local_pollen[0];
	out_local_dip_1[1] = local_ovule[0] * local_pollen[1] + local_ovule[1] * local_pollen[0]; 
	out_local_dip_1[2] = local_ovule[0] * local_pollen[2] + local_ovule[2] * local_pollen[0];
	out_local_dip_1[3] = local_ovule[0] * local_pollen[3] + local_ovule[3] * local_pollen[0];
	out_local_dip_1[4] = local_ovule[1] * local_pollen[1];
	out_local_dip_1[5] = local_ovule[1] * local_pollen[2] + local_ovule[2] * local_pollen[1];
	out_local_dip_1[6] = local_ovule[1] * local_pollen[3] + local_ovule[3] * local_pollen[1];
	out_local_dip_1[7] = local_ovule[2] * local_pollen[2];
	out_local_dip_1[8] = local_ovule[2] * local_pollen[3] + local_ovule[3] * local_pollen[2];
	out_local_dip_1[9] = local_ovule[3] * local_pollen[3];

	double out_local_SUM_1(0.0);
	for (int i(0); i < 10; ++i) { out_local_SUM_1 += out_local_dip_1[i]; }
	double out_local_FREQ_1[10];
	for (int i(0); i < 10; ++i) { out_local_FREQ_1[i] = out_local_dip_1[i] / out_local_SUM_1; }

	// Seed genotypes produced through outcrossing immigrants
	double immigrant_pollen[4] = {};
	GAMETE_PROD_STO(dip_IND_2, Me_Mu_Matrix_2,immigrant_pollen);
	
	double out_immigrant_dip_1[10] = {};
	out_immigrant_dip_1[0] = local_ovule[0] * immigrant_pollen[0];
	out_immigrant_dip_1[1] = local_ovule[0] * immigrant_pollen[1] + local_ovule[1] * immigrant_pollen[0] ;
	out_immigrant_dip_1[2] = local_ovule[0] * immigrant_pollen[2] + local_ovule[2] * immigrant_pollen[0];
	out_immigrant_dip_1[3] = local_ovule[0] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[0];
	out_immigrant_dip_1[4] = local_ovule[1] * immigrant_pollen[1];
	out_immigrant_dip_1[5] = local_ovule[1] * immigrant_pollen[2] + local_ovule[2] * immigrant_pollen[1];
	out_immigrant_dip_1[6] = local_ovule[1] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[1];
	out_immigrant_dip_1[7] = local_ovule[2] * immigrant_pollen[2];
	out_immigrant_dip_1[8] = local_ovule[2] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[2];
	out_immigrant_dip_1[9] = local_ovule[3] * immigrant_pollen[3];

	double out_im_SUM_1(0.0);
	for (int i(0); i < 10; ++i) { out_im_SUM_1 += out_immigrant_dip_1[i]; }
	double out_im_FREQ_1[10];
	for (int i(0); i < 10; ++i) { out_im_FREQ_1[i] = out_immigrant_dip_1[i] / out_im_SUM_1; }

	// Total seed Genotypes
	for (int i(0); i < 10; ++i) {
		dip_FREQ_1[i] =((self_r_1 * self_dip_FREQ[i]) + ((1 - self_r_1) * (1 - m_h_1) * out_local_FREQ_1[i]) + 
		((1 - self_r_1) * m_h_1 * out_im_FREQ_1[i]));
	}
}
void REPRODUCTION_POP2_STO(const double self_r_2, const unsigned int dip_IND_2[10], const double Me_Mu_Matrix_2[][4], double* dip_FREQ_2,const double m_h_2, const unsigned int dip_IND_1[10], const double Me_Mu_Matrix_1[][4])
{	// Seed genotypes produced through selfing
	double self_dip[10] = {};
	for (int i(0); i < 10; ++i) {
		self_dip[0] += dip_IND_2[i] * Me_Mu_Matrix_2[i][0] * Me_Mu_Matrix_2[i][0];
		self_dip[1] += dip_IND_2[i] * Me_Mu_Matrix_2[i][0] * Me_Mu_Matrix_2[i][1] * 2;
		self_dip[2] += dip_IND_2[i] * Me_Mu_Matrix_2[i][0] * Me_Mu_Matrix_2[i][2] * 2;
		self_dip[3] += dip_IND_2[i] * Me_Mu_Matrix_2[i][0] * Me_Mu_Matrix_2[i][3] * 2;
		self_dip[4] += dip_IND_2[i] * Me_Mu_Matrix_2[i][1] * Me_Mu_Matrix_2[i][1];
		self_dip[5] += dip_IND_2[i] * Me_Mu_Matrix_2[i][1] * Me_Mu_Matrix_2[i][2] * 2;
		self_dip[6] += dip_IND_2[i] * Me_Mu_Matrix_2[i][1] * Me_Mu_Matrix_2[i][3] * 2;
		self_dip[7] += dip_IND_2[i] * Me_Mu_Matrix_2[i][2] * Me_Mu_Matrix_2[i][2];
		self_dip[8] += dip_IND_2[i] * Me_Mu_Matrix_2[i][2] * Me_Mu_Matrix_2[i][3] * 2;
		self_dip[9] += dip_IND_2[i] * Me_Mu_Matrix_2[i][3] * Me_Mu_Matrix_2[i][3];
	}

	double self_dip_SUM(0.0);
	for (int i(0); i < 10; ++i) { self_dip_SUM += self_dip[i]; }
	double self_dip_FREQ[10];
	for (int i(0); i < 10; ++i) { self_dip_FREQ[i] = self_dip[i] / self_dip_SUM; }

	// Seed genotypes produced through outcrossing local
	double local_ovule[4] = {};
	double local_pollen[4] = {};
	GAMETE_PROD_STO(dip_IND_2, Me_Mu_Matrix_2, local_ovule);
	GAMETE_PROD_STO(dip_IND_2, Me_Mu_Matrix_2, local_pollen);

	double out_local_dip_2[10] = {};
	out_local_dip_2[0] = local_ovule[0] * local_pollen[0];
	out_local_dip_2[1] = local_ovule[0] * local_pollen[1] + local_ovule[1] * local_pollen[0]; 
	out_local_dip_2[2] = local_ovule[0] * local_pollen[2] + local_ovule[2] * local_pollen[0];
	out_local_dip_2[3] = local_ovule[0] * local_pollen[3] + local_ovule[3] * local_pollen[0];
	out_local_dip_2[4] = local_ovule[1] * local_pollen[1];
	out_local_dip_2[5] = local_ovule[1] * local_pollen[2] + local_ovule[2] * local_pollen[1];
	out_local_dip_2[6] = local_ovule[1] * local_pollen[3] + local_ovule[3] * local_pollen[1];
	out_local_dip_2[7] = local_ovule[2] * local_pollen[2];
	out_local_dip_2[8] = local_ovule[2] * local_pollen[3] + local_ovule[3] * local_pollen[2];
	out_local_dip_2[9] = local_ovule[3] * local_pollen[3];

	double out_local_SUM_2(0.0);
	for (int i(0); i < 10; ++i) { out_local_SUM_2 += out_local_dip_2[i]; }
	double out_local_FREQ_2[10];
	for (int i(0); i < 10; ++i) { out_local_FREQ_2[i] = out_local_dip_2[i] / out_local_SUM_2; }

	// Seed genotypes produced through outcrossing immigrants
	double immigrant_pollen[4] = {};
	GAMETE_PROD_STO(dip_IND_1, Me_Mu_Matrix_1,immigrant_pollen);
	
	double out_immigrant_dip_2[10] = {};
	out_immigrant_dip_2[0] = local_ovule[0] * immigrant_pollen[0];
	out_immigrant_dip_2[1] = local_ovule[0] * immigrant_pollen[1] + local_ovule[1] * immigrant_pollen[0] ;
	out_immigrant_dip_2[2] = local_ovule[0] * immigrant_pollen[2] + local_ovule[2] * immigrant_pollen[0];
	out_immigrant_dip_2[3] = local_ovule[0] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[0];
	out_immigrant_dip_2[4] = local_ovule[1] * immigrant_pollen[1];
	out_immigrant_dip_2[5] = local_ovule[1] * immigrant_pollen[2] + local_ovule[2] * immigrant_pollen[1];
	out_immigrant_dip_2[6] = local_ovule[1] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[1];
	out_immigrant_dip_2[7] = local_ovule[2] * immigrant_pollen[2];
	out_immigrant_dip_2[8] = local_ovule[2] * immigrant_pollen[3] + local_ovule[3] * immigrant_pollen[2];
	out_immigrant_dip_2[9] = local_ovule[3] * immigrant_pollen[3];

	double out_im_SUM_2(0.0);
	for (int i(0); i < 10; ++i) { out_im_SUM_2 += out_immigrant_dip_2[i]; }
	double out_im_FREQ_2[10];
	for (int i(0); i < 10; ++i) { out_im_FREQ_2[i] = out_immigrant_dip_2[i] / out_im_SUM_2; }

	// Total seed Genotypes
	for (int i(0); i < 10; ++i) {
		dip_FREQ_2[i] =((self_r_2 * self_dip_FREQ[i]) + ((1 - self_r_2) * (1 - m_h_2) * out_local_FREQ_2[i]) + 
		((1 - self_r_2) * m_h_2 * out_im_FREQ_2[i]));
	}
}
// Compute genotypic frequencies after selection, given genotypic frequencies before selection, and the fitness landscape
void SELECTION(const double dip_FREQ_pre_sel[10], const double Fitness[10], double* dip_FREQ_post_sel)
{
	double dip_post_sel[10] = {};
	dip_post_sel[0] = Fitness[0] * dip_FREQ_pre_sel[0];
	dip_post_sel[1] = Fitness[1] * dip_FREQ_pre_sel[1];
	dip_post_sel[2] = Fitness[2] * dip_FREQ_pre_sel[2];
	dip_post_sel[3] = Fitness[3] * dip_FREQ_pre_sel[3]; 
	dip_post_sel[4] = Fitness[4] * dip_FREQ_pre_sel[4];
	dip_post_sel[5] = Fitness[5] * dip_FREQ_pre_sel[5];
	dip_post_sel[6] = Fitness[6] * dip_FREQ_pre_sel[6];
	dip_post_sel[7] = Fitness[7] * dip_FREQ_pre_sel[7];
	dip_post_sel[8] = Fitness[8] * dip_FREQ_pre_sel[8];
	dip_post_sel[9] = Fitness[9] * dip_FREQ_pre_sel[9];

	double SUM(0.0);
	for (int i(0); i < 10; ++i) { SUM += dip_post_sel[i]; }
	for (int i(0); i < 10; ++i) { dip_FREQ_post_sel[i] = dip_post_sel[i] / SUM; }
}

// Diploid migration 
void SEED_MIGRATION_POP1(const double m_d_2, const double selected_1[10], const double selected_2[10], double* final_seeds_1 )
{
	for (int j(0); j < 10; ++j){
		final_seeds_1[j] = (1 - m_d_2) * selected_1[j] + m_d_2 * selected_2[j]; //Selected seeds from island migrate to continent
	}
}
void SEED_MIGRATION_POP2(const double m_d_1, const double selected_1[10], const double selected_2[10], double* final_seeds_2 )
{
	for (int j(0); j < 10; ++j){
		final_seeds_2[j] = (1 - m_d_1) * selected_2[j] + m_d_1 * selected_1[j]; //Selected seeds from continent migrate to island
	}
}

// Drift
void DRIFT(const gsl_rng* r, const unsigned int N, const double dip_FREQ[], unsigned int* dip_IND)
{
	gsl_ran_multinomial(r, 10, N, dip_FREQ, dip_IND);
}
void DRIFT_WITH_DIRICHLET(const gsl_rng* r, const double Dirichlet_Cst, const unsigned int N, double* dip_FREQ, unsigned int* dip_IND)
{
	for (int i(0); i < 10; ++i) { /*mutiply seed frequencies by the Dirichlet multiplier*/
		dip_FREQ[i] = Dirichlet_Cst * dip_FREQ[i];
	}

	double adults_FREQ_BIS[10] = {}; /*create an empty array for the Dirichlet output*/

	gsl_ran_dirichlet(r, 10, dip_FREQ, adults_FREQ_BIS); /*Dirichlet sampling*/

	gsl_ran_multinomial(r, 10, N, adults_FREQ_BIS, dip_IND);
}
// Define Dirichlet functions given N
void DIRICHLET_CURVED(const unsigned int N, double Dirichlet[11])
{
	Dirichlet[0] = 0.0769231 * ((double)-1250 + (double)1237 * N);
	Dirichlet[1] = 0.037037 * ((double)-10000 + (double)9973 * N);
	Dirichlet[2] = 1e5 * N;
	Dirichlet[3] = 0.333333 * ((double)-1250 + (double)1247 * N);
	Dirichlet[4] = 0.00917431 * ((double)-10000 + (double)9891 * N);
	Dirichlet[5] = 0.00355872 * ((double)-10000 + (double)9719 * N);
	Dirichlet[6] = 0.00341297 * ((double)-5000 + (double)4707 * N);
	Dirichlet[7] = 0.000887311 * ((double)-10000 + (double)8873 * N);
	Dirichlet[8] = 0.000465766 * ((double)-10000 + (double)7853 * N);
	Dirichlet[9] = 0.000232937 * ((double)-10000 + (double)5707 * N);
	Dirichlet[10] = 0.000115701 * ((double)-10000 + (double)1357 * N);
}
// Compute allele frequencies
void ALLELE_FREQ_COMP(const double* dip_FREQ, double* allele_FREQ) 
{
	allele_FREQ[0] = ((2 * (double)dip_FREQ[0]) + (2 * (double)dip_FREQ[1]) + (double)dip_FREQ[2] + (double)dip_FREQ[3] + (2 * (double)dip_FREQ[4]) + (double)dip_FREQ[5] + (double)dip_FREQ[6]) / 2.0; // A
	allele_FREQ[1] = ((double)dip_FREQ[2] + (double)dip_FREQ[3] + (double)dip_FREQ[5] + (double)dip_FREQ[6] + (2 * (double)dip_FREQ[7]) + (2 * (double)dip_FREQ[8]) + (2 * (double)dip_FREQ[9])) / 2.0; // a
	allele_FREQ[2] = ((2 * (double)dip_FREQ[0]) + (double)dip_FREQ[1] + (2 * (double)dip_FREQ[2]) + (double)dip_FREQ[3] + (double)dip_FREQ[5] + (2 * (double)dip_FREQ[7]) + (double)dip_FREQ[8]) / 2.0; // B
	allele_FREQ[3] = ((double)dip_FREQ[1] + (double)dip_FREQ[3] + (2 * (double)dip_FREQ[4]) + (double)dip_FREQ[5] + (2 * (double)dip_FREQ[6]) + (double)dip_FREQ[8] + (2 * (double)dip_FREQ[9])) / 2.0; // b
}

void ALLELE_FREQ_COMP_STO(const unsigned int* dip_IND, double N, double* allele_FREQ) 
{
	allele_FREQ[0] = ((2 * (double)dip_IND[0]) + (2 * (double)dip_IND[1]) + (double)dip_IND[2] + (double)dip_IND[3] + (2 * (double)dip_IND[4]) + (double)dip_IND[5] + (double)dip_IND[6]) / (2 * N); // A
	allele_FREQ[1] = ((double)dip_IND[2] + (double)dip_IND[3] + (double)dip_IND[5] + (double)dip_IND[6] + (2 * (double)dip_IND[7]) + (2 * (double)dip_IND[8]) + (2 * (double)dip_IND[9])) / (2 * N); // a
	allele_FREQ[2] = ((2 * (double)dip_IND[0]) + (double)dip_IND[1] + (2 * (double)dip_IND[2]) + (double)dip_IND[3] + (double)dip_IND[5] + (2 * (double)dip_IND[7]) + (double)dip_IND[8]) / (2 * N); // B
	allele_FREQ[3] = ((double)dip_IND[1] + (double)dip_IND[3] + (2 * (double)dip_IND[4]) + (double)dip_IND[5] + (2 * (double)dip_IND[6]) + (double)dip_IND[8] + (2 * (double)dip_IND[9])) / (2 * N); // b
}
