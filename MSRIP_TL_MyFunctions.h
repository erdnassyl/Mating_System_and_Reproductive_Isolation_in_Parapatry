// Mating System and Reproductive Isolation in Allopatry
// Two-Loci C++ Models � FUNCTIONS
// Lucas Marie-Orleach

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
void FITNESS_LANDSCAPE_BM(const double sa, const double ha, const double sb, const double hb, const double s_B, const double h_B, const double k_B, double* Fitness)
{
	Fitness[0] = 1;
	Fitness[1] = (1 + (hb * sb));
	Fitness[2] = (1 + (ha * sa));
	Fitness[3] = (1 + (ha * sa)) * (1 + (hb * sb)) * (1 + (k_B * h_B * s_B));
	Fitness[4] = (1 + sb);
	Fitness[5] = (1 + (ha * sa)) * (1 + (hb * sb)) * (1 + (k_B * h_B * s_B));
	Fitness[6] = (1 + (ha * sa)) * (1 + sb) * (1 + (h_B * s_B));
	Fitness[7] = (1 + sa);
	Fitness[8] = (1 + sa) * (1 + (hb * sb)) * (1 + (h_B * s_B));
	Fitness[9] = (1 + sa) * (1 + sb) * (1 + s_B);
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
// Define Reproduction without pollen migration function given selfing rate, genotypes of adult individuals, the Meiosis_Mutation matrix, and the Seeds's genotypic frequency 
void REPRODUCTION(const double self_r_1, const double dip_IND_1[10], const double Me_Mu_Matrix_1[][4], double* dip_FREQ_1,const double m_h_1, const double dip_IND_2[10], const double Me_Mu_Matrix_2[][4])
	
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
void SEED_MIGRATION(const double m_d_1,const double m_d_2, const double Fitness_1[10], const double Fitness_2[10], const double pre_sel_1[10],const double pre_sel_2[10], double* final_seeds_1, double* final_seeds_2 )
{
	double selected_local[10]={0};
	double selected_immigrant[10]={0};
	SELECTION(pre_sel_1,Fitness_1,selected_local);
	SELECTION(pre_sel_2,Fitness_2,selected_immigrant);

	for (int j(0); j < 10; ++j){
		final_seeds_1[j] = (1 - m_d_1) * selected_local[j] + m_d_1 * selected_immigrant[j];
		final_seeds_2[j] = (1 - m_d_2) * selected_immigrant[j] + m_d_2 * selected_local[j];
	}
}
// Compute allele frequencies
void ALLELE_FREQ_COMP(const double* dip_FREQ, double* allele_FREQ) 
{
	allele_FREQ[0] = ((2 * (double)dip_FREQ[0]) + (2 * (double)dip_FREQ[1]) + (double)dip_FREQ[2] + (double)dip_FREQ[3] + (2 * (double)dip_FREQ[4]) + (double)dip_FREQ[5] + (double)dip_FREQ[6]) / 2.0; // A
	allele_FREQ[1] = ((double)dip_FREQ[2] + (double)dip_FREQ[3] + (double)dip_FREQ[5] + (double)dip_FREQ[6] + (2 * (double)dip_FREQ[7]) + (2 * (double)dip_FREQ[8]) + (2 * (double)dip_FREQ[9])) / 2.0; // a
	allele_FREQ[2] = ((2 * (double)dip_FREQ[0]) + (double)dip_FREQ[1] + (2 * (double)dip_FREQ[2]) + (double)dip_FREQ[3] + (double)dip_FREQ[5] + (2 * (double)dip_FREQ[7]) + (double)dip_FREQ[8]) / 2.0; // B
	allele_FREQ[3] = ((double)dip_FREQ[1] + (double)dip_FREQ[3] + (2 * (double)dip_FREQ[4]) + (double)dip_FREQ[5] + (2 * (double)dip_FREQ[6]) + (double)dip_FREQ[8] + (2 * (double)dip_FREQ[9])) / 2.0; // b
}
