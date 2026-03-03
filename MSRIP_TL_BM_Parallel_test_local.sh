#!/bin/bash

# ================================================================
# Script local pour PC - adaptation de MSRIP_TL_BM_Parallel.sh
# =================================================================

# Nettoyage des fichiers précédents
rm -f ./msri.exe
rm -f ./ColumnHeader*.csv
rm -f ./Output*.csv

# --- Compilation ---
echo "Compilation en cours..."
g++ -O3 -Wall -Wextra -std=c++11 -I/mingw64/include -L/mingw64/lib -o msri.exe ./MSRIP_TL_BM_Main.cpp -lgsl -lgslcblas -lm

if [ $? -ne 0 ]; then
    echo "ERREUR : la compilation a échoué, arrêt du script."
    exit 1
fi
echo "Compilation réussie !"

# Paramètres fixes
threshold=10000
iteration=1
span=0
interval=10

touch Output_TL_BM.csv

# Boucle de simulations 
echo "Lancement des simulations..."

for N_1 in 100
do
 for N_2 in 100
 do
  for self_r_1 in 0 
  do
   for self_r_2 in 0 
   do
    for mu_1 in 2.5e-4
    do
     for mu_2 in 2.5e-4
     do
      for ha_1 in 0.5
      do
       for sa_1 in 0.1
       do
        for hb_1 in 0.5
        do
         for sb_1 in 0
         do
          for ha_2 in 0.5
          do
           for sa_2 in 0.1
           do
            for hb_2 in 0.5
            do
             for sb_2 in 0
             do
              for rec_1 in 0.5
              do
               for rec_2 in 0.5
               do
                for h_B_1 in 0.5
                do
                 for s_B_1 in 0.1
                 do
                  for s_B_exp_1 in -4
                  do
                   for h_B_2 in 0.5
                   do
                    for s_B_2 in 0.1
                    do
                     for s_B_exp_2 in -4
                     do
                      for m_h_1 in 0.4
                      do
                       for m_h_2 in 0
                       do
                        for m_d_1 in 0
                        do
                         for m_d_2 in 0
                         do
                          echo "Running..."
                          ./msri.exe ${threshold} ${iteration} ${span} ${interval} \
                                     ${N_1} ${N_2} \
                                     ${self_r_1} ${self_r_2} \
                                     ${mu_1} ${mu_2} \
                                     ${ha_1} ${sa_1} ${hb_1} ${sb_1} \
                                     ${ha_2} ${sa_2} ${hb_2} ${sb_2} \
                                     ${rec_1} ${rec_2} \
                                     ${h_B_1} ${s_B_1} ${s_B_exp_1} \
                                     ${h_B_2} ${s_B_2} ${s_B_exp_2} \
                                     ${m_h_1} ${m_h_2} ${m_d_1} ${m_d_2}
                         done
                        done
                       done
                      done
                     done
                    done
                   done
                  done
                 done
                done
               done
              done
             done
            done
           done
          done
         done
        done
       done
      done
     done
    done
   done
  done
 done
done

echo "Simulations terminées !"

# En-tête CSV 
printf "threshold,iteration,span,interval,N_1,N_2,self_r_1,self_r_2,mu_Aa_1,mu_aA_1,mu_Bb_1,mu_bB_1,mu_Aa_2,mu_aA_2,mu_Bb_2,mu_bB_2,ha_1,sa_1,hb_1,sb_1,ha_2,sa_2,hb_2,sb_2,rec_1,rec_2,h_B_1,k_B_1,s_B_1,h_B_2,k_B_2,s_B_2,m_h_1,m_h_2,m_d_1,m_d_2,gen," >> ColumnHeader_BM.csv
printf "A_FREQ_1,a_FREQ_1,B_FREQ_1,b_FREQ_1,A_FREQ_2,a_FREQ_2,B_FREQ_2,b_FREQ_2," >> ColumnHeader_BM.csv
printf "AB/AB_0,AB/Ab_0,AB/aB_0,AB/ab_0,Ab/Ab_0,Ab/aB_0,Ab/ab_0,aB/aB_0,aB/ab_0,ab/ab_0_pop1," >> ColumnHeader_BM.csv
printf "AB/AB_0,AB/Ab_0,AB/aB_0,AB/ab_0,Ab/Ab_0,Ab/aB_0,Ab/ab_0,aB/aB_0,aB/ab_0,ab/ab_0_pop2" >> ColumnHeader_BM.csv

counter=1
while [ $counter -le $span ]
do
  printf ",AB/AB_-"$((interval*counter))"_pop1,AB/Ab_-"$((interval*counter))"_pop1,AB/aB_-"$((interval*counter))"_pop1,AB/ab_-"$((interval*counter))"_pop1,Ab/Ab_-"$((interval*counter))"_pop1,Ab/aB_-"$((interval*counter))"_pop1,Ab/ab_-"$((interval*counter))"_pop1,aB/aB_-"$((interval*counter))"_pop1,aB/ab_-"$((interval*counter))"_pop1,ab/ab_-"$((interval*counter))"_pop1" >> ColumnHeader_BM.csv
  printf ",AB/AB_-"$((interval*counter))"_pop2,AB/Ab_-"$((interval*counter))"_pop2,AB/aB_-"$((interval*counter))"_pop2,AB/ab_-"$((interval*counter))"_pop2,Ab/Ab_-"$((interval*counter))"_pop2,Ab/aB_-"$((interval*counter))"_pop2,Ab/ab_-"$((interval*counter))"_pop2,aB/aB_-"$((interval*counter))"_pop2,aB/ab_-"$((interval*counter))"_pop2,ab/ab_-"$((interval*counter))"_pop2" >> ColumnHeader_BM.csv
  counter=$(($counter+1))
done

printf "\n" >> ColumnHeader_BM.csv

echo "Fichier ColumnHeader_BM.csv créé."