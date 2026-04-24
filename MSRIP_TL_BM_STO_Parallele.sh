#!/bin/bash

# Script local pour PC - adaptation de MSRIP_TL_BM_Parallel.sh
# Script pour Secondary Contact case

# Nettoyage des fichiers précédents
rm -f ./msri.exe
rm -f ./ColumnHeader*.csv
rm -f ./Output*.csv

# Compilation 
echo "Compilation loading..."
g++ -O3 -Wall -Wextra -std=c++11 -I/mingw64/include -L/mingw64/lib -o msri.exe ./MSRIP_TL_BM_STO_Main.cpp -lgsl -lgslcblas -lm

if [ $? -ne 0 ]; then
    echo "ERROR : compilation failed."
    exit 1
fi

# Paramètres fixes
threshold=9999
iteration=10
span=100
interval=10

touch Output_TL_BM_SC_STO.csv

echo "Simulations starting..."

for N1 in 1000; do
  for N2 in 1000; do
    for self_r_1 in 0; do
      for self_r_2 in 0; do
        for mu_1 in 0; do
          for mu_2 in 0; do
            for ha_1 in 0.5; do
              for sa_1 in 0; do
                for hb_1 in 0.5; do
                  for sb_1 in 0; do
                    for epsilon_1_1 in 0; do        # epsilon_1_1 = epsilon_2_1
                      for epsilon_3_1 in 0; do       # epsilon_3_1
                        for epsilon_4_1 in 0; do     # epsilon_4_1
                          for ha_2 in 0.5; do
                            for sa_2 in 0.1; do
                              for hb_2 in 0.5; do
                                for sb_2 in -0.1; do
                                  for epsilon_1_2 in -0.1; do      # epsilon_1_2 = epsilon_2_2
                                    for epsilon_3_2 in -0.2; do    # epsilon_3_2
                                      for epsilon_4_2 in -0.4; do  # epsilon_4_2
                                        for rec_1 in 0; do
                                          for rec_2 in 0.001 0.01 0.1 0.5; do
                                            for m_h_1 in 0; do    # flux du continent vers l'île
                                              for m_h_2 in 0; do  # flux de l'île vers le continent
                                                for m_d_1 in $(seq 0.01 0.001 0.2); do
                                                  for m_d_2 in 0; do
                                                    echo "Running m_d_1=${m_d_1}, self_r_2=${self_r_2}..."
                                                    ./msri.exe \
                                                      ${threshold} ${iteration} ${span} ${interval} \
                                                      ${N1} ${N2} \
                                                      ${self_r_1} ${self_r_2} \
                                                      ${mu_1} ${mu_2} \
                                                      ${ha_1} ${sa_1} ${hb_1} ${sb_1} ${epsilon_1_1} ${epsilon_3_1} ${epsilon_4_1} \
                                                      ${ha_2} ${sa_2} ${hb_2} ${sb_2} ${epsilon_1_2} ${epsilon_3_2} ${epsilon_4_2} \
                                                      ${rec_1} ${rec_2} \
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

# CSV Header
printf "threshold,iteration,span,interval,N1,N2,self_r_1,self_r_2," > ColumnHeader_TL_BM_SC_STO.csv
printf "ha_2,sa_2,hb_2,sb_2,epsilon_1_2,epsilon_3_2,epsilon_4_2,rec_2," >> ColumnHeader_TL_BM_SC_STO.csv
printf "m_h_1,m_h_2,m_d_1,m_d_2,gen," >> ColumnHeader_TL_BM_SC_STO.csv
printf "A_FREQ_1,a_FREQ_1,B_FREQ_1,b_FREQ_1," >> ColumnHeader_TL_BM_SC_STO.csv
printf "A_FREQ_2,a_FREQ_2,B_FREQ_2,b_FREQ_2," >> ColumnHeader_TL_BM_SC_STO.csv
printf "ABAB_0_pop1,ABAB_0_pop2,ABAb_0_pop1,ABAb_0_pop2,ABaB_0_pop1,ABaB_0_pop2,ABab_0_pop1,ABab_0_pop2,AbAb_0_pop1,AbAb_0_pop2,AbaB_0_pop1,AbaB_0_pop2,Abab_0_pop1,Abab_0_pop2,aBaB_0_pop1,aBaB_0_pop2,aBab_0_pop1,aBab_0_pop2,abab_0_pop1,abab_0_pop2" >> ColumnHeader_TL_BM_SC_STO.csv

counter=1
while [ $counter -lt $span ]
do
  printf ",ABAB_"$((interval*counter))"_pop1,ABAB_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",ABAb_"$((interval*counter))"_pop1,ABAb_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",ABaB_"$((interval*counter))"_pop1,ABaB_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",ABab_"$((interval*counter))"_pop1,ABab_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",AbAb_"$((interval*counter))"_pop1,AbAb_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",AbaB_"$((interval*counter))"_pop1,AbaB_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",Abab_"$((interval*counter))"_pop1,Abab_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",aBaB_"$((interval*counter))"_pop1,aBaB_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",aBab_"$((interval*counter))"_pop1,aBab_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  printf ",abab_"$((interval*counter))"_pop1,abab_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC_STO.csv
  counter=$(($counter+1))
done

printf "\n" >> ColumnHeader_TL_BM_SC_STO.csv

echo "ColumnHeader_TL_BM_SC_STO.csv created."
