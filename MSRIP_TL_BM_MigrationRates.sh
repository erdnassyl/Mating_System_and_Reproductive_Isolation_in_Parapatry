#!/bin/bash
# Script local pour PC - adaptation de MSRIP_TL_BM_Parallel.sh
# Script pour Secondary Contact case

# Nettoyage des fichiers précédents
rm -f ./msri.exe
rm -f ./ColumnHeader*.csv
rm -f ./Output*.csv

# Compilation 
echo "Compilation loading..."
g++ -O3 -Wall -Wextra -std=c++11 -I/mingw64/include -L/mingw64/lib -o msri.exe ./MSRIP_TL_BM_SC_Main.cpp -lgsl -lgslcblas -lm
if [ $? -ne 0 ]; then
    echo "ERROR : compilation failed."
    exit 1
fi

# Paramètres fixes
threshold=1000000
iteration=1
span=0
interval=10

# Initialise des variables mémoire
prev_a_freq_2=1  
prev_m_d_1=0

# Fichier de sortie
touch Output_TL_BM_SC.csv
SAUT_FILE="migration_threshold.csv"
echo "self_r_2,rec_2,sa,sb,m_d_1_critique,a_freq_2" > "$SAUT_FILE"

echo "Simulations starting..."
for self_r_1 in 0; do
  for self_r_2 in $(seq 0.8 0.1 0.9); do
    for mu_1 in 0; do
      for mu_2 in 0; do
        for ha_1 in 0.5; do
          for sa_1 in 0; do
            for hb_1 in 0.5; do
              for sb_1 in 0; do
                for epsilon_1_1 in 0; do
                  for epsilon_3_1 in 0; do
                    for epsilon_4_1 in 0; do
                      for ha_2 in 0.5; do
                        for sa_2 in 0.1; do
                          for hb_2 in 0.5; do
                            for sb_2 in -0.1; do
                              for epsilon_1_2 in 0; do
                                for epsilon_3_2 in -0.2; do
                                  for epsilon_4_2 in -0.4; do
                                    for rec_1 in 0; do
                                      for rec_2 in 0.001 0.01 0.1 0.5; do

                                        # Réinitialisation pour chaque nouvelle combinaison
                                        prev_a_freq_2=1
                                        prev_m_d_1=0

                                        for m_h_1 in 0; do
                                          for m_h_2 in 0; do
                                            for m_d_1 in $(seq 0.02 0.0001 0.4); do
                                              for m_d_2 in 0; do

                                                echo "Running m_d_1=$m_d_1 rec_2=$rec_2 self_r_2=$self_r_2..."

                                                # Vide le CSV avant chaque run
                                                > Output_TL_BM_SC.csv

                                                ./msri.exe \
                                                  ${threshold} ${iteration} ${span} ${interval} \
                                                  ${self_r_1} ${self_r_2} \
                                                  ${mu_1} ${mu_2} \
                                                  ${ha_1} ${sa_1} ${hb_1} ${sb_1} ${epsilon_1_1} ${epsilon_3_1} ${epsilon_4_1} \
                                                  ${ha_2} ${sa_2} ${hb_2} ${sb_2} ${epsilon_1_2} ${epsilon_3_2} ${epsilon_4_2} \
                                                  ${rec_1} ${rec_2} \
                                                  ${m_h_1} ${m_h_2} ${m_d_1} ${m_d_2}

                                                a_freq_2=$(tail -n 1 Output_TL_BM_SC.csv | awk -F',' '{print $25}')

                                                # Détection du saut : prev >= 0.001 ET actuel < 0.001
                                                is_jump=$(awk -v actuel="$a_freq_2" -v prev="$prev_a_freq_2" 'BEGIN {
                                                  if (prev >= 0.001 && actuel < 0.001) print 1; else print 0
                                                }')

                                                if [ "$is_jump" -eq 1 ]; then
                                                  echo "${self_r_2},${rec_2},${sa_2},${sb_2},${prev_m_d_1},${prev_a_freq_2}" >> "$SAUT_FILE"
                                                  echo "Saut détecté à m=$m_d_1 ! Enregistrement de m=$prev_m_d_1"
                                                  break 3
                                                fi

                                                prev_a_freq_2=$a_freq_2
                                                prev_m_d_1=$m_d_1

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
printf "threshold,iteration,span,interval,self_r_1,self_r_2," > ColumnHeader_TL_BM_SC.csv
printf "ha_2,sa_2,hb_2,sb_2,epsilon_1_2,epsilon_3_2,epsilon_4_2,rec_2," >> ColumnHeader_TL_BM_SC.csv
printf "m_h_1,m_h_2,m_d_1,m_d_2,gen," >> ColumnHeader_TL_BM_SC.csv
printf "A_FREQ_1,a_FREQ_1,B_FREQ_1,b_FREQ_1," >> ColumnHeader_TL_BM_SC.csv
printf "A_FREQ_2,a_FREQ_2,B_FREQ_2,b_FREQ_2," >> ColumnHeader_TL_BM_SC.csv
printf "ABAB_0_pop1,ABAB_0_pop2,ABAb_0_pop1,ABAb_0_pop2,ABaB_0_pop1,ABaB_0_pop2,ABab_0_pop1,ABab_0_pop2,AbAb_0_pop1,AbAb_0_pop2,AbaB_0_pop1,AbaB_0_pop2,Abab_0_pop1,Abab_0_pop2,aBaB_0_pop1,aBaB_0_pop2,aBab_0_pop1,aBab_0_pop2,abab_0_pop1,abab_0_pop2" >> ColumnHeader_TL_BM_SC.csv

counter=0
while [ $counter -lt $span ]; do
  printf ",ABAB_"$((interval*counter))"_pop1,ABAB_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",ABAb_"$((interval*counter))"_pop1,ABAb_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",ABaB_"$((interval*counter))"_pop1,ABaB_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",ABab_"$((interval*counter))"_pop1,ABab_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",AbAb_"$((interval*counter))"_pop1,AbAb_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",AbaB_"$((interval*counter))"_pop1,AbaB_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",Abab_"$((interval*counter))"_pop1,Abab_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",aBaB_"$((interval*counter))"_pop1,aBaB_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",aBab_"$((interval*counter))"_pop1,aBab_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  printf ",abab_"$((interval*counter))"_pop1,abab_"$((interval*counter))"_pop2" >> ColumnHeader_TL_BM_SC.csv
  counter=$(($counter+1))
done
printf "\n" >> ColumnHeader_TL_BM_SC.csv
echo "ColumnHeader_TL_BM_SC.csv created."
