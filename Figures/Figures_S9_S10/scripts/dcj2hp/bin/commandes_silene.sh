echo "ax"
java -Xmx512M -server -cp .:.. mcmc/DCJCoupledPathMCMC s ~/2023_silene/data_cooked/chr_ancestor_sex ~/2023_silene/data_cooked/chr_X 5000 100 500 20 0.01 > ~/2023_silene/data_cooked/silene_DCJ2HP_OUTPUT_ax
echo "ay"
java -Xmx512M -server -cp .:.. mcmc/DCJCoupledPathMCMC s ~/2023_silene/data_cooked/chr_ancestor_sex ~/2023_silene/data_cooked/chr_Y 5000 100 500 20 0.01 > ~/2023_silene/data_cooked/silene_DCJ2HP_OUTPUT_ay
echo "ac"
java -Xmx512M -server -cp .:.. mcmc/DCJCoupledPathMCMC s ~/2023_silene/data_cooked/chr_ancestor_speciation ~/2023_silene/data_cooked/chr_conica 5000 100 500 20 0.01 > ~/2023_silene/data_cooked/silene_DCJ2HP_OUTPUT_ac
echo "av"
java -Xmx512M -server -cp .:.. mcmc/DCJCoupledPathMCMC s ~/2023_silene/data_cooked/chr_ancestor_speciation ~/2023_silene/data_cooked/chr_vulgaris 5000 100 500 20 0.01 > ~/2023_silene/data_cooked/silene_DCJ2HP_OUTPUT_av
echo "aa"
java -Xmx512M -server -cp .:.. mcmc/DCJCoupledPathMCMC s ~/2023_silene/data_cooked/chr_ancestor_speciation ~/2023_silene/data_cooked/chr_ancestor_sex 5000 100 500 20 0.01 > ~/2023_silene/data_cooked/silene_DCJ2HP_OUTPUT_aa
