# run AraSim with the veff setup file

./AraSim test/veff_birefringence/birefringence_setup 0 outputs/. 2>&1 | tee test/veff_birefringence/veff_birefringence_test_output.txt 

# then, do a comparison to check for consistency

EXPECTED_GLOBAL_PASS=2
EXPECTED_TOTAL_WEIGHT=1.99728
EXPECTED_TOTAL_WEIGHT_SIGMA=0.0001
python3 test/check_sim.py test/veff_birefringence/veff_birefringence_test_output.txt $EXPECTED_GLOBAL_PASS $EXPECTED_TOTAL_WEIGHT $EXPECTED_TOTAL_WEIGHT_SIGMA