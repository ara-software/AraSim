# run AraSim with the veff setup file

./AraSim test/veff_ara3/ara3_setup 0 outputs/. 2>&1 | tee test/veff_ara3/veff_ara3_test_output.txt 

# then, do a comparison to check for consistency

EXPECTED_GLOBAL_PASS=1
EXPECTED_TOTAL_WEIGHT=5.20232e-15
EXPECTED_TOTAL_WEIGHT_SIGMA=0.0001
python3 test/check_sim.py test/veff_ara3/veff_ara3_test_output.txt $EXPECTED_GLOBAL_PASS $EXPECTED_TOTAL_WEIGHT $EXPECTED_TOTAL_WEIGHT_SIGMA