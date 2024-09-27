# run AraSim with the veff setup file

# ./AraSim test/veff_pa/PA_setup 0 outputs/. 2>&1 | tee test/veff_pa/veff_pa_test_output.txt 

# then, do a comparison to check for consistency

EXPECTED_GLOBAL_PASS=11
EXPECTED_TOTAL_WEIGHT=4.85117
EXPECTED_TOTAL_WEIGHT_SIGMA=0.0001
python3 test/check_sim.py test/veff_pa/veff_pa_test_output.txt $EXPECTED_GLOBAL_PASS $EXPECTED_TOTAL_WEIGHT $EXPECTED_TOTAL_WEIGHT_SIGMA