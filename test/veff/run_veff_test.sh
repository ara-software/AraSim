# run AraSim with the veff setup file

./AraSim test/veff/setup_veff_test.txt 0 outputs/. 2>&1 | tee test/veff/veff_test_output.txt 

# then, do a comparison to check for consistency

python3 test/veff/check_veff.py test/veff/veff_test_output.txt