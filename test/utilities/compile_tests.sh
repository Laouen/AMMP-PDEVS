for f in *.xml; do
	echo "Compiling model ${f}" 
	pmgbp_generate_model -f ${f} -d ../../
	cd ../../
	make all
	make clean_all
	cd test/utilities
done