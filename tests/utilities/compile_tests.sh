for f in *.xml; do
	echo "Compiling model ${f}" 
	cd ../../model_generator
	python generate_model.py -f ../tests/utilities/${f}
	cd ../
	make all
	make clean_all
	cd tests/utilities
done