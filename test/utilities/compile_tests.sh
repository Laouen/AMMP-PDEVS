for f in *.xml; do
	echo "Compiling model ${f}" 
	cd ../../model_generator
	python generate_model.py -f ../test/utilities/${f}
	cd ../
	make all
	make clean_all
	cd test/utilities
done