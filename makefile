clean:
	rm -rf __pycache__
	rm -rf */__pycache__

docclean:
	rm -rf doc/build/*

veryclean: docclean, clean
