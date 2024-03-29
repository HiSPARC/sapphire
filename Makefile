.PHONY: test unittests flaketest doctest update_data

test: unittests flaketest doctest

unittests:
	coverage run -m unittest

flaketest:
	flake8 sapphire

doctest:
	sphinx-build -anW doc doc/_build/html

update_data:
ifeq ($(strip $(shell git status --porcelain | wc -l)), 0)
	@echo "Updating local data. Creating test data to match local data and committing."
	update_local_data
	create_and_store_test_data
	git add -A
	git commit -m "Updated local and test data for timestamp `date +%s`"
	@echo "Done."
else
	$(error Working tree is not clean, please commit all changes.)
endif
