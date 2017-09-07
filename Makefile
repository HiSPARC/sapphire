.PHONY: test unittests flaketest doctest update_data

test: unittests flaketest doctest

unittests:
	coverage run --branch --source=sapphire setup.py test

flaketest:
	flake8 --exclude=sapphire/tests/ sapphire
	flake8 --ignore=E501 sapphire/tests/

doctest:
	sphinx-build -anW doc doc/_build/html

update_data:
ifeq ($(strip $(shell git status --porcelain | wc -l)), 0)
	@echo "Updating local data. Creating test data to match local data and committing."
	sapphire/data/update_local_data
	sapphire/tests/create_and_store_test_data
	git add -A
	git commit -m "Updated local and test data for timestamp `date +%s`"
	@echo "Done."
else
	$(error Working tree is not clean, please commit all changes.)
endif
