.PHONY: test
test: unittests flaketest doctest

.PHONY: unittests
unittests:
	coverage run -m unittest

.PHONY: linttest
linttest:
	ruff check .
	ruff format --check .

.PHONY: lintfix
lintfix:
	ruff check --fix-only .
	ruff format .
	ruff check --fix-only .
	ruff format .

.PHONY: doctest
doctest:
	sphinx-build -anW doc doc/_build/html

.PHONY: update_data
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
