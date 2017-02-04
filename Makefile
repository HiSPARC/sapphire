.PHONY: test unittests flaketest docstest gh-pages update_data

test: unittests flaketest doctest

unittests:
	coverage run --branch --source=sapphire setup.py test

flaketest:
	flake8 --ignore=Z --exclude=sapphire/transformations/geographic.py,sapphire/tests/,sapphire/corsika/qsub_corsika.py sapphire
	flake8 --ignore=E501 sapphire/tests/
	flake8 --ignore=E501 sapphire/corsika/qsub_corsika.py
	flake8 --exit-zero --ignore=Z sapphire/transformations/geographic.py

doctest:
	sphinx-build -anW doc doc/_build/html

gh-pages:
ifeq ($(strip $(shell git status --porcelain | wc -l)), 0)
	git checkout gh-pages
	git rm -rf .
	git clean -dxf
	git checkout HEAD .nojekyll
	git checkout master doc sapphire
	make -C doc/ html
	mv -fv doc/_build/html/* .
	rm -rf doc/ sapphire/
	git add -A
	git commit -m "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`"
	git checkout master
else
	$(error Working tree is not clean, please commit all changes.)
endif

update_data:
ifeq ($(strip $(shell git status --porcelain | wc -l)), 0)
	@echo "Updating local data. Creating test data to match local data and committing."
	#sapphire/data/update_local_data
	sapphire/tests/create_and_store_test_data
	git add -A
	git commit -m "Updated local and test data for `date +%s`"
	@echo "Done."
else
	$(error Working tree is not clean, please commit all changes.)
endif
