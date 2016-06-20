.PHONY: gh-pages test

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

test:
	python setup.py test
	flake8 --ignore=Z --exclude=sapphire/transformations/geographic.py,sapphire/tests/,sapphire/corsika/qsub_corsika.py sapphire
	flake8 --ignore=E501 sapphire/tests/
	flake8 --ignore=E501 sapphire/corsika/qsub_corsika.py
	flake8 --exit-zero --ignore=Z sapphire/transformations/geographic.py
	sphinx-build -anW doc doc/_build/html
