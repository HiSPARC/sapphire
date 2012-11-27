.PHONY: gh-pages

gh-pages:
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
