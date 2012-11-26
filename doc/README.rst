Sphinx documentation builds in gh-pages branch
==============================================

This directory contains Sphinx documentation for SAPPHiRE.  On Github,
each repository can have a gh-pages branch to contain html pages for your
project.  SAPPHiRE's documentation lives in
http://hisparc.github.com/sapphire/.

It is not straight-forward to keep the documentation in your master branch
while publishing the builds in the gh-pages branch.  From a few sources,
I've gathered the following approach.


Creating the gh-pages branch
----------------------------

Follow these instructions::

    $ git checkout --orphan gh-pages
    $ git rm -rf .
    $ touch .nojekyll
    $ git add .
    $ git commit -m "Initial commit"


Building gh-pages docs
----------------------

Issue the following commands from the repository root::

    $ git checkout gh-pages
    $ git rm -rf .
    $ git checkout HEAD .nojekyll
    $ git checkout master doc
    $ make -C doc html
    $ mv -fv doc/_build/html/* .
    $ rm -rf doc/
    $ git add -A ??
    $ git commit -m "Generated gh-pages for `git log master -1 --pretty=short
    --abbrev-commit`"
    $ git checkout master


Easy deployment
---------------

The commands from the previous section are added to the Makefile as the
gh-pages target.


Sources
-------

https://github.com/matthew-brett/gh-sphinx-template
http://blog.nikhilism.com/2012/08/automatic-github-pages-generation-from.html
https://help.github.com/articles/creating-project-pages-manually
