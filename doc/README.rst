Sphinx documentation builds in gh-pages branch
==============================================

This contains Sphinx documentation for SAPPHiRE.  On Github, each
repository can have a gh-pages branch to contain html pages for your
project.  SAPPHiRE's documentation lives in
http://docs.hisparc.nl/sapphire/.

It is not straight-forward to keep the documentation in your master branch
while publishing the builds in the gh-pages branch.  From a few sources,
I've gathered the following approach.


Easy deployment
---------------

The following commands are added to the Makefile as the gh-pages target.
To update the docs simply run the following command from the repository
root::

   $ make gh-pages


Creating the gh-pages branch
----------------------------

The gh-pages branch should be separate from the normale code branches.
It will be created as an 'orphaned' branch, which will not contain any
of the master branch history. Follow these instructions::

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
    $ git clean -dxf
    $ git checkout HEAD .nojekyll
    $ git checkout master doc sapphire
    $ make -C doc/ html
    $ mv -fv doc/_build/html/* .
    $ rm -rf doc/ sapphire/
    $ git add -A
    $ git commit -m "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`"
    $ git checkout master


Sources
-------

* https://github.com/matthew-brett/gh-sphinx-template
* http://blog.nikhilism.com/2012/08/automatic-github-pages-generation-from.html
* https://help.github.com/articles/creating-project-pages-manually
