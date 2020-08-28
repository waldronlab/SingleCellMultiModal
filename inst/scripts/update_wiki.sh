SCMM="$HOME/gh/SingleCellMultiModal"

WIKI="$HOME/wiki/SingleCellMultiModal.wiki"

cd $SCMM

export R_LIBS_USER="/media/$USER/1D24A0EA4286043C1/bioc-devel/"

RDEV="$HOME/src/svn/r-release/R/bin/R --no-save --no-restore-data"

$RDEV CMD INSTALL $SCMM

$RDEV -e "rmarkdown::render('inst/scripts/Contributing-Guidelines.Rmd', output_file = '$WIKI/Contributing-Guidelines.md')"

cd $WIKI

git diff

git pull origin master
git commit -am "update wiki"
git push origin master

