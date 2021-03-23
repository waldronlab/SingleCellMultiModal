SCMM="$HOME/gh/SingleCellMultiModal"

WIKI="$HOME/wiki/SingleCellMultiModal.wiki"

RVER="devel"

cd $SCMM

export R_LIBS_USER="/media/$USER/1D24A0EA4286043C1/bioc-$RVER/"

RCMD="$HOME/src/svn/r-$RVER/R/bin/R --no-save --no-restore-data"

$RCMD CMD INSTALL $SCMM

$RCMD -e "rmarkdown::render('inst/scripts/Contributing-Guidelines.Rmd', output_file = '$WIKI/Contributing-Guidelines.md')"

cd $WIKI

git diff

git pull origin master
git commit -am "update wiki"
git push origin master

