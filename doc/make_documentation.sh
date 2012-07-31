#!/bin/bash

HOME_DIR=/home/users/alfthan
PUBLIC_HTML_DIR=${HOME_DIR}/public_html
DOXYGEN=/home/users/alfthan/doxygen-1.8.1.2/bin/doxygen


#lockfile, only one script at a time
LOCKFILE=${HOME_DIR}/.make_doc_is_working    


VLASIATOR_DIR=${HOME_DIR}/vlasiator/trunk
VLASIATOR_COMMITLOG=${PUBLIC_HTML_DIR}/vlasiator/commitlog.html
VLASIATOR_DOXYGENDOCS=${PUBLIC_HTML_DIR}/vlasiator/doc

DCCRG_DIR=${HOME_DIR}/vlasiator/dccrg
DCCRG_DOXYGENDOCS=${PUBLIC_HTML_DIR}/vlasiator/dccrg_doc

PHIPROF_DIR=${HOME_DIR}/vlasiator/phiprof
PHIPROF_DOXYGENDOCS=${PUBLIC_HTML_DIR}/vlasiator/phiprof_doc


if [ ! -e $LOCKFILE ]
then
touch $LOCKFILE

#first vlasiator 
    cd $VLASIATOR_DIR
#Update local copy of svn and store if it was already up-to-date
    isUpToDate=$( svn up | grep  "At revision" |wc -l )
    echo VLASIATOR up to date $isUpToDate
    if [ "$1" = "--force" ]
    then
	isUpToDate=0
    fi

    if [ ${isUpToDate} -ne 1 ] 
    then
#create svn log
	svn log --xml --verbose > log.xml                
	xsltproc -o $VLASIATOR_COMMITLOG  ../svn2cl-0.13/svn2html.xsl log.xml	
	sed  -i 's@QUESPACE-[0-9]*@<a href="https://agile.fmi.fi/browse/&">&</a>@g' $VLASIATOR_COMMITLOG  
	
#create doxygen docs
	$DOXYGEN
	rm -rf $VLASIATOR_DOXYGENDOCS
	cp -r doc/html $VLASIATOR_DOXYGENDOCS
	
#make files readable
	chmod -R a+rX ~/public_html/vlasiator
    fi


#then dccrg
    cd $DCCRG_DIR
#Update local copy of svn and store if it was already up-to-date
    isUpToDate=$( git pull |grep "Already up-to-date."|wc -l )
    echo DCCRG up to date $isUpToDate
    if [ "$1" = "--force" ]
    then
	isUpToDate=0
    fi

    if [ ${isUpToDate} -ne 1 ] 
    then
#no commitlog is prodcued at the moment
	
#create doxygen docs
	$DOXYGEN
	rm -rf $DCCRG_DOXYGENDOCS
	cp -r documentation/html $DCCRG_DOXYGENDOCS
	
#make files readable
	chmod -R a+rX $PUBLIC_HTML_DIR
    fi

#then dccrg
    cd $PHIPROF_DIR
#Update local copy of svn and store if it was already up-to-date
    isUpToDate=$( git pull |grep "Already up-to-date."|wc -l )
    echo PHIPROF up to date $isUpToDate

    if [ "$1" = "--force" ]
    then
	isUpToDate=0
    fi

    if [ ${isUpToDate} -ne 1 ] 
    then
#no commitlog is prodcued at the moment
	
#create doxygen docs
	$DOXYGEN
	rm -rf $PHIPROF_DOXYGENDOCS
	cp -r documentation/html $PHIPROF_DOXYGENDOCS
	
#make files readable
	chmod -R a+rX $PUBLIC_HTML_DIR
    fi

    
    rm $LOCKFILE
else
    echo "Warning, other $0 script already running!"
fi


