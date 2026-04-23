#!/bin/zsh

PR=$1  # Pull request ID
SOURCEDIR=~uganse/vlasiator
GITHUB_USERNAME=ursg

# Check command line argument
if [ -z $PR ]; then
   echo "Usage: ./build_and_test_pull_request.sh <PR number>"
   exit 1
fi

# Proper shell error handling and cleanup
set -euo pipefail
function cleanup {
   echo "### Script failed. Cleaning up and aborting. ###"
   git checkout $OLDBRANCH
   popd
}
trap cleanup EXIT

# Go to source dir, check out correct branch
echo "--- Going to $SOURCEDIR ---"
pushd $SOURCEDIR

OLDBRANCH="${$(git symbolic-ref HEAD)/refs\/heads\//}"
echo ">>> [Old branch was $OLDBRANCH]"
if [[ `git status --porcelain --untracked-files=no` ]]; then
   echo "Your vlasiator source dir $SOURCEDIR has local changes."
   echo "Stash or commit them, then retry. Aborting."
   popd
   exit 1
fi

echo "--- Determining pull request $PR's git HEAD ---"
HEAD=$(curl -s -u $GITHUB_USERNAME:`cat github_token` -X GET https://api.github.com/repos/fmihpc/vlasiator/pulls/$PR  | jq -r .head.sha)
echo ">>> It's $HEAD."

echo -e "{\"body\": \"Automatic testpackage output for pull request #$PR.\n\nBased on commit $HEAD.\n"  > /tmp/githubcomment_$HEAD.txt
echo -e "Process started by user \`\`\`$USER\`\`\` on host \`\`\`$HOST\`\`\` at `date --rfc-822`.\n" >> /tmp/githubcomment_$HEAD.txt

git fetch origin -q
echo "--- Checking out $HEAD ---"
git checkout -q $HEAD

echo "--- Compiling Vlasiator ---"
function cleanup {
   echo "### Compilation failed. See testpackage_build_$HEAD.log for details. Cleaning up and aborting. ###"
   echo -e ":x: Compilation failed.\n<details><summary>Compilation error log</summary>\n<p>\n\`\`\`" >> /tmp/githubcomment_$HEAD.txt
   cat testpackage_build_$HEAD.log >> /tmp/githubcomment_$HEAD.txt
   echo -e "\`\`\`\n</p>\n</details>\n" >> /tmp/githubcomment_$HEAD.txt
   git checkout $OLDBRANCH
   popd
}

make -s clean 
srun -M vorna --job-name tp_compile --interactive --nodes=1 -n 1 -c 16 --mem=40G -p test -o testpackage_build_$HEAD.log -t 0:10:0 make -s -j 16 testpackage tools

WARNINGS=`grep -c warning: testpackage_build_$HEAD.log`
ERRORS=`grep -c error: testpackage_build_$HEAD.log || true`

# Note we don't need to output error count here, because that would anyway have
# stopped the script already.
echo ">>> Compilation succeeded with $WARNINGS warnings."

echo -e ":heavy_check_mark: Compilation succeeded with $WARNINGS warnings.\n" >> /tmp/githubcomment_$HEAD.txt

# Go back to testpackage dir
git checkout -q $OLDBRANCH
popd

cp $SOURCEDIR/vlasiator .
cp $SOURCEDIR/vlsvdiff_DP .


function cleanup {
   echo "### Running testpackage failed. Cleaning up and aborting. ###"
   echo -e ":x: Running testpackage failed.\nError log:\n\`\`\`" >> /tmp/githubcomment_$HEAD.txt
   cat testpackage_run_$HEAD.log >> /tmp/githubcomment_$HEAD.txt
   echo -e "\`\`\`\n" >> /tmp/githubcomment_$HEAD.txt
}

# Run it actually.
echo "--- Now running testpackage ---"
sbatch -W -o ./testpackage_run_$HEAD.log ./small_test_vorna.sh

echo ">>> Testpackage run successfully!"
echo -e ":heavy_check_mark: Testpackage ran successfully.\n" >> /tmp/githubcomment_$HEAD.txt

echo "--- Parsing testpackage results ---"
CURRENT_TEST=""
MAXERR=0.  # Absolute error
MAXREL=0.  # Relative error
MAXVAR=""  # Variable with said error
while read line; do
   echo "$line" >> /tmp/curtest.log

   if echo "$line" | grep -q '^running '; then
      if [ ! -z $CURRENT_TEST ]; then
         if [ 1 -eq $(( $MAXERR == 9999 )) ]; then
            echo -e "\e[31;1m>>> Test $CURRENT_TEST failed to run.\e[0m"
            echo -e ":red_circle: $CURRENT_TEST: Failed to run or died with an error."  >> /tmp/githubcomment_$HEAD.txt
         elif [ 1 -eq $(( $MAXERR > 0 )) ]; then
            echo -e "\e[33;1m>>> Nonzero diffs in test $CURRENT_TEST\e[0m"
            echo -e ":large_orange_diamond: $CURRENT_TEST: Nonzero diffs: \`$MAXVAR\` has absolute error $MAXERR, relative error $MAXREL."  >> /tmp/githubcomment_$HEAD.txt
         else
            echo -e "\e[32m>>> Zero diffs in test $CURRENT_TEST\e[0m"
            echo -e ":heavy_check_mark: $CURRENT_TEST: Zero diffs." >> /tmp/githubcomment_$HEAD.txt 
         fi

      fi
      CURRENT_TEST=`echo "$line" | tail -c +9`
      MAXERR=0.
      echo "$line" > /tmp/curtest.log

   elif echo "$line" | egrep -q '^(fg_|vg_|proton/vg_)\w+\s+[0-9.e+-]+\s+[0-9.e+-]+\s*$'; then
      VAR=`echo "$line" | awk '{print $1}'`
      ABS=`echo "$line" | awk '{print $2}'`
      REL=`echo "$line" | awk '{print $3}'`

      if [ 1 -eq $(( $ABS > $MAXERR )) ]; then
         MAXERR=$ABS
         MAXVAR=$VAR
         MAXREL=$REL
      fi

   elif echo "$line" | egrep -qi "aborted|segmentation|failed"; then
      MAXERR=9999
      MAXVAR="Failure"
   fi
done < testpackage_run_$HEAD.log
if [ 1 -eq $(( $MAXERR == 9999 )) ]; then
   echo -e "\e[31;1m>>> Test $CURRENT_TEST failed to run.\e[0m"
   echo -e ":red_circle: $CURRENT_TEST: Failed to run or died with an error."  >> /tmp/githubcomment_$HEAD.txt
elif [ 1 -eq $(( $MAXERR > 0 )) ]; then
   echo -e "\e[33;1m>>> Nonzero diffs in test $CURRENT_TEST\e[0m"
   echo -e ":large_orange_diamond: $CURRENT_TEST: Nonzero diffs: \`$MAXVAR\` has absolute error $MAXERR, relative error $MAXREL."  >> /tmp/githubcomment_$HEAD.txt
else
   echo -e "\e[32m>>> Zero diffs in test $CURRENT_TEST\e[0m"
   echo -e ":heav_check_mark: $CURRENT_TEST: Zero diffs." >> /tmp/githubcomment_$HEAD.txt
fi

echo "\nFull testpackage run output is available in \``pwd`/testpackage_run_$HEAD.log\` on host \`$HOST\`" >> /tmp/githubcomment_$HEAD.txt
echo -n "\"}" >> /tmp/githubcomment_$HEAD.txt

# Upload to github.
#echo "Github comment created at /tmp/githubcomment_$HEAD.txt"

# The uploaded comment needs its newlines quoted
sed -zi 's/\n/\\r\\n/g' /tmp/githubcomment_$HEAD.txt

echo "--- Sending result to github ---"
curl -s -u $GITHUB_USERNAME:`cat $SOURCEDIR/github_token` -X POST -H "Accept: application/vnd.github.v3+json" https://api.github.com/repos/fmihpc/vlasiator/issues/$PR/comments --data-binary "@/tmp/githubcomment_$HEAD.txt"

# Clean up
echo "--- DONE. ---"
trap - EXIT

