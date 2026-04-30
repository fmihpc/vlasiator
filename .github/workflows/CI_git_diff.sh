
#!/usr/bin/env bash
git_diff=$(git diff --name-only origin/$@...)
declare -i full_count=$(echo $git_diff | wc -w)
declare -a unallowed_strings=(".cpp\b" ".hpp\b" ".c\b" ".sh\b" )

#Check if git diff contains any changes to the actual code
for unallowed_string in ${unallowed_strings[@]}
do
  if [[ $(echo $git_diff | grep -Po $unallowed_string) ]]; then
    exit 0
  fi
done

declare -a file_check=("\w+.md" "INSTALL" "CITATION.cff" "Doxyfile." "\bdoc/" "\b.vscode/" "\bdoxytmp/" "\bscaling/" "\bvisit_variables/")
declare -i allowed_count=0

#checking if only allowed files are changed. 
for allowed_string in ${file_check[@]}
do
  diff="$(echo $git_diff | grep -Po $allowed_string | wc -w)"
  allowed_count=$(( $allowed_count + $diff ))
done
if [[ $allowed_count == $full_count ]]; then
  echo "no_test"
fi

