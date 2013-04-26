#!/usr/bin/zsh
if test ${1} = ""; then
  echo "usage:"
  echo "online_to_Cpp.sh start increment end"
  exit 1
fi 

for i in `seq ${1} ${2} ${3}`; do
  num=`printf "%06d" ${i}`
  grep -e '^1  1' onlinemeas.${num} | awk '{print $3, $4}' > Cpp.data.${num}.online
done
