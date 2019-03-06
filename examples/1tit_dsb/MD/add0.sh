#/bin/bash
here=$(pwd)
file=$(ls -1 table_b*.xvg)
for i in $file;do
echo $i
#sed -i -e '1i      0.002 1.000000e+36 1.000000e+29' $i
gsed -i -e '1i            0.004 0.000000 0.000000' $i
gsed -i -e '1i            0.002 0.000000 0.000000' $i
gsed -i -e '1i            0.000 0.000000 0.000000' $i

done
