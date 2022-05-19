
dirname=$1
echo "Dir = $dirname"

touch $dirname.compression.ratio

> $dirname.compression.ratio

for eachfile in $dirname/*.h5
do
	echo $eachfile
	ratio=`h5dump -A  -p $eachfile  |grep COMPRESSION | tr -s ' ' | cut -d ' ' -f 4 | cut -d '(' -f 2 | cut -d ':' -f 1`
	eachfile_base=`basename $eachfile`
	echo "$eachfile_base $ratio" >> $dirname.compression.ratio
done

