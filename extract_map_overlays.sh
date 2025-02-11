## Shell script to pull map overlay kml files from zip files

for infile in `ls *.zip`
do
	## construct filenames and pathways (safe_folder, outfilename, etc)
	infile_base=`echo $infile | cut -d "_" -f6`
	safe_folder=`echo $infile | sed 's/.zip/.SAFE/'`
	outfilename=${infile_base}_map_overlay.kml
	echo "Extracting: " $outfilename

	unzip -p $infile $safe_folder/preview/map-overlay.kml > $outfilename
	echo ""
done
