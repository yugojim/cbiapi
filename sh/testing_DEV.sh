
CASEID=$1
ROOTPWD=sdfWER234SDF
CBIO_PATH=/Users/jxhuang555/Library/CloudStorage/OneDrive-WiAdvanceTechnologyCorporation/WiA_Prj/2308-VeteransHosp.cBio北榮基因/devops/src_cBioAPI/testing-docker
echo "--- script started ---"
echo $ROOTPWD > sudo -s 
cd $CBIO_PATH

touch ./cbioApiStatus/${CASEID}.1.running
touch ./cbioApiStatus/${CASEID}.2.etl

sleep 10s

echo "-- End of script ---"
touch ./cbioApiStatus/${CASEID}.9.finish

exit
