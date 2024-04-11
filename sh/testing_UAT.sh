
CASEID=$1
ROOTPWD=sdfWER234SDF
CBIO_PATH=/home/sdfWER234/cbioportal-docker-compose
echo "--- script started ---"
echo $ROOTPWD > sudo -s 
cd $CBIO_PATH

touch ./cbioApiStatus/${CASEID}.1.running
touch ./cbioApiStatus/${CASEID}.2.etl

sleep 10s

echo "-- End of script ---"
touch ./cbioApiStatus/${CASEID}.9.finish

exit
