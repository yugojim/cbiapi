
CASEID=$1
ROOTPWD=sdfWER234SDF
CBIO_PATH=/home/gshuang555/cbioportal-docker-compose
echo "--- script started ---"
echo $ROOTPWD > sudo -s 
cd $CBIO_PATH

touch ./cbioApiStatus/${CASEID}.1.running

sleep 10s

echo "-- End of script ---"
touch ./cbioApiStatus/${CASEID}.9.finish

exit
