

CASEID=$1
ROOTPWD=sdfWER234SDF
CBIO_PATH=/home/gshuang555/cbioportal-docker-compose

echo "--- script started ---"
echo $ROOTPWD > sudo -s 

cd $CBIO_PATH
touch ./cbioApiStatus/${CASEID}.1.ready_to_start

cd $CBIO_PATH
python3 ./f2_cbioportal_preprocess.py ${CASEID}.ngs.csv
touch ./cbioApiStatus/${CASEID}.2.finish_running_etl_process

sleep 60s

cd $CBIO_PATH
docker-compose run cbioportal metaImport.py -u http://cbioportal:8080 -s study/data_upload/ -o

echo "-- End of script ---"
touch ./cbioApiStatus/${CASEID}.9.finish_uploading

exit
