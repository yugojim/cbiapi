

CASEID=$1
ROOTPWD=EvanChou123456789
echo "--- script started ---"
echo $ROOTPWD > sudo -s 
touch ./cbioApiStatus/${CASEID}.1.running

cd /home/EvanChou/cbioportal-docker-compose
python3 ./cbioportal_preprocess.py
touch ./cbioApiStatus/${CASEID}.2.etl

sleep 2s

cd /home/EvanChou/cbioportal-docker-compose
docker-compose run cbioportal metaImport.py -u http://cbioportal:8080 -s study/data_upload/ -o

echo "-- End of script ---"
touch ./cbioApiStatus/${CASEID}.9.finish

exit
