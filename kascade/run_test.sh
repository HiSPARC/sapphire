echo "Removing stale files"
rm -f generator.h5
rm -f generator-buffer.sql.gz
rm -f generator-kascade.dat.gz

echo "Generating events"
python event_generator.py

echo "Dumping database and gzipping files"
mysqldump5 -u buffer --password=Buffer4hisp! buffer > generator-buffer.sql
gzip generator-buffer.sql
gzip generator-kascade.dat

echo "Reading database dump into hd5 file"
cd read_sqldump
python read_sqldump.py

echo "Reading kascade file into hdf5 file"
cd ..
python read_kascade.py

echo "Analysing hisparc / kascade data"
python analysis_kascade.py
