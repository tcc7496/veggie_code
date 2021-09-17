: 'shell script to batch process atmospheric correction. Must be run from folder with input files'

for file in $(ls -1); do
    /Users/taracunningham/Downloads/Sen2Cor-02.09.00-Darwin64/bin/L2A_Process $file --output_dir /Users/taracunningham/projects/dissertation/sen2processing/processing/l2a
    echo "Atmospheric correction comleted for $file"
done
echo "Atmospheric correction completed for all files"
