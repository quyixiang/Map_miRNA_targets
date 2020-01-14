file_path=/Users/quyixiang/Desktop/get_sequence
result_path=/Users/quyixiang/Documents/result


python $file_path/miRNA-csv_extraction.py $1

mv $result_path/$1/*.csv $result_path/$1/$1.csv

if [ ! -f "$result_path/$1/$1.csv" ]; then
    rm -r "$result_path"/"$1"
else

    cat $result_path/$1/$1.csv | tail -n +15 | awk 'BEGIN{FS=",";} {print $1 "," $5}' > $result_path/$1/$1_csv_result.txt
    test=`cat $result_path/$1/$1_csv_result.txt | awk 'NR==1' | awk -F ',' '{if($2 < 0.7) print "0"; else print "1"}'`
    if [ "$test" -eq "0" ]; then
        rm -r "$result_path"/"$1"
    else
        mkdir $result_path/$1/sequence_result

        for line in `cat $result_path/$1/$1_csv_result.txt`

        do
        line1=`echo "$line" | awk 'BEGIN{FS=",";} {print $1}'`
        line2=`echo "$line" | awk 'BEGIN{FS=",";} {print $2}'`
        
        hold="0.7"
        state=`echo "$line2<$hold" | bc`
        #不等于1就说明大于等于0.7
        if [ "$state" != "1" ]; then
            

            resu=`python $file_path/esearch.py $line1`

            python $file_path/efetch.py $resu > $result_path/$1/sequence_result/$line1.txt

            sleep 0.2

        fi

        done
    fi
fi



