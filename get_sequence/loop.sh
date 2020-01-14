file_path=/Users/quyixiang/Desktop/get_sequence

for line in `cat $file_path/miRNA_name.txt`
do
line=${line%?}
bash $file_path/readline.sh $line

done