### Step 1

```bash
cat <RNAfile>.csv | tail -n +15 | awk 'BEGIN{FS=",";} {print$1}' > result.txt
```

##### The purpose of this step is to extract the name of the gene we need. And the name of the genes will be saved in one file whose name is result.txt

### Step 2

##### Make sure that the files 'efetch.py', 'research.py' and 'readline.sh' are saved in one folder. Just run ./readline.sh and it will work.



If it didn't work for you, you should install biopython package.(Just run ``pip install biopython``.)
