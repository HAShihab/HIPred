
rm -rf ./tmp
mkdir ./tmp

stdbuf -i0 -o0 -e0 Rscript big_matrix.R ./data/Annotations/median.csv.gz ./tmp/hipred_median.csv ./tmp/hipred_median.txt

# how does mean/max perform ...
stdbuf -i0 -o0 -e0 Rscript big_matrix.R ./data/Annotations/mean.csv.gz ./tmp/hipred_mean.csv ./tmp/hipred_mean.txt
stdbuf -i0 -o0 -e0 Rscript big_matrix.R ./data/Annotations/max.csv.gz ./tmp/hipred_max.txt ./tmp/hipred_max.txt


## TODO - append code for MKL and Stacking (after review)


# derive genome-wide predictions
Rscript predict_genome.R

# make a prediction/correlation matrix
python predict_merge.py
python predict_corr.py

python benchmark.py > ./tmp/benchmark.txt

convert ./tmp/OMIM\ HI.png -pointsize 42 label:'A' +swap -append ./tmp/a.png && mv ./tmp/a.png ./tmp/OMIM\ HI.png
convert ./tmp/OMIM\ HI\ de\ novo.png -pointsize 42 label:'B' +swap -append ./tmp/b.png && mv ./tmp/b.png ./tmp/OMIM\ HI\ de\ novo.png
convert ./tmp/MGI\ Lethality.png -pointsize 42 label:'C' +swap -append ./tmp/c.png && mv ./tmp/c.png ./tmp/MGI\ Lethality.png
convert ./tmp/MGI\ Seizure.png -pointsize 42 label:'D' +swap -append ./tmp/d.png && mv ./tmp/d.png ./tmp/MGI\ Seizure.png
convert ./tmp/ASD1.png -pointsize 42 label:'E' +swap -append ./tmp/e.png && mv ./tmp/e.png ./tmp/ASD1.png
convert ./tmp/ASD2.png -pointsize 42 label:'F' +swap -append ./tmp/f.png && mv ./tmp/f.png ./tmp/ASD2.png

# bias check, wrt. well-studied genes
for x in `zcat ./data/Annotations/mean.csv.gz | cut -d "," -f 1 | grep -v "^id"`; do echo "$x" `./edirect/esearch -db pubmed -query "$x" | ./edirect/xtract -pattern ENTREZ_DIRECT -element Count`; done | tee ./tmp/PubMed_Counts.txt

Rscript bias_check.R
