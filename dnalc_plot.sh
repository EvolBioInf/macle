#/bin/bash
PLOTFILE=plotdata.tmp
read keyword
if [ "$keyword" != "DNALC_PLOT" ]; then
  >&2 echo "ERROR: Unexpected input data to dnalc_plot.sh! Did you set the correct flag in dnalc?"
  exit
fi

read script
script=$(echo $script | sed 's/$PLOTFILE/'"$PLOTFILE/g")

cat | grep -ve '-1' > $PLOTFILE
script+="bind 'q' 'exit gnuplot'"

mkfifo $$.gnuplot-pipe
gnuplot -p <$$.gnuplot-pipe & pid=$! exec 3>$$.gnuplot-pipe
echo "$script" >&3
while [ -n "$(pidof gnuplot)" ]; do
    sleep .25s
done
exec 3>&-

rm $$.gnuplot-pipe
rm $PLOTFILE
