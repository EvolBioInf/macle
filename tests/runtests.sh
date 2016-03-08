#!/bin/bash
logfile=build/tests.log
echo -e "\e[1;37mRunning unit tests:\e[0m"
for i in build/*_tests
do
    if test -f $i
    then
        vglog="/tmp/valgrind-$(basename $i).log"
        if [ -n "$VALGRIND" ]; then
          VG="$VALGRIND --error-exitcode=2 --log-file=$vglog"
        else
          VG=""
        fi
        $VG ./$i 2>> $logfile
        code=$?
        if [ $code == 0 ]; then
            echo $i PASS
        elif [ $code == 1 ]; then
            echo -e "\e[31mERROR\e[0m in test $i: here's $logfile"
            echo "--------"
            tail $logfile
            exit 1
        else
            echo -e "\e[31mMEMORY LEAK\e[0m in test $i: here's $vglog"
            echo "--------"
            tail $vglog
            exit 2
        fi
    fi
done
echo ""
