#!/bin/bash
echo -e "\e[1;37mRunning unit tests:\e[0m"
logfile=build/tests.log
for i in build/*_tests
do
    if test -f $i
    then
        if $VALGRIND ./$i 2>> $logfile
        then
            echo $i PASS
        else
            echo -e "\e[31mERROR\e[0m in test $i: here's $logfile"
            echo "--------"
            tail $logfile
            exit 1
        fi
    fi
done
echo ""
