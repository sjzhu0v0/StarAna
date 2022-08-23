if [ $# -lt 1 ]; then
        get_file_list.pl -keys path,filename -cond filename~st_physics,trgsetupname=production_pAu200_2015,production=P16id,filetype=daq_reco_picoDst,storage!=HPSS -limit 0 -delim / -distinct > file.list
        sed -i 's/^/root\:\/\/xrdstar\.rcf\.bnl\.gov\:1095\//' file.list
else
        if [ ! -f file.list ]; then
                ./getfilelist.sh
        fi
        ./getrunid.sh file.list > getfilelist.temp1
        sed -n '1,'${1}'p' getfilelist.temp1 > getfilelist.temp2
        echo -n "{{"
        index=0
        for i in `cat getfilelist.temp2`
        do
                echo -n "${i},$index},{"
                index=$((index+1))
                grep $i file.list >> file.list${1}
        done
        rm getfilelist.temp*
fi