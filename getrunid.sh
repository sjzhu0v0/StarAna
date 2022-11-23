sed '{s/^.*st_ssdmb_//g}' < ${1} > getrunid.temp1
sed '{s/_raw.*$//g}' < getrunid.temp1 > getrunid.temp2
sed '{s/adc_//g}' < getrunid.temp2 | sort | uniq

rm getrunid.temp*
