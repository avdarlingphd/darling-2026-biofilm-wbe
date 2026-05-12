cat /n/netscratch/hhealy_lab/avdarling/kraken_coverage/*/coverage_metrics.tsv | grep -v "^sample" | awk '
BEGIN {total=0; breadth1=0; breadth10=0; depth_sum=0; detected=0; breadth1_sum=0; breadth10_sum=0}
{
    total++
    if ($7 > 0) {
        breadth1++
        breadth1_sum += $7
        depth_sum += $9
        detected++
    }
    if ($8 > 0) {
        breadth10++
        breadth10_sum += $8
    }
}
END {
    print "Total combinations: " total
    print ">=1x coverage: " breadth1 " (" int(breadth1/total*100) "%)"
    print ">=10x coverage: " breadth10 " (" int(breadth10/total*100) "%)"
    print "Mean breadth_1x (detected only): " breadth1_sum/detected
    print "Mean breadth_10x (detected only): " breadth10_sum/breadth10
    print "Mean depth (detected only): " depth_sum/detected
}'
