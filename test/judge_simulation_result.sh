#!/bin/bash
exit_code=0
commit_id=${CI_COMMIT_SHA:0:8}

data=$1
cp_judge_flag=$2
mcp_judge_flag=$3
max_purity_delta=0.05
max_ploidy_delta=0.2
min_cp_recall=0.7
min_cp_precision=0.7
min_mcp_recall=0.7
min_mcp_precision=0.7
current_result_path=${data}/result_${commit_id}*.tsv

echo "Start to judge result of commit ID $commit_id of $data.."
cat ${current_result_path}

while read line; do
    echo $line;
    true_purity=`echo ${line}|awk '{print $3}'`
    true_ploidy=`echo ${line}|awk '{print $4}'`
    echo "Assessing one result with true purity=${true_purity}, ploidy=${true_ploidy}"
    accurity_purity=`echo ${line}|awk '{print $5}'`
    accurity_ploidy=`echo ${line}|awk '{print $6}'`
    echo "Purity estimate is $accurity_purity, ploidy estimate is $accurity_ploidy."
    if [ "z$accurity_purity" ==  "z" ]; then
        echo "ERROR: accurity_purity is empty."
        exit 2
    fi
    if [ "z$accurity_ploidy" ==  "z" ]; then
        echo "ERROR: accurity_ploidy is empty."
        exit 2
    fi

    delta_purity=`echo ${line}|awk '{print $7}'`
    delta_ploidy=`echo ${line}|awk '{print $8}'`
    cp_recall=`echo ${line}|awk '{print $9}'`
    cp_precision=`echo ${line}|awk '{print $10}'`
    mcp_recall=`echo ${line}|awk '{print $11}'`
    mcp_precision=`echo ${line}|awk '{print $12}'`
    echo "CNV recall is $cp_recall, CNV precision is $cp_precision."
    echo "major allele copy number recall is $mcp_recall, major allele copy number precision is $mcp_precision."
    if [ "z$cp_recall" ==  "z" ]; then
        echo "ERROR: cp_recall is empty."
        exit 2
    fi
    if [ "z$cp_precision" ==  "z" ]; then
        echo "ERROR: cp_precision is empty."
        exit 2
    fi
    if [ "z$mcp_recall" ==  "z" ]; then
        echo "ERROR: mcp_recall is empty."
        exit 2
    fi
    if [ "z$mcp_precision" ==  "z" ]; then
        echo "ERROR: mcp_precision is empty."
        exit 2
    fi

    if [ `echo "${delta_purity} > ${max_purity_delta}"|bc` -ne 0 ] ; then
        echo "Purity estimate ${accurity_purity} deviated from the truth ${true_purity} by ${delta_purity}, more than 0.05."
        exit_code=1
    fi
    if [ `echo "${delta_ploidy} > ${max_ploidy_delta}"|bc` -ne 0 ]; then
        echo "Ploidy estimate ${accurity_ploidy} deviated from the truth ${true_ploidy} by ${delta_ploidy}, more than 0.2."
        exit_code=1
    fi
    if [ ${cp_judge_flag} != 0 -a `echo "${cp_recall} < ${min_cp_recall}"|bc` -ne 0 ]; then
        echo "CNV Recall is below ${min_cp_recall}: ${cp_recall}"
        exit_code=1
    fi
    if [ ${cp_judge_flag} != 0 -a `echo "${cp_precision} < ${min_cp_precision}"|bc` -ne 0 ]; then
        echo "CNV precision is below ${min_cp_precision}: ${cp_precision}"
        exit_code=1
    fi
    if [ ${mcp_judge_flag} != 0 -a `echo "${mcp_recall} < ${min_mcp_recall}"|bc` -ne 0 ]; then
        echo "CNV Recall is below ${min_mcp_recall}: ${mcp_recall}"
        exit_code=1
    fi
    if [ ${mcp_judge_flag} != 0 -a `echo "${mcp_precision} < ${min_mcp_precision}"|bc` -ne 0 ]; then
        echo "CNV precision is below ${min_mcp_precision}: ${mcp_precision}"
        exit_code=1
    fi
done < <(tail -n1 ${current_result_path})

echo "exit code is ${exit_code}."
exit ${exit_code}
