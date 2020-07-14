#!/bin/bash
exit_code=0
commit_id=${CI_BUILD_REF:0:8}

echo "Start to plot results of commit ID $commit_id.."
destfile=/y/hlab/AccurityTestData/gitlab_run/run_result.tsv

figure_output_folder=/y/hlab/AccurityTestData/gitlab_run/run_result_png
echo "Plotting result to $figure_output_folder ..."
cmd="test/plotResult.py -f ${destfile} -o ${figure_output_folder}"
echo ${cmd}
eval ${cmd}
plotExitCode=$?
if [ $plotExitCode != 0 ]; then
	echo "ERROR in plotting result."
	exit $plotExitCode
fi

echo "Copying $figure_output_folder content to ./ for artifact uploading ..."
cp -apr $figure_output_folder ./
echo "Copying run_result_png done."

echo "Done!"
exit ${exit_code}
