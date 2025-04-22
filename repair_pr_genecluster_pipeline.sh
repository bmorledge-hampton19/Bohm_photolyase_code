#!/bin/bash
# repair_genecluster_pipeline.sh
#Code adapted from code on Taylor lab github site

a="repair_"
export VAL=$1
export TWO=$2
export ZERO=$3
export REPAIR=$a${TWO}
mkdir ${REPAIR}

cd ${REPAIR}

#copy wig files to directory
cp ../${TWO}/${TWO}_dipy_bkgd_plus.wig ${TWO}_dipy_bkgd_plus.wig
cp ../${TWO}/${TWO}_dipy_bkgd_minus.wig ${TWO}_dipy_bkgd_minus.wig
cp ../${ZERO}/${ZERO}_dipy_bkgd_plus.wig ${ZERO}_dipy_bkgd_plus.wig
cp ../${ZERO}/${ZERO}_dipy_bkgd_minus.wig ${ZERO}_dipy_bkgd_minus.wig

# process data to get fraction remaining
printf "${TWO}_dipy_bkgd_plus.wig\n${TWO}_dipy_bkgd_minus.wig\n${ZERO}_dipy_bkgd_plus.wig\n${ZERO}_dipy_bkgd_minus.wig\n" | perl ../cpdrepair_plotorf_highres.pl >${REPAIR}_dipy_bkgd_tcr_15bpx2bins.txt

printf "${REPAIR}_dipy_bkgd_tcr_15bpx2bins.txt\n${VAL}\n" | perl ../arbnorm_clusters.pl

printf "${REPAIR}_dipy_bkgd_tcr_15bpx2bins_norm${VAL}.txt\n0.55\n" | perl ../arbcenter_clusters.pl

printf "${REPAIR}_dipy_bkgd_tcr_15bpx2bins_norm${VAL}_center0.55.txt" | perl ../orfplot_trxall.pl

printf "${REPAIR}_dipy_bkgd_tcr_15bpx2bins_norm${VAL}_center0.55_trxsorted.txt" | perl ../split_nts_ts.pl

perl ../format_CDT.pl <${REPAIR}_dipy_bkgd_tcr_15bpx2bins_norm${VAL}_center0.55_trxsorted_ts.txt >${REPAIR}_dipy_bkgd_tcr_15bpx2bins_norm${VAL}_center0.55_trxsorted_ts_cluster.cdt

perl ../format_CDT.pl <${REPAIR}_dipy_bkgd_tcr_15bpx2bins_norm${VAL}_center0.55_trxsorted_nts.txt >${REPAIR}_dipy_bkgd_tcr_15bpx2bins_norm${VAL}_center0.55_trxsorted_nts_cluster.cdt

