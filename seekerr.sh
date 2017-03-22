#!/bin/bash
# 이름: seekerr.sh
# 목적: bsub 를 이용해 submit 한 job 이 모두 끝난 후에 결과물에 이상이 없는지 확인하고 이상이 있는 파일의 목록을 생성한다.
# 사용법: $ ./seekerr.sh setX
# 저자: 장우영(16.8.23)

echo "Starting error checker..."
PWD=`pwd`
OUTPUTDIR=/eos/ams/user/w/wyjang/He3He4/Output
LOGDIR=${PWD}/Output/$1/log
ERRDIR=${PWD}/Output/$1/err
echo "$0 $1"
# 파일리스트로부터 전체 파일 목록을 가져와서 실제로 생성된 파일의 목록과 비교한다.
echo "[1/12] Making file-list by : $ ls FileList/$1 > list.tmp"
ls ${PWD}/FileList/$1 > list.tmp
echo "Done."
echo "[2/12] Making output-list by : $ ls ${OUTPUTDIR}/$1/*.root > outlist.tmp"
ls ${OUTPUTDIR}/$1/*.root > outlist.tmp
echo "Done."
#cat outlist.tmp | sed -e "s#${OUTPUTDIR}/set[0-9]\{1,\}/##g" -e "s#.root##g" > outlist2.tmp
cat outlist.tmp | sed -e "s#${OUTPUTDIR}/$1/##g" -e "s#.root##g" > outlist2.tmp
# 정규표현식 해설 : set[0-9]\{1,\}
# set 은 단순 문자열이고 [0-9]는 임의의 숫자를 뜻하는데 뒤에 {1,}으로 인해 한 자리 이상의 숫자를 가리킨다.
echo "[3/12] Comparing file-list and output-list to build error-list"
diff list.tmp outlist2.tmp | grep [0-9]\{10\} > tmp_errorlist # 누락된 파일들을 에러목록에 추가.
echo "Done."

echo "[4/12] Removing temporary files..."
echo "$ rm -rf list.tmp"
rm -rf list.tmp
echo "$ rm -rf outlist.tmp"
rm -rf outlist.tmp
echo "Done"

# 에러파일 목록과 데이터파일 목록 간에 차이가 있는지 확인한다.
echo "[5/12] Making error-list by : $ ls ${ERRDIR}/*.err | sed -e "s#${ERRDIR}/##g" -e "s#.err##g" > errlist.tmp"
ls ${ERRDIR}/*.err | sed -e "s#${ERRDIR}/##g" -e "s#.err##g" > errlist.tmp
echo "Done"
echo "[6/12] Comparing error-list with file-list"
diff errlist.tmp outlist2.tmp | grep [0-9]\{10\} >> tmp_errorlist
echo "Done"
rm -rf errlist.tmp

# 로그파일 목록과 데이터파일 목록 간에 차이가 있는지 확인한다.
echo "[7/12] Making log-list by : $ ls ${LOGDIR}/*.log | sed -e "s#${LOGDIR}/##g" -e "s#.log##g" > loglist.tmp"
ls ${LOGDIR}/*.log | sed -e "s#${LOGDIR}/##g" -e "s#.log##g" > loglist.tmp
echo "Done"
echo "[8/12] Comparing log-list with file-list"
diff loglist.tmp outlist2.tmp | grep [0-9]\{10\} >> tmp_errorlist
echo "Done"
rm -rf loglist.tmp

# 배드런으로 확인된 파일들을 검색하여 목록을 작성한다.
# 에러가 난 결과 파일을 검색하여 목록을 작성한다.
outfile="errorlist_$1"
cd ${OUTPUTDIR}/$1
# ================ Begin of Working in ${OUTPUTDIR}/$1 ===================================
#grep 'session not found' ./err/* > tmp_errorlist
# 정상적으로 종료된 런이나 배드런으로 판명되어 종료된 런을 제외하고 나머지 파일들에 대해 에러 여부를 검색한다.
echo "[9/12] Building list of runs that wasn't closed properly."
grep --files-without-match 'The program exited successfully' ${LOGDIR}/* | sed -e "s#${LOGDIR}/##g" -e "s#.log##g" -e "s#./##g" >> tmp_errorlist # grep 의 --files-without-match 는 뒤이어 오는 패턴이 존재하지 않는 파일들의 목록을 리턴한다. 이 경우 'The program exited successfully' 라는 문구가 존재하지 않는 로그파일의 목록을 찾아낸다.
grep 'core dumped|bad run' ${ERRDIR}/* >> ${PWD}/tmp_errorlist
cat tmp_errorlist | sed -e "s#./err/##g" -e "s/.err//g" -e "s/:[0-9]\{6\}//g" -e "s#:/afs/cern.ch/user/w/wyjang/work/KNUTree/exec.sh:##g" > tmp2_errorlist
# ================ End of Working in ${OUTPUTDIR}/$1 ===================================
cd /afs/cern.ch/work/w/wyjang/KNUTree
awk '$1 {print $1}' ${OUTPUTDIR}/$1/tmp2_errorlist > $outfile
echo "[10/12] Removing temporarily generated files"
rm ${OUTPUTDIR}/$1/tmp_errorlist
rm ${OUTPUTDIR}/$1/tmp2_errorlist
echo "Done"
echo "Quarantine is finished! The result is written at <$outfile> file."
echo "[11/12] Now building the new file list."

# 작성된 목록을 바탕으로 파일리스트를 생성한다.
outdir=FileList/rerun_$1
prefix=root://eosams.cern.ch/
datapath=/eos/ams/Data/AMS02/2014/ISS.B950/pass6
surfix="?svcClass=amsuser"
echo "File list directory is $outdir"

if [ ! -d "FileList" ]; then
  mkdir FileList
fi
if [ ! -d "$outdir" ]; then
  mkdir $outdir
else
  rm -rf $outdir
  mkdir $outdir
fi

echo "[12/12] Copying list files into $outdir"
for i in `cat $outfile`
do
  if [ -e "FileList/$1/$i" ]; then
    cp FileList/$1/$i $outdir
  fi
  #echo $i | gawk "BEGIN {i=\"$i\"}{pre=\"$prefix\"}{datapath=\"$datapath/\"}{prefix=\"$prefix\"}{surfix=\"$surfix\"}{printf(\"%s%s%s%s\n\", pre, datapath, i, surfix, \$0)}"
done

echo "Work complete!"
