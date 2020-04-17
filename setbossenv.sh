#!/bin/bash

ver=$1
HOMEPATH="/workfs/bes/$USER"
ENVIRONMENT="cmthome-$ver"

# starting ...

if [ $ver = "" ]
then
	exit -1
fi

if [ -d $HOMEPATH/$ENVIRONMENT ] && [ -d $HOMEPATH/workarea-$ver ]
then
	echo "You will delete BOSS-$ver, do you want to rebuidt? (y/n)"
	read var
	if [ $var == "n" ] || [ $var == "no" ]
	then
		exit -1
	fi
fi

pushd $HOMEPATH

if [ -d cmthome-$ver ]
then
	rm -rf cmthome-$ver
fi
mkdir cmthome-$ver

if [ -d workarea-$ver ]
then
	rm -rf workarea-$ver
fi
mkdir workarea-$ver

popd

pushd $HOMEPATH/cmthome-$ver

cp -r /afs/ihep.ac.cn/bes3/offline/Boss/cmthome/$ENVIRONMENT/*  . 

#
sed -e "s/^#path/path/" -e "s/^#macro WorkArea/macro WorkArea/g" requirements >tmp_requirements
sed -e "s#/ihepbatch/bes/maqm/workarea#$HOMEPATH/workarea-$ver#g" tmp_requirements > requirements
rm -f tmp_requirements

#
if [ -e setupCVS.sh~ ]
then
   rm -f setupCVS.sh~
fi
sed -e "s/maqm/$USER/" setupCVS.sh |sed -e "s/\$user/$USER/g" > setupCVS.sh~
#
mv -f setupCVS.sh~  setupCVS.sh
#
source setupCMT.sh
cmt config
source setupCVS.sh
source setup.sh

#TestRelease
echo "TestRelease ... ... "
cp -r $BesArea/TestRelease $HOMEPATH/workarea-$ver

popd

pushd $HOMEPATH/workarea-$ver/TestRelease/*/cmt
rm -rf ../x86*
cmt broadcast cmt config
source setup.sh 
cmt broadcast make

echo "boss environment is OK!"

popd

echo -e "source /workfs/bes/$USER/cmthome-$ver/setupCVS.sh\nsource /workfs/bes/$USER/cmthome-$ver/setupCMT.sh\nsource /workfs/bes/$USER/cmthome-$ver/setup.sh\nsource /workfs/bes/$USER/workarea-$ver/TestRelease/*/cmt/setup.sh\n\nalias wk=\"cd /workfs/bes/$USER/workarea-$ver\"" > bossenv$ver

