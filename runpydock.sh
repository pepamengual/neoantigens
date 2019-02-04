pydock='/home/pepamengual/pydock/pyDock3/pyDock3'
folder="all"
cd $folder
for i in *.pdb;
do
filename=$(echo $i | rev | cut -c 5- | rev)
### .ini file
printf '[receptor] 
pdb     = XXX
mol     = YYY
newmol  = YYY

[ligand]
pdb     = XXX
mol     = ZZZ
newmol  = ZZZ' > $filename.ini
###

sed -i "s/XXX/$i/g" $filename.ini
firstchain=$(grep "^ATOM" $i | cut -c22 | sort -u | head -n1)
secondchain=$(grep "^ATOM" $i | cut -c22 | sort -u | tail -n1)
sed -i "s/YYY/$firstchain/g" $filename.ini
sed -i "s/ZZZ/$secondchain/g" $filename.ini


$pydock $filename bindEy

done
