#!/usr/bin/perl

open(F1,"file1.txt");
open(F2,"file2.txt");

@array1=();
@array2=();
$count1=0;
$count2=0;
$sum=0;

while(<F1>)
{
	chomp($_);
	@a=split(/\t/,$_);
	$array1[$count1][0]=$a[0];
	$array1[$count1][1]=$a[1];
	$count1++;
}
close F1; 


while(<F2>)
{
	chomp($_);
	@a=split(/\t/,$_);
	$array2[$count2][0]=$a[0];
	$array2[$count2][1]=$a[1];
	$count2++;
}
close F1;

for($i=0;$i<$count1;$i++)
{
	$diff1=($array1[$i][1]-$array1[$i][0]);
	$diff2=($array2[$i][1]-$array2[$i][0]);
	print abs($diff1-$diff2). "\n";
	$sum+=abs($diff1-$diff2); 	
}
print $sum;
