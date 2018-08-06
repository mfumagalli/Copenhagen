
# https://www.asciitable.com/

#  ~/juliapro/JuliaPro-0.6.2.1/julia

include("/home/mfumagal/Software/ngsJulia/templates.jl");
include("/home/mfumagal/Software/ngsJulia/generics.jl");
include("/home/mfumagal/Software/ngsPoly/functions.jl");

alleles = ['A','C','G','T'];

mySite = Site("chrom", 1, 'A');

myReads = Reads("AAAG", "5555");

calcGenoLogLike1(myReads, mySite);

myReads = Reads("AAAG", "5555");

gl = calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 3);
print("AA: ",gl[1],"\nAG: ",gl[2],"\nGG: ",gl[3])


myReads = Reads("AAAG", "5550");

gl = calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 3);
print("AA: ",gl[1],"\nAG: ",gl[2],"\nGG: ",gl[3])


myReads = Reads("AAAG", "555K");

gl = calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 3);
print("AA: ",gl[1],"\nAG: ",gl[2],"\nGG: ",gl[3])


myReads = Reads("AAAAAAAAAG", "5555555550");

gl = calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 3);
print("AA: ",gl[1],"\nAG: ",gl[2],"\nGG: ",gl[3])




