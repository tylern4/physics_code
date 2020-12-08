using DataFrames;
using CSV;

corr = CSV.read("/Users/tylern/Data/e1d/mom_corr_elastic.dat", DataFrame);

print(corr);