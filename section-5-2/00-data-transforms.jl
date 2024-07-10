using DrWatson
@quickactivate "gov-cons-anticipation"

using DataFrames, CSV
using DataFramesMeta
using Dates

function date_to_quarter(data::Date)
  y = year(data)
  q = quarterofyear(data)
  return y + (q-1)/4
end

include(scriptsdir("utilities.jl"))

################################################################################
# Quarterly Defense Spending
# Note that the spending is seasonally adjusted and annualised
# Measured in billions
################################################################################

path_fred_quarterly = datadir("raw", "Defense_quarterly.csv")
fred_quarterly = DataFrame(CSV.File(path_fred_quarterly))
first(fred_quarterly, 5)

fred_quarterly = @chain fred_quarterly begin
  @transform begin
    :quarter = date_to_quarter.(:observation_date)
    :year = year.(:observation_date)
  end
  @rename :govdef_q = :A997RC1Q027SBEA_20240125
  @rtransform :govdef_q = :govdef_q / 4  # de-annualising
  @select :quarter :year :govdef_q
end
first(fred_quarterly, 5)

################################################################################
# Annual Defense Spending 
# Measured in billions
################################################################################

path_fred_annual = datadir("raw", "Defense_annual.csv")
fred_annual = DataFrame(CSV.File(path_fred_annual))
first(fred_annual, 5)

fred_annual = @chain fred_annual begin
  @transform :year = year.(:observation_date)
  @rename :govdef_a = :A997RC1A027NBEA_20240125 
  @select :year :govdef_a
end
first(fred_annual, 5)

################################################################################
# Merging the FRED data
################################################################################

fred = leftjoin(fred_quarterly, fred_annual, on = :year)
fred = @chain fred begin
  groupby(:year)
  # checking whether quarterly data makes sense
  @transform :govdef_a_sum = sum(:govdef_q)  
  # how much of yearly data is attributed to each quarter
  @transform :share_quarter = :govdef_q ./ :govdef_a  
end
first(fred, 5)

################################################################################
# Military Expenditure data from Our World in Data
# Data is as percentage of GDP
################################################################################

path_defense_owid = datadir("raw", "military-expenditure-as-a-share-of-gdp-long.csv")
defense_owid = DataFrame(CSV.File(path_defense_owid))
# renaming the long column name
names(defense_owid)
col_name = names(defense_owid)[end]
rename!(defense_owid, col_name => :mil_exp_perc_gdp)

defense_owid = @chain defense_owid begin
  @rsubset :Code == "USA"
  @rename :year = :Year
  @select :year :mil_exp_perc_gdp
end
first(defense_owid, 5)

################################################################################
# Ramey 2018 data
################################################################################

path_ramey = datadir("raw", "RZDAT.csv")
ramey = DataFrame(CSV.File(path_ramey))
ramey = @chain ramey begin
  @select :quarter :ngov :ngdp 
  @rtransform :year = floor(Int, :quarter)
end
first(ramey, 5)

ramey_annual = @chain ramey begin
  groupby(:year)
  @combine begin
    # GDP is annualised rate. I assume the same is true for Gov spending
    :ngov = sum(:ngov ./ 4)
    :ngdp = sum(:ngdp ./ 4)
  end
end

################################################################################
# Merging Ramey and OWID
################################################################################

mil_spending = leftjoin(defense_owid, ramey_annual, on=:year)
mil_spending = @chain mil_spending begin
  @orderby(:year)
  @rsubset :year >= 1890
  @rtransform :mil_exp = :ngdp * :mil_exp_perc_gdp / 100
end

################################################################################
# Merging with FRED data
################################################################################

mil_spending = leftjoin(mil_spending, fred_annual; on=:year)
mil_spending = @orderby(mil_spending, :year)

################################################################################
# Plot showing that the obtained defense spending from OWID largely 
# agrees with the official FRED data on a yearly frequency. 
################################################################################

using CairoMakie
CairoMakie.activate!()

# figure in levels
fig = Figure();
ax = Axis(fig[1, 1]; title="FRED and OWID data largely agree on a yearly frequency.");
lines!(ax, mil_spending.year, mil_spending.mil_exp; color=:red, label="OWID");
lines!(ax, mil_spending.year, mil_spending.govdef_a; color=:black, label="FRED");
Legend(fig[2, :], ax; orientation=:horizontal, framevisible=false);
save(plotsdir("FRED-OWID-comparison.pdf"), fig);

function diff(x::AbstractArray)
  return vcat(missing, x[2:end] - x[1:(end-1)])
end

# figure in diff
fig = Figure();
ax = Axis(fig[1, 1]; title="FRED and OWID data largely agree on a yearly frequency.");
lines!(ax, mil_spending.year, diff(mil_spending.mil_exp); color=:red, label="OWID");
lines!(ax, mil_spending.year, diff(mil_spending.govdef_a); color=:black, label="FRED");
Legend(fig[2, :], ax; orientation=:horizontal, framevisible=false);
save(plotsdir("FRED-OWID-comparison-diff.pdf"), fig);

################################################################################
# Getting quarterly historical data
################################################################################

# first getting share of yearly government spending attributed to each quarter
ramey_quarter_share = @chain ramey begin
  @rename :ngov_q = :ngov
  @rename :ngdp_q = :ngdp
  leftjoin(ramey_annual; on=:year)
  @rtransform :quarter_share = :ngov_q / 4 / :ngov
  @orderby :quarter
end

# calculating quarterly military spending by applying the quarterly 
# government spending shares to yearly military spending
spending_quarter = @chain ramey_quarter_share begin
  leftjoin(select(mil_spending, :year, :mil_exp); on=:year)
  @rtransform :mil_spending_q = :quarter_share * :mil_exp 
  @select :quarter :ngov_q :ngdp_q :mil_spending_q
  @rsubset :quarter >= 1890
  @orderby :quarter
end

# Comparing this to official fred data visiually
# The two are a very close match
tmp = leftjoin(spending_quarter, fred_quarterly; on=:quarter)
tmp = @orderby(tmp, :quarter)
fig = Figure();
ax = Axis(fig[1, 1]);
lines!(ax, tmp.quarter, tmp.mil_spending_q; color=:red);
lines!(ax, tmp.quarter, tmp.govdef_q; color=:black);
save(plotsdir("FRED-OWID-comparison-quarterly.pdf"), fig);

# comparing in differences
# The match is worse but still reasonable.
fig = Figure();
ax = Axis(fig[1, 1]);
lines!(ax, tmp.quarter, diff(tmp.mil_spending_q); color=:red, label="OWID");
lines!(ax, tmp.quarter, diff(tmp.govdef_q); color=:black, label="FRED");
Legend(fig[2, :], ax; orientation=:horizontal, framevisible=false);

# The quarterly military spending series has a correlation above 0.99
# with the official FRED quarterly spending series. Correlation in differences
# is only 0.40.
using Statistics
@chain tmp begin
  @select :mil_spending_q :govdef_q 
  dropmissing
  Matrix
  cor
end
@chain tmp begin
  @select :mil_spending_q :govdef_q 
  @transform :mil_spending_q = diff(:mil_spending_q)
  @transform :govdef_q = diff(:govdef_q)
  dropmissing
  Matrix
  cor
end

# how much are the shares correlated with each other? 
# Highly with 0.82 correlation
share = @chain ramey_quarter_share begin
  @select :quarter :quarter_share 
  @rename :share_quarter_ramey = :quarter_share 
  leftjoin(select(fred, :quarter, :share_quarter => :share_quarter_fred); on=:quarter)
  @orderby :quarter
end
@chain share begin
  @select :share_quarter_ramey :share_quarter_fred
  dropmissing
  Matrix
  cor
end

################################################################################
# Saving
################################################################################

spending_quarter = select(spending_quarter, :quarter, :mil_spending_q => :govdef_owid)
first(spending_quarter, 5)
spending_quarter = @chain spending_quarter begin
  leftjoin(fred_quarterly; on=:quarter)
  @rename :govdef_fred = :govdef_q
  @rtransform :govdef_spliced = ismissing(:govdef_fred) ? :govdef_owid : :govdef_fred
  @select :quarter :govdef_owid :govdef_fred :govdef_spliced
  @orderby :quarter
end

CSV.write(datadir("spending_quarter.csv"), spending_quarter)

################################################################################
# Creating dataset for LP estimation
################################################################################

ramey2018 = DataFrame(CSV.File(datadir("raw", "RZDAT.csv")))
govdef = DataFrame(CSV.File(datadir("spending_quarter.csv")))

data = outerjoin(govdef, ramey2018, on=:quarter)
data = select(data, :quarter, :)
data = @orderby(data, :quarter)

data = @chain data begin
    @select(
        :quarter,
        :pop,  # population
        :ynorm = :rgdp_pott6,  # real potential GDP, based on 6th degree polynomial fit from 1889:1 - 2015:4, ommiting Great Depression and WWII
        :news,  # military news instrument
        :pgdp,  # GDP implicit price deflator 
        :ngov,  # nominal government purchases
        :ngdp,  # nominal GDP
        :rgdp,  # real GDP
        :nfedcurrreceipts_nipa,  # nominal Federal current receipts, NIPA accrual basis 
        :pubfeddebt_treas,  # nominal federal debt in the hands of the public, cash basis
        :govdef_owid,
        :govdef_spliced,
        :tbill
    )
    @transform :newsy = :news ./ (lag(:ynorm) .* lag(:pgdp))
    @transform :rgov = :ngov ./ :pgdp
    @transform :rgovdef_owid = :govdef_owid ./ :pgdp
    @transform :rgovdef_spliced = :govdef_spliced ./ :pgdp
    @transform :rtax = :nfedcurrreceipts_nipa ./ :pgdp
    @transform :taxy = :nfedcurrreceipts_nipa ./ :ngdp
    @transform :debty = :pubfeddebt_treas ./ lag(:ngdp)
    @transform :lpgdp = log.(:pgdp)
    @transform :ly = log.(:rgdp)
    @transform :infl = 400 .* diff(:lpgdp)
    @transform :y = :rgdp ./ :ynorm
    @transform :g = :rgov ./ :ynorm
    @transform :gdef_owid = :rgovdef_owid ./ :ynorm
    @transform :gdef_spliced = :rgovdef_spliced ./ :ynorm

end

data_baseline = select(data, :quarter, :newsy, :y, :g, :gdef_spliced)
data_baseline = @rsubset(data_baseline, :quarter <= 2015 && :quarter >= 1890)
data_baseline = @chain data_baseline begin
    @rtransform begin
        :newsy = Float64(:newsy)
        :y = Float64(:y)
        :g = Float64(:g)
        :gdef = Float64(:gdef_spliced)
    end
    @select :gdef :g :y :newsy :quarter
end

data_lp = select(data_baseline, Not(:quarter))

CSV.write(datadir("data_lp.csv"), data_lp)
