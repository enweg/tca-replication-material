using DataFrames
using TSFrames
using Dates


"""
    _replace_NaN_string(df::DataFrame)
    _replace_NaN_string(vec::AbstractVector)

Replace "NaN" strings with `NaN64`. 
"""
function _replace_NaN_string(vec::AbstractVector)
    return map(x -> x=="NaN" ? NaN64 : x, vec)
end
function _replace_NaN_string(df::DataFrame)
    for i in axes(df, 2)
        df[!, i] = _replace_NaN_string(df[!, i])
    end
    return df
end

"""

Transform a quarter string to a valid date. 

## Arguments

- `string::String`: Some string containing year quarter information
- `pattern::Regex`: Must be a regular expression that captures two groups. The
  first must be the year and the second must be the quarter

## Keyword Arguments

- `end_of_quarter::Bool=false`: Should the last day of the quarter be used as the
  date. If false, the first day of the quarter will be used.

## Example

- transform the string `"2020Q1"` into a valid date

```
string = "2020Q1"
date = _quarter_string_to_date(string)
```

- transform the string `"2020-Q1"` into a valid date
```
string = "2020-Q1"
pattern = r"([0-9]{4})-Q([1-4]{1})"
date = _quarter_string_to_date(string, pattern)
```

"""
function _quarter_string_to_date(string::String, pattern::Regex=r"([0-9]{4})Q([1-4]{1})"; end_of_quarter::Bool=false)
    m = match(pattern, string)
    year, quarter = parse.(Int, m.captures)
    month = quarter*3
    date = Dates.Date("$year/$month", dateformat"yyyy/mm")
    # date = Dates.lastdayofmonth(date)
    date = end_of_quarter ? Dates.lastdayofquarter(date) : Dates.firstdayofquarter(date)
    # if !end_of_quarter
    #     date = date - Dates.Quarter(1) + Day(1)
    # end
    return date
end

"""
    _drop_nan(ts::TSFrame)
    _drop_nan(df::DataFrame)

Drop NaN values from a DataFrame or a TSFrame. 
"""
function _drop_nan(df::DataFrame)
    return filter(row -> all(x -> !(x isa Number && isnan(x)), row), df)
end
function _drop_nan(ts::TSFrame)
    return TSFrame(_drop_nan(ts.coredata))
end

function _drop_missing(df::DataFrame)
    return filter(row -> all(x -> !ismissing(x), row), df)
end
function _drop_missing(ts::TSFrame)
    return TSFrame(_drop_missing(ts.coredata))
end

function _lag(vec::AbstractVector, l::Int=1, pad=missing)
    return vcat(fill(pad, l), vec[begin:end-l])
end

function has_missing_in_between(vec::AbstractVector)
    where_missing = findall(x -> ismissing(x), vec)
    length(where_missing) == 0 && return false
    return !all((where_missing - _lag(where_missing))[2:end] .== 1)
end

function DataFrames.select(ts::TSFrame, args...; kwargs...)
    return TSFrame(DataFrames.select(ts.coredata, :Index, args...; kwargs...))
end

"""
Replace all `missing` entries in a vector with zero.
"""
function _replace_missing_zero(x::AbstractVector{Union{Missing, T}}) where {T<:Real}
    x[ismissing.(x)] .= zero(T)
    return x
end
