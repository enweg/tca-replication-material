using DataFrames
using TSFrames
using Dates

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
