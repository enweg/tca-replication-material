################################################################################
# Utility functions
################################################################################

lag(x) = vcat(missing, x[1:(end-1)])
diff(x) = x .- lag(x)
