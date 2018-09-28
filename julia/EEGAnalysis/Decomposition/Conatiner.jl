module Container

struct EEGAnalysisContainer
    
    # meta
    subject::String  # subject name or id
    expname::String  # experiment name or id
    
    data_dir::String 
    result_dir::String
    
    # recording
    fs::Int32  # sampling frequency
    iti::Float32
    
    # analysis parameters
    roi::Tuple{Float32,Float32}  # region of interest
    
    # data
    channels::Array{Float32, 2}
    times::Array{Float32, 1}
    markers::Dict{String, Array{Float32, 1}}

end

end