module EDFReader

using MAT
using JLD

struct EDFData
    
    # Header
    version::Int64                      # 8 ascii : version of this data format (0) 
    patient_info::String                # 80 ascii : local patient identification (mind item 3 of the additional EDF+ specs)
    record_info::String                 # 80 ascii : local recording identification (mind item 4 of the additional EDF+ specs)
    start_date::String                  # 8 ascii : startdate of recording (dd.mm.yy) (mind item 2 of the additional EDF+ specs)
    start_time::String                  # 8 ascii : starttime of recording (hh.mm.ss) 
    header_length::Int64                # 8 ascii : number of bytes in header record 
                                        # 44 ascii : reserved 
    recordnum::Int64                    # 8 ascii : number of data records (-1 if unknown, obey item 10 of the additional EDF+ specs) 
    sampleduration::Float32             # 8 ascii : duration of a data record, in seconds 
    nchannel::Int64                     # 4 ascii : number of signals (ns) in data record 
    channelLabels::Array{String, 1}     # ns * 16 ascii : ns * label (e.g. EEG Fpz-Cz or Body temp) (mind item 9 of the additional EDF+ specs)
    channelType::Array{String, 1}       # ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode) 
    physical_dim::Array{String, 1}      # ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC) 
    physical_min::Array{Float32, 2}     # ns * 8 ascii : ns * physical minimum (e.g. -500 or 34) 
    physical_max::Array{Float32, 2}     # ns * 8 ascii : ns * physical maximum (e.g. 500 or 40) 
    digital_min::Array{Int16, 2}        # ns * 8 ascii : ns * digital minimum (e.g. -2048) 
    digital_max::Array{Int16, 2}        # ns * 8 ascii : ns * digital maximum (e.g. 2047) 
    prefiltering::Array{String, 1}      # ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz) 
    samples::Array{Int16, 1}            # ns * 8 ascii : ns * nr of samples in each data record 
    reserved_samples::Array{Int16, 1}   # ns * 32 ascii : ns * reserved
    
    # Data
    data::Array{Float32, 2}
    reserved_data::Array{Float32, 2}
end

function readEDFFile(filename::String, analog::Bool=false)
    rawfile = open(filename)
    
    #Header
    version = read(rawfile, UInt8, 8)       .|> Char |> String |> parse
    patient_info = read(rawfile, UInt8, 80) .|> Char |> String |> strip
    record_info = read(rawfile, UInt8, 80)  .|> Char |> String |> strip
    start_date = read(rawfile, UInt8, 8)    .|> Char |> String |> strip
    start_time = read(rawfile, UInt8, 8)    .|> Char |> String |> strip
    header_length = read(rawfile, UInt8, 8) .|> Char |> String |> parse
    reserved = read(rawfile, UInt8, 44)
    recordnum = read(rawfile, UInt8, 8)     .|> Char |> String |> parse
    sampleduration = read(rawfile, UInt8, 8).|> Char |> String |> parse
    nchannel = read(rawfile, UInt8, 4)      .|> Char |> String |> parse
    
    channelLabels = Array{String,1}(nchannel)
    for chidx = 1:nchannel
        channelLabels[chidx] = read(rawfile, UInt8, 16) .|> Char |> String |> strip
    end

    channelType = Array{String,1}(nchannel)
    for chidx = 1:nchannel
        channelType[chidx] = read(rawfile, UInt8, 80) .|> Char |> String |> strip
    end
    
    physical_dim = Array{String,1}(nchannel)
    for chidx = 1:nchannel
        physical_dim[chidx] = read(rawfile, UInt8, 8) .|> Char |> String |> strip
    end
    
    physical_min = zeros(Float32, (1, nchannel))
    for chidx = 1:nchannel
        physical_min[chidx] = read(rawfile, UInt8, 8) .|> Char |> String |> parse
    end
    
    physical_max = zeros(Float32, (1, nchannel))
    for chidx = 1:nchannel
        physical_max[chidx] = read(rawfile, UInt8, 8) .|> Char |> String |> parse
    end
    
    digital_min = zeros(Int16, (1, nchannel))
    for chidx = 1:nchannel
        digital_min[chidx] = read(rawfile, UInt8, 8) .|> Char |> String |> parse
    end
    
    digital_max = zeros(Int16, (1, nchannel))
    for chidx = 1:nchannel
        digital_max[chidx] = read(rawfile, UInt8, 8) .|> Char |> String |> parse
    end
    
    prefiltering = Array{String, 1}(nchannel)
    for chidx = 1:nchannel
        prefiltering[chidx] = read(rawfile, UInt8, 80) .|> Char |> String |> strip
    end
    
    samples = Array{Int16, 1}(nchannel)
    for chidx = 1:nchannel
        samples[chidx] = read(rawfile, UInt8, 8) .|> Char |> String |> parse
    end
    
    reserved_samples = Array{Int16, 1}(nchannel)
    for chidx = 1:nchannel
        reserved_samples[chidx] = read(rawfile, UInt8, 32) .|> Char |> String |> parse |> (x)-> x==nothing?0:x
    end
    
    # reserved = read(rawfile, UInt8, header_length-256+nchannel*256)
    
    #Data
    data = zeros(Float32, (recordnum*samples[1], nchannel))
    reserved_data = zeros(Float32, (recordnum*reserved_samples[1], nchannel))
    
    step = samples |> sum
    for ri = 1:recordnum
        record_data = read(rawfile, Int16, step)
        data[samples[1]*(ri-1)+1:samples[1]*ri,:] = reshape(record_data, (samples[1], nchannel))
    end
    
    close(rawfile)
    
    if analog
        data = (data .- digital_min) ./ (digital_max - digital_min) .* 
               (physical_max - physical_min) .+ physical_min
    end
    
    return EDFData( version, patient_info, record_info,
                    start_date, start_time, header_length,
                    recordnum, sampleduration, nchannel,
                    channelLabels, channelType, physical_dim, 
                    physical_min, physical_max, digital_min, digital_max,
                    prefiltering, samples, reserved_samples,
                    data, reserved_data)
end


function edf2mat(filename::String, edfdata::EDFData)
    matdict = Dict([ String(key) => getfield(edfdata, key) for key in fieldnames(edfdata)])
    matwrite(filename, matdict)
end


function edf2jld(filename::String, edfdata::EDFData)
    matdict = Dict([ String(key) => getfield(edfdata, key) for key in fieldnames(edfdata)])
    JLD.save(filename, matdict)
end

end