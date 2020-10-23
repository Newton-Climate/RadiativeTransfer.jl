#####
##### Types for describing an Instrument Operations (convolution, resampling)
#####

"""
    type InstrumentOperator
Abstract InstrumentOperator type 
"""
abstract type InstrumentOperator end



"""
    struct KernelInstrument{FT}

A struct which provides all parameters for the convolution with a kernel

# Fields
$(DocStringExtensions.FIELDS)
"""
struct KernelInstrument{FT} <: InstrumentOperator
    "convolution Kernel" 
    kernel::FT
    "Output spectra grid"
    ν_out::Array
end;