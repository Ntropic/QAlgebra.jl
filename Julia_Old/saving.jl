using JLD2
using Dates
using CairoMakie

function save_io(func::Function, args...; printing::Bool=false, redo::Bool=false, kwargs...)
    #Saves the inputs and outputs of a function to a file.
    # Generate filename from function name
    # file contains "input" and "output" variables
    # input is a vector of dictionaries containing the input arguments
    # output is a vector of output arguments in arbitrary format
    run_func::Bool = true
    check_input::Bool = false
    filename = joinpath("Results", string(func) * ".jld2")
    if !isdir("Results")
        mkdir("Results")
        run_func = true
    end
    # does the file exist?
    exists::Bool = isfile(filename)
    # if it does, load it using JLD2
    if exists
        # check 
        input = nothing
        jldopen(filename, "r") do file
            # check if input exists in file
            all_keys = keys(file)
            if "input" in all_keys
                # check if args and kwargs exist in file
                check_input = true
                input = file["input"]
            end
        end
    end
    # generate tags for input vector names
    # make a copy of the input arguments in case they get changed by func
    args = deepcopy(args)
    kwargs = deepcopy(kwargs)
    # find method of func that will be dispatched with args (kwargs doesn't pplay a role for this)
    # Use `which` to find the corresponding method
    dispatched_method = which(func, Tuple(typeof.(args)))
    argument_names = Base.method_argnames(dispatched_method)[2:end]
    input_variables::Dict{String,Any} = Dict()
    for (arg_name, arg) in zip(argument_names, args)
        input_variables[string(arg_name)] = arg
    end
    for (key, value) in kwargs
        input_variables[string(key)] = value
    end
    index = -1
    if check_input
        # check in input if input_variables are already there
        # if they are, don't run func (unless redo is true)
        for (i, in_dict) in enumerate(input)
            if in_dict == input_variables
                if !redo
                    run_func = false
                end
                index = i
                break
            end
        end
    end

    if !run_func
        # Load the data
        if printing
            println("Loading data from file $filename")
        end
        jldopen(filename, "r") do file
            results = file["output"][index]
        end
    else
        if printing
            println("Constructing data by running function $func")
        end
        predate = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
        # Run the function
        results = func(args...; kwargs...)
        postdate = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
        # Store inputs in the file
        if exists # file exists
            jldopen(filename, "r+") do file
                if index < 0
                    push!(file["input"], input_variables)
                    push!(file["output"], results)
                    push!(file["predate"], predate)
                    push!(file["postdate"], postdate)
                else
                    file["output"][index] = results
                    file["predate"][index] = predate
                    file["postdate"][index] = postdate
                end
            end
        else
            jldopen(filename, "w") do file
                file["input"] = [input_variables]
                file["output"] = [results]
                file["predate"] = [predate]
                file["postdate"] = [postdate]
            end
        end
    end
    return results
end
#### Test 
#function fun_name(a, b, c=1; d=1)
#    A = a + b + c
#    B = a * b * c+ d
#    return A, B
#end
## run the function using save_io
#A, B = save_io(fun_name, 2, 3, d=3, printing=true)

# Function that saves plots to a file in Plots.jl folder 
function save_plot(fig, filename::String="plot"; redo::Bool=true, folder::String="Plots", kwargs...)
    # Saves the plot to a file.
    # Arguments: 
    # pl: the plot to be saved
    # filename [String]: the name of the file to be saved (default: "plot")
    # redo [Bool]: if true, the file will be saved even if it already exists (default: true)
    # kwargs: additional keyword arguments to be passed to savefig
    #          - e.g. dpi=300, quality=95, etc.
    # Returns: nothing


    # generate Plots directory if missing 
    # add .png if file ending is missing [does not end in .jpg, .png, .pdf, .svg, .pdf]
    if !isdir(folder)
        mkdir(folder)
    end
    if !endswith(filename, ".jpg") && !endswith(filename, ".png") && !endswith(filename, ".pdf") && !endswith(filename, ".svg")
        filename *= ".png"
    end
    filefolder = joinpath(folder, filename)
    if isfile(filefolder) && !redo
        return
    end
    # save the plot
    save(filefolder, fig; kwargs...)
    return fig
end
function save_fig(fig, filename::String="plot"; redo::Bool=true, folder::String="Plots", kwargs...)
    # Saves the plot to a file.
    return save_plot(fig, filename; redo=redo, folder=folder, kwargs...)
end

function save_solution(solution::DifferentialSolution, filename::String="solution"; redo::Bool=true)
    # Saves the solution to a file.
    # Arguments: 
    # solution: the solution to be saved
    # filename [String]: the name of the file to be saved (default: "solution")
    # redo [Bool]: if true, the file will be saved even if it already exists (default: true)
    # Returns: nothing
    # generate Plots directory if missing 
    if !isdir("Results")
        mkdir("Results")
    end
    if !endswith(filename, ".jld2")
        filename *= ".jld2"
    end
    filefolder = joinpath("Results", filename)
    if isfile(filefolder) && !redo
        return
    end
    # save the solution
    jldopen(filefolder, "w") do file
        file["solution"] = solution
    end
    return
end

function load_solution(filename::String="solution")
    # Loads the solution from a file.
    # Arguments: 
    # filename [String]: the name of the file to be loaded (default: "solution")
    # Returns: the solution
    # generate Plots directory if missing 
    if !isdir("Results")
        mkdir("Results")
    end
    if !endswith(filename, ".jld2")
        filename *= ".jld2"
    end
    filefolder = joinpath("Results", filename)
    
    # load the solution
    solution = jldopen(filefolder, "r") do file
        return file["solution"]
    end
    return solution
end