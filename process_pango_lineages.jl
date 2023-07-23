import FASTX
import BioSequences
using  GenomicAnnotations
import BioGenerics
import JSON
using ArgParse
using DataFrames
using Combinatorics
using DataStructures
using CSV
using DelimitedFiles
using Dates
ENV["COLUMNS"]=2000

function get_closest_ancestor(lineage_alias,pseudoreference_dict)
    parents= []
    for i in range(1,stop=length(split(lineage_alias,".")),step=1)
        new_parent = join(split(lineage_alias,".")[1:i],".")
        append!(parents,[new_parent])
    end
    parents = reverse(parents)[2:end]
    for parent in parents
        if haskey(pseudoreference_dict,parent)
            return parent
        end
    end
    return "Root"
end;

function convert_lineage_name(lineage,lineage_aliases)

    if haskey(lineage_aliases, split(lineage,".")[1]) 
        alias = [lineage_aliases[split(lineage,".")[1]]]
        remainder = split(lineage,".")[2:end]
        return join(vcat(alias,remainder) ,".")
    elseif length(split(lineage,".")) >4 && haskey(lineage_aliases,join(split(lineage,".")[1:5]))
        alias = [lineage_aliases[split(lineage,".")[1:5]]]
        remainder = split(lineage,".")[5:end]
        return join(vcat(alias,remainder) ,".")
    end
    return lineage
end;

function is_recombinant_lineage(lineage,lineage_aliases)
    if haskey(lineage_aliases, split(lineage,".")[1]) 
        alias = lineage_aliases[split(lineage,".")[1]]
        println(alias isa String)
        if alias isa String == false
            return true
        end
    end
    return false
end;

function make_pseudoreference_dict(fasta)
    pseudoreference_dict = Dict()
    for record in fasta
        pseudoreference_dict[BioGenerics.seqname(record)] = BioGenerics.sequence(record)
    end
    return pseudoreference_dict
end;

function process_recombinant(recombinant,lineage_aliases)
    if occursin(".",recombinant)
        return [],[]
    end
    recomb_list = Vector{String}(lineage_aliases[recombinant])
    processed_recomb_list = []
    for r in recomb_list
        r = replace(r,"*" => "" )*"."*recombinant
        if occursin("/",r)
            split_r = split(r,".")
            start_of_r = split_r[1:end-1]
            end_of_r = split(split_r[end],"/" )
            for r_end in end_of_r
                println(join(start_of_r,".")*"."*r_end)
                push!(processed_recomb_list,join(start_of_r,".")*"."*r_end)
            end
        else
            push!(processed_recomb_list,r)
        end
    end
    alias_column = [recombinant for i in 1:length(processed_recomb_list)]
    return processed_recomb_list,alias_column
end

function dealiase(lineage_list,lineage_aliases)
    dealiased = false
    lineages_aliased = lineage_list
    while dealiased == false
        old_lineages_aliased = lineages_aliased
        lineages_aliased = [string(convert_lineage_name(lineage,lineage_aliases)) for lineage in old_lineages_aliased]
        if isequal(old_lineages_aliased, lineages_aliased)
            dealiased=true
        end
    end
    return lineages_aliased
end


function add_node_to_tree(node_row,tree)
    found_node = false
    #Get Root Node
    current_node = collect(keys(tree))[1]
    while found_node == false
        node_to_place = node_row.Lineage
        node_to_place_split = split(node_to_place,".")
        for child in current_node["children"] 
            if child.full_name == node_row.parent
                new_node =  Dict{String,Any}("name"=>node_row.Alias,"full_name"=>node_row.Lineage,"parent"=>node_row.Parent,"alias_parent"=>node_row.Parent_Alias,"children"=>[])
                push!(child.children,node_row)
                found_node = true
                break
            end
        end
    end
end

function alias_lineages(lineage,lineage_aliases)
    if count(".",lineage) <= 3 
        return lineage
    else
        alias_level = trunc(Int,count(".",lineage)/3)*3+1
        full_name = join(split(lineage,".")[1:alias_level],".")
        remainder = split(lineage,".")[alias_level+1:end]
        if length(remainder) == 0
            full_name =  join(split(lineage,".")[1:alias_level-3],".")
            remainder = split(lineage,".")[(alias_level-3)+1:end]
        end
    end
    
    for (key,value) in lineage_aliases
        if value isa String
            if full_name == value
                alias = key
#                 println(join(vcat(alias,remainder) ,"."))
                return lineage =  join(vcat(alias,remainder) ,".")
            end
        end
    end
end



function build_tree(node_dataframe,extend=Dict{String,Dict{String,Any}}(),annotations=DataFrame([]) )
    
    #Look for recombinant empty nodes
    recombinant_subtrees = [extend[key] for key in collect(keys(extend)) if extend[key]["parent"] == "Recombinant"]
    
    tree_dict = extend
    for lineage_row in eachrow(node_dataframe)
        
        estimated_lineage_count = 0
        if nrow(annotations) != 0
            lineage_annotation = filter(:Lineage => n -> n == lineage_row["Alias"], annotations)
            if nrow(lineage_annotation) != 0
                estimated_lineage_count = lineage_annotation[1,"# assigned"]
            end
        end

        #If new node create empty node
        if haskey(tree_dict,lineage_row.Lineage) == false
            tree_dict[lineage_row.Lineage] = Dict{String,Any}("name"=>lineage_row.Alias,"full_name"=>lineage_row.Lineage,"parent"=>lineage_row.Parent,"alias_parent"=>lineage_row.Parent_Alias ,"children"=>[],"estimated_count"=>estimated_lineage_count)
            #Add recombinant subtrees
            if occursin("X",lineage_row.Alias) == true
                for node in recombinant_subtrees
                    if node["name"] == lineage_row.Alias
                        tree_dict[lineage_row.Lineage]["children"] = node["children"]
                        continue
                    end
                end     
            end
        end
        height = lineage_row.Level
        child_height = lineage_row.Level+1
        possible_children = filter(node_dataframe -> node_dataframe.Level == child_height, node_dataframe)[!,:Lineage]

        for child_node in possible_children
            child_row_df = filter(node_dataframe -> node_dataframe.Lineage == child_node, node_dataframe)[1,:]
            if occursin(string(lineage_row.Lineage)*".",string(child_row_df.Lineage)) 
                push!(tree_dict[lineage_row.Lineage]["children"],tree_dict[child_row_df.Lineage])
            end     
        end
        
    end
    return tree_dict
end
# cd("/home4/2191618l/Github/PangoTreeFiles/")

cog_metadata_url = "https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz"
download(cog_metadata_url,"/home4/2191618l/Github/PangoTreeFiles/cog_metadata.csv.gz")


cog_uk = DataFrame(CSV.File("/home4/2191618l/Github/Pango_Tree/cog_metadata.csv",delim=","))
println("Lineage Meta Counts Loaded")

count_dataframe = combine(groupby(cog_uk, [:usher_lineage]), nrow => :count)
rename!(count_dataframe,["Lineage","# assigned"])
sort!(count_dataframe,"# assigned",rev=true)

max_count = filter( r -> r["# assigned"] ==  maximum(count_dataframe[!,"# assigned"]),count_dataframe  )[!,2]

pango_lineage_notes_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
pango_aliases_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"

download(pango_lineage_notes_url,"/home4/2191618l/Github/PangoTreeFiles/local_lineage_notes.tsv")
download(pango_aliases_url,"/home4/2191618l/Github/PangoTreeFiles/local_aliases.json")

lineages = DataFrame(CSV.File("/home4/2191618l/Github/PangoTreeFiles/local_lineage_notes.tsv",delim="\t"))[!,1]
lineage_aliases = JSON.parse(open("/home4/2191618l/Github/PangoTreeFiles/local_aliases.json"))
filter!(p->(p[2] !=""),lineage_aliases)

recombinant_lineages = [recombinant for recombinant in lineages if occursin("X", recombinant) == true]
recombinant_sub_lineages = [recombinant for recombinant in recombinant_lineages if occursin(".", recombinant) == true]
lineages = [string(lineage) for lineage in lineages if occursin("*", lineage) == false && occursin("X", lineage) == false]
println("Separate Recombinants out")

recombinant_aliased_column = []
recombinant_lineage_column = []
for recombinant in recombinant_lineages
    recombinant_lineage,recombinant_alias = process_recombinant(recombinant,lineage_aliases)
    if length(recombinant_lineage) == 0
        continue
    else
        append!(recombinant_aliased_column,recombinant_lineage)
        append!(recombinant_lineage_column,recombinant_alias)
    end
end
recombinant_aliased_column = dealiase(recombinant_aliased_column,lineage_aliases)
recombinant_ordered_lineages = sort!(DataFrame(Lineage=recombinant_aliased_column,Alias=recombinant_lineage_column ,Level = [count(i->(i=='.'), lineage) for lineage in recombinant_aliased_column]),:Level)
recombinant_ordered_lineages[!,:Parent]= [occursin(".",lineage) ? join(split(lineage,".")[1:end-1],".") : "Root"  for lineage in recombinant_ordered_lineages[!,:Lineage]]
recombinant_ordered_lineages[!,:Parent_Alias] = [alias_lineages(l,lineage_aliases) for l in recombinant_ordered_lineages[!,:Parent] ]
println("Process Recombinants")
recombinant_ordered_lineages


recombinant_sub_ordered_lineages = sort!(DataFrame(Lineage=recombinant_sub_lineages,Alias=recombinant_sub_lineages ,Level = [count(i->(i=='.'), lineage) for lineage in recombinant_sub_lineages]),:Level)
recombinant_sub_ordered_lineages[!,:Parent]= [occursin(".",lineage) ? join(split(lineage,".")[1:end-1],".") : "Recombinant"  for lineage in recombinant_sub_ordered_lineages[!,:Lineage]]
sort!(recombinant_sub_ordered_lineages,:Level,rev=true)
recombinant_sub_ordered_lineages = vcat(recombinant_sub_ordered_lineages,DataFrame(Lineage=recombinant_lineages, Alias=recombinant_lineages, Level = 0, Parent="Recombinant"))
recombinant_sub_ordered_lineages[!,:Parent_Alias] = recombinant_sub_ordered_lineages[!,:Parent]
println("Processed recombinant subtrees")


lineages_aliased = dealiase(lineages,lineage_aliases)
println("Lineages Aliased")
ordered_lineages = sort!(DataFrame(Lineage = lineages_aliased,Alias=lineages,Level = [count(i->(i=='.'), lineage) for lineage in lineages_aliased]),:Level)
ordered_lineages[!,:Parent]= [occursin(".",lineage) ? join(split(lineage,".")[1:end-1],".") : "Root"  for lineage in ordered_lineages[!,:Lineage]]
ordered_lineages[!,:Parent_Alias]= [alias_lineages(l,lineage_aliases) for l in ordered_lineages[!,:Parent] ]
println("Lineages in Order")

ordered_lineages = vcat(ordered_lineages,recombinant_ordered_lineages)
sort!(ordered_lineages,:Level,rev=true)


#Add Non-Recombinants to Tree
node_dict = build_tree(ordered_lineages,build_tree(recombinant_sub_ordered_lineages),count_dataframe)

output =Dict{String,Any}("name"=>"Root","alias"=>"Root","children"=>[node_dict["A"],node_dict["B"]])

# Add Recombinants to Node Dict

using JSON
open("/home4/2191618l/Github/PangoTreeFiles/julia_lineages.json","w") do f
    JSON.print(f,output)
end

node_dict_list = collect(values(node_dict))
using JSON
open("/home4/2191618l/Github/PangoTreeFiles/node_list.json","w") do f
    JSON.print(f,node_dict_list)
end

# Make Recombinant Graphs
recombinant_parents = Dict([])
for recombinant_row in eachrow(recombinant_ordered_lineages)
    if haskey(recombinant_parents,recombinant_row.Parent_Alias) == false
        recombinant_parents[recombinant_row.Parent_Alias] = [recombinant_row.Alias]
    else
        push!(recombinant_parents[recombinant_row.Parent_Alias],recombinant_row.Alias)
    end
end

recombinant_children = Dict([])

for node in node_dict_list
#     println(node)
#     println(node["name"] in recombinant_lineages)
    if node["name"] in recombinant_lineages
        if (occursin(".",node["name"]))
            continue
        end
        for kid in node["children"]
            kid = kid["name"]
            if haskey(recombinant_children,node["name"]) == false
                recombinant_children[node["name"]] = Set([kid])
            else
                push!(recombinant_children[node["name"]],kid)
            end
        end
    end
    
end

println(recombinant_parents)

graph = Dict([("nodes",[]),("links",[])])
for (key,values) in merge(recombinant_parents)
    println(key)
    if key === nothing
        println("Key labelled as nothing")
        println(values)
    else
        if occursin("X",key) && occursin(".",key) ==false
            println(key)
            push!(graph["nodes"],Dict([("id",key),("colour","orange")]))
            for val in values
            push!(graph["links"], Dict([("source",value),("target",key)])) 
            end
        elseif occursin(".",key)
            push!(graph["nodes"],Dict([("id",key),("colour","red")]))
            for val in values
            push!(graph["links"], Dict([("source",value),("target",key)])) 
            end
        end
    end
end

for (key,values) in merge(recombinant_children)
    if key === nothing
        println("Key labelled as nothing")
        println(values)
    else
        for child in values
            push!(graph["links"], Dict([("source",key),("target",child)]))
            node_for_child=false
            for node in graph["nodes"]
                if child== node["id"]
                    node_for_child == true
                end
            end
            if node_for_child == false
                if occursin(".",child)
                    push!(graph["nodes"],Dict([("id",child),("colour","blue")]))
                elseif occursin("X",child) && occursin(".",child) == false
                    push!(graph["nodes"],Dict([("id",child),("colour","orange")]))
                end
            end
        end
    end
end
open("/home4/2191618l/Github/PangoTreeFiles/recombinant_graphs.json","w") do f
    JSON.print(f,graph)
end

variables =  Dict{String,Any}("max_estimated_count"=>max_count)
open("/home4/2191618l/Github/PangoTreeFiles/variables.json","w") do f
    JSON.print(f,variables)
end