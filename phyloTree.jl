
# Get seed mask based on weight and palindromic seed family
function get_seed_mask(weight::Int64)
    no_seeds = zeros(UInt32, 0)
    seed_mask_3 =  [
        0x0,0xb,   
        0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0
    ]
    seed_mask_4 =  [
        0x0,0x3b,  
        0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0
    ]
    seed_mask_5 =  [
        0x0,0x6b,  
        0x0,0x139,    
        0x0,0x193,
        0x0,0x6b, 
        0x0,0x0,  0x0,0x0
    ]
    seed_mask_6 =  [
        0x0,0x58D, #0b10110001101,
        0x0,0x653, #0b11001010011,
        0x0,0x1AB, #0b110101011,
        0x0,0xdb,  #0b11011011,
        0x0,0x0,  0x0,0x0
    ]
    seed_mask_7 =  [
        0x0,0x1953, #0b1100101010011
        0x0,0x588d, #0b101100010001101
        0x0,0x688b, #0b110100010001011
        0x0,0x17d,  #0b101111101,
        0x0,0x164d, #0b1011001001101,
        0x0,0x0,  0x0,0x0
    ]   
    seed_mask_8 =  [
        0x0,0x3927, #0b11100100100111,
        0x0,0x1CA7, #0b1110010100111,
        0x0,0x6553, #0b110010101010011,
        0x0,0xb6d,  #0b101101101101,
        0x0,0x0,  0x0,0x0
    ]
    seed_mask_9 =  [
        0x0,0x7497,  #0b111010010010111,
        0x0,0x1c927, #0b11100100100100111,
        0x0,0x72a7,  #0b111001010100111,
        0x0,0x6fb,   #0b11011111011,
        0x0,0x16ed,  #0b1011011101101,
        0x0,0x0
    ]
    seed_mask_10 = [
        0x0,0x1d297, #0b11101001010010111,
        0x0,0x3A497, #0b111010010010010111,
        0x0,0xE997,  #0b1110100110010111,
        0x0,0x6D5B,  #0b110110101011011,
        0x0,0x0,  0x0,0x0
    ]
    seed_mask_11 = [
        0x0,0x7954f,  #0b11110010101001111,
        0x0,0x75257,  #0b1110101001001010111,
        0x0,0x1c9527, #0b111001001010100100111,
        0x0,0x5bed,   #0b101101111101101,
        0x0,0x5b26d,  #0b1011011001001101101,
        0x0,0x0
    ]
    seed_mask_12 = [
        0x0,0x7954f,  #0b1111001010101001111,
        0x0,0x3D32F,  #0b111101001100101111,
        0x0,0x768B7,  #0b1110110100010110111,
        0x0,0x5B56D,  #0b1011011010101101101,
        0x0,0x0,  0x0,0x0
    ]
    seed_mask_13 = [
        0x0,0x792a4f, #0b11110010010101001001111,
        0x0,0x1d64d7, #0b111010110010011010111,
        0x0,0x1d3597, #0b111010011010110010111,
        0x0,0x1b7db,  #0b11011011111011011,
        0x0,0x75ad7,  #0b1110101101011010111,
        0x0,0x0
    ]
    seed_mask_14 = [
        0x0,0x1e6acf, #0b111100110101011001111,
        0x0,0xF59AF,  #0b11110101100110101111,
        0x0,0x3D4CAF, #0b1111010100110010101111,
        0x0,0x35AD6B, #0b1101011010110101101011,
        0x0,0x0,  0x0,0x0        
    ]
    seed_mask_15 = [
        0x0,0x7ac9af, #0b11110101100100110101111
        0x0,0x7b2a6f, #0b11110110010101001101111
        0x0,0x79aacf, #0b11110011010101011001111
        0x0,0x16df6d, #0b101101101111101101101
        0x0,0x6b5d6b, #0b11010110101110101101011
        0x0,0x0
    ]
    seed_mask_16 = [
        0x0,0xf599af, #0b111101011001100110101111,
        0x0,0xEE5A77, #0b111011100101101001110111,
        0x0,0x7CD59F, #0b11111001101010110011111,
        0x0,0xEB5AD7, #0b111010110101101011010111,
        0x0,0x0,  0x0,0x0
    ]
    seed_mask_17 = [
        0x0,0x6dbedb, #0b11011011011111011011011,
        0x0,0x0, 0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0       
    ]
    seed_mask_18 = [
        0x0,0x3E6B59F, #0b11111001101011010110011111,
        0x0,0x3EB335F, #0b11111010110011001101011111,
        0x0,0x7B3566F, #0b111101100110101011001101111,
        0x0,0x0,  0x0,0x0,  0x0,0x0
    ]
    seed_mask_19 = [
        0x0,0x7b974ef, #0b111101110010111010011101111
        0x0,0x7d6735f, #0b111110101100111001101011111
        0x0,0x1edd74f, #0b1111011011101011101101111
        0x0,0x0,  0x0,0x0,  0x0,0x0
    ]
    seed_mask_20 = [
        0x0,0x1F59B35F,#0b11111010110011011001101011111,
        0x0,0x3EDCEDF, #0b11111011011100111011011111,
        0x0,0xFAE675F, #0b1111101011100110011101011111,
        0x0,0x0,  0x0,0x0,  0x0,0x0
    ]
    seed_mask_21 = [
        0x0,0x7ddaddf, #0b111110111011010110111011111,
        0x0,0xaeb3f,   #0b11111100110101110101100111111,
        0x0,0x7eb76bf, #0b111111010110111011010111111,
        0x0,0x0,  0x0,0x0,  0x0,0x0
    ]
    seed_mask_22 = [0x0,0x003fffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_23 = [0x0,0x007fffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_24 = [0x0,0x00ffffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_25 = [0x0,0x01ffffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_26 = [0x0,0x03ffffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_27 = [0x0,0x07ffffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_28 = [0x0,0x0fffffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_29 = [0x0,0x1fffffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_30 = [0x0,0x3fffffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    seed_mask_31 = [0x0,0x7fffffff,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0,  0x0,0x0]
    
    seed_mask = [
        no_seeds,
        no_seeds,
        no_seeds,
        seed_mask_3,
        seed_mask_4,
        seed_mask_5,
        seed_mask_6,
        seed_mask_7,
        seed_mask_8,
        seed_mask_9,
        seed_mask_10,
        seed_mask_11,
        seed_mask_12,
        seed_mask_13,
        seed_mask_14,
        seed_mask_15,
        seed_mask_16,
        seed_mask_17,
        seed_mask_18,
        seed_mask_19,
        seed_mask_20,
        seed_mask_21,
        seed_mask_22,
        seed_mask_23,
        seed_mask_24,
        seed_mask_25,
        seed_mask_26,
        seed_mask_27,
        seed_mask_28,
        seed_mask_29,
        seed_mask_30,
        seed_mask_31,
    ]
    return seed_mask[weight]
end;

# All Ones
function get_solid_seed(weight::Int64)
    seed = 1
    seed <<= weight
    return seed
end;

function get_seed(weight::Int64, seed_rank::Int64)
    mask = get_seed_mask(weight + 1)
    low = mask[seed_rank*2 + 2]
    if (low == 0)
        return get_solid_seed(weight)
    end
    high = mask[seed_rank*2 + 1]
    seed = 0
    seed |= high
    seed <<=32
    seed |= low
    return seed
end

# From the paper and code, mer_size is calculated by lg_2(totalLen/num_seq)/1.5
# If mer_size is even, increment by 1
# mer_size is bounded by [5, 31]

function get_default_mer_size(total_len::Int64, num_seq::Int64)
    mer_size = round(Int64, (log2(round(Int64, total_len/num_seq)))/1.5)
    if iseven(mer_size)
        mer_size += 1
    end
    if mer_size < 5
        mer_size = 5
    end
    if mer_size > 31
        mer_size = 31
    end
    return mer_size
end

function get_default_seed(sequences, mer_size, coding_seed)
    totalLen = sum(length.(sequences))
    numSeq = length(sequences)
    if mer_size == 0
        mer_size = get_default_mer_size(totalLen, numSeq)
    end
    default_seed = get_seed(mer_size, coding_seed)
end;


function get_seed_length(seed::Int64)
    leftIndex = 64 - leading_zeros(seed)
    rightIndex = trailing_zeros(seed)
    seedLength = leftIndex - rightIndex
end

function apply_seed_mask(str::String, default_seed::Int64)
    bin_seed_mask = digits(default_seed, base=2)
    mer = ""
    for i in 1:get_seed_length(default_seed)
        if bin_seed_mask[i] == 1
            mer = mer*str[i]
        end
    end
    return mer
end

rev_map = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c');

#  Gets lexicographically smaller among the forward and reverse mers
function getFwdOrRevMer(fwd_mer::String)
    rev_mer = String("")
    for nucleotide in fwd_mer
        rev_mer = rev_mer*rev_map[nucleotide]
    end
    if isless(fwd_mer, rev_mer)
        return [fwd_mer, 0]
    end
    return [rev_mer, 1]
end

# Returns sorted and unique mer_list for a sequence
function get_sorted_mer_map(mer_list)
    # Sort the mer_list lexicographically
    sort!(mer_list, by=x->x[1]);
    mer_map = Dict{String, Array}()
    for val in mer_list
        lc_val = lowercase(val[1])  
        if haskey(mer_map, lc_val)
            push!(mer_map[lc_val], val[2:end])
        else
            push!(mer_map, lc_val=>[val[2:end]])
        end            
    end
    # Remove repeats
    for val in mer_map
        if (length(mer_map[val[1]])) >= 2
            delete!(mer_map, val[1])
        end
    end
    return mer_map
end

# Returns a hashmap from "mer" => list of [seq_index, mer_direction, mer_start, mer_end]].
# eg. of an entry: "atgtc" => Array{Any,1}[[4, 0, 1, 7] 
# Where [4, 0, 1, 7] means it's the mer of the 4th sequence, 
# having direction 0 or forward strand, starts at positon 1 and ends at positon 7
# Hashmap is used to keep track of number of matches.

function create_match_list(mer_map_list)
    match_list_map = Dict{String, Array}()
    for seq_index in 1:length(mer_map_list)
        mer_map = mer_map_list[seq_index]
        for key_val in mer_map
            key = key_val[1]
            val = key_val[2][1]
            if haskey(match_list_map, key)
                push!(match_list_map[key], vcat(seq_index, val))
            else
                push!(match_list_map, key=>[vcat(seq_index, val)])
            end
        end
    end 
    
    # Remove entries that have no matches
    for key_val in match_list_map
        if length(key_val[2]) < 2
            delete!(match_list_map, key_val[1])
        end
    end
    return match_list_map
end

# Returns forward or reverse complement nucleotide depending on the direction
function get_nucleotide(direction::Int64, nucleotide)
    if direction == 0
        return nucleotide
    else
        return rev_map[nucleotide]
    end
end

FWD_MASK = 0x1
BWD_MASK = 0x2
function extend_matches(sequences, match_list_map)
    for value in collect(values(match_list_map))
#         println(value)
        newValue = []
        
        fwd_step = 0
        bwd_step = 0 
        ext_dir = FWD_MASK | BWD_MASK
        new_start = ""
        new_end = ""
        
        while (ext_dir != 0)
            for mer_ind in 1:length(value)
                mer = value[mer_ind]
                seq_ind = mer[1]
                direction = mer[2]
                start_ind = mer[3]
                end_ind = mer[4]            
                cur_seq = sequences[seq_ind]
            
                # INITIALIZE
                if (mer_ind == 1)
                    # Check and increment in the backward direction
                    if (ext_dir & BWD_MASK == BWD_MASK && (start_ind - bwd_step > 1))
                        bwd_step += 1
                    else
                        ext_dir &= FWD_MASK
                    end
                    # Check and increment in the forward direction
                    if (ext_dir & FWD_MASK == FWD_MASK && (end_ind + fwd_step < length(cur_seq)))
                        fwd_step += 1
                    else
                        ext_dir &= BWD_MASK
                    end
                    # Assign new_start and new_end depending on the direction
                    new_start = get_nucleotide(direction, cur_seq[start_ind - bwd_step])
                    new_end = get_nucleotide(direction, cur_seq[end_ind + fwd_step])
                    
                    # After initialization move on to the next iteration
                    continue
                end
            
                # BACK EXTEND
                if (ext_dir & BWD_MASK == BWD_MASK)
                    if !(start_ind - bwd_step >= 1 && 
                            get_nucleotide(direction, cur_seq[start_ind - bwd_step]) == new_start)
                        bwd_step -= 1
                        new_start = get_nucleotide(direction, cur_seq[start_ind - bwd_step])
                        ext_dir &= FWD_MASK
                    end
                end

                # FWD EXTEND
                if (ext_dir & FWD_MASK == FWD_MASK)
                    if !(end_ind + fwd_step <= length(cur_seq) && 
                            get_nucleotide(direction, cur_seq[end_ind + fwd_step]) == new_end)
                        fwd_step -= 1
                        new_end = get_nucleotide(direction, cur_seq[start_ind + fwd_step])
                        ext_dir &= BWD_MASK
                    end  
                end
            
                # If we are not extending in either direction, break
                if (ext_dir == 0)
                    # println("Not extending in either direction")
                    break
                end 
            end
        end
        # Make the changes
        if fwd_step != 0 || bwd_step != 0
#             println("fwd_step is: ", fwd_step)
#             println("bwd_step is: ", bwd_step)
            for mer_ind in 1:length(value)
                value[mer_ind][3] = value[mer_ind][3] - bwd_step
                value[mer_ind][4] = value[mer_ind][4] + fwd_step
            end
        end
    end
end

# Returns a list of hashmaps ("Mer" => [Direction, Start-Position, End-Position])
# Direction is 0 if the mer is on forward strand and 1 if the mer is on the reverse complement strand
# Each entry is a hashmap of mers in each sequence
# So length of the returned value = Number of sequences

function get_mer_map_list(sequences::Array{String, 1}, default_seed)
    mer_map_list = Array{Any,1}()
    for seq in sequences
        mer_list = []
        mer_map = Dict{String, Array}()
        len = length(seq) - get_seed_length(default_seed) + 1
        for pos in 1:len
            str = seq[pos:pos+get_seed_length(default_seed)-1]
            fwd_mer = apply_seed_mask(str, default_seed) 
            mer = vcat(getFwdOrRevMer(fwd_mer), [pos, pos + get_seed_length(default_seed) - 1])
            push!(mer_list, mer)
        end  
        mer_map = get_sorted_mer_map(mer_list)
        push!(mer_map_list, mer_map)
    end  
    return mer_map_list
end

function get_pairwise_matches(sequences::Array{String, 1})
    DEFAULT_SEED_RANK = 3
    mer_size = 0
    default_seed = get_default_seed(sequences, mer_size, DEFAULT_SEED_RANK)

    # Apply default_seed to each position of each sequence to get a list of mers for each sequence
    mer_map_list = get_mer_map_list(sequences, default_seed)

    # match_list_map is a hashmap from "mer" => list of [seq_index, mer_direction, mer_start, mer_end]].
    # eg. of an entry: "atgtc" => Array{Any,1}[[4, 0, 1, 7] 
    # Where [4, 0, 1, 7] means it's the mer of the 4th sequence, 
    # having direction 0 or forward strand, starts at positon 1 and ends at positon 7
    match_list_map = Dict{String,Array}()
    match_list_map = create_match_list(mer_map_list)

    # Now, extend matches
    extend_matches(sequences, match_list_map)

    extended_match_list = Array{Array{Union{Nothing, Array{Int64,1}},1},1}()
    for value in collect(values(match_list_map))
        row_match_list = Array{Union{Nothing, Array{Int64,1}},1}()
        for match_index in 1:length(sequences)
            has_match_index = false
            for mer in value
                if match_index == mer[1]
                    push!(row_match_list, mer[2:end])
                    has_match_index = true
                end
            end
            if !has_match_index
                push!(row_match_list, nothing)
            end
        end
        push!(extended_match_list, row_match_list)
    end
    return extended_match_list
end

# Sorts matches to be longest to shortest
function sortMatches(matchA::Array{Union{Nothing, Array{Int64,1}},1}, 
        matchB::Array{Union{Nothing, Array{Int64,1}},1}) 
    firstA = findfirst(i->i!=nothing, matchA)
    firstB = findfirst(i->i!=nothing, matchB)
    lenA = matchA[firstA][3] - matchA[firstA][2]
    lenB = matchB[firstB][3] - matchB[firstB][2]
    return(lenA < lenB)
end

# Returns matches that exist in both pairwise sequences
function keepMatches(i::Int64, j::Int64, matches::Array{Array{Union{Nothing, Array{Int64,1}},1},1}) 
    cur_matches = Array{Union{Nothing, Array{Int64,1}},1}[]
    for match in matches
        if match[i]!=nothing && match[j]!=nothing
            push!(cur_matches, [deepcopy(match[i]), deepcopy(match[j])])
        end
    end
    return(cur_matches)
end


function trimMatches(a::Array{Union{Nothing, Array{Int64,1}},1}, 
        b::Array{Union{Nothing, Array{Int64,1}},1}) 
    #print("type: ", a, "\n")
    # Checking if the second match overlapped with the first match in either sequence
    while ( b[1][2] <= a[1][3] && b[1][2] >= a[1][2] ) || ( b[2][2] <= a[1][3] && b[2][2] >= a[1][3] )
        # move the start position up until it no longer overlaps with the other match
        b[1][2] += 1
        b[2][2] += 1
    end
    
    while ( b[1][3] <= a[1][3] && b[1][3] >= a[1][2] ) || ( b[2][3] <= a[1][3] && b[2][3] >= a[1][3] )
        # move the end position down until it no longer overlaps with the other match
        b[1][3] -= 1
        b[2][3] -= 1
    end
    
    if b[1][3] - b[1][2] <= 0 # a match has been completely removed
        b[1] = nothing
        b[2] = nothing
    end
end


# Calculates coverage -- can be weighted on which nucleotide has been substituted
function ntWeight(g1::String, g2::String, matches::Array{Array{Union{Nothing, Array{Int64,1}},1},1})
    coverageA = 0
    coverageB = 0
    
    strong = 1
    weak = 0.9
    purine_transition = .5
    pyrimidine_transition = .5
    transversion = 0
    
    for m in 1:length(matches)
        if matches[m][1] != nothing
            ntSeq2 = matches[m][2][2] # starting point nucleotide in the second sequence
            for ntSeq1 in matches[m][1][2]:matches[m][1][3] # starting point nucleotide in the first sequence
                if (g1[ntSeq1]=='a' && g2[ntSeq2]=='a') || (g1[ntSeq1]=='t' && g2[ntSeq2]=='t') # strong
                    coverageA += strong
                    coverageB += strong
                elseif (g1[ntSeq1]=='g' && g2[ntSeq2]=='g') || (g1[ntSeq1]=='c' && g2[ntSeq2]=="c") # weak
                    coverageA += weak
                    coverageB += weak
                elseif (g1[ntSeq1]=='a' && g2[ntSeq2]=='g') || (g1[ntSeq1]=='g' && g2[ntSeq2]=='a') # purine
                    coverageA += purine_transition
                    coverageB += purine_transition
                elseif (g1[ntSeq1]=='c' && g2[ntSeq2]=='t') || (g1[ntSeq1]=='t' && g2[ntSeq2]=='c') #pyrimidine
                    coverageA += pyrimidine_transition
                    coverageB += pyrimidine_transition
                else # transversion substitution
                    coverageA += transversion
                    coverageB += transversion
                end                
                ntSeq2 += 1                
            end
        end
    end    
    return coverageA, coverageB
    
end


########## Begin Processing Matches to Find Edit Distances ##########

function findDistances(sequences, matches)    
    # Make sure that the incoming sequences are all lowercase
    seqnum=1
    for seq in sequences
        sequences[seqnum] = lowercase(seq)
        seqnum += 1
    end    
    #initialize output matrix
    edit_distance_matrix = zeros(Float64, (length(sequences), length(sequences))) 

    # Step 1) Sort the matches by descending length 
    sort!(matches, lt=sortMatches, rev=true)

    # Pairwise Comparison
    for i in 1:length(sequences)                
        for j in (i+1):length(sequences)
            # Step 2) Subset the matches array to only contain matches that exist 
            #         for the current pair of sequences being compared
            cur_matches = keepMatches(i,j,matches)
            
            for k in 1:length(cur_matches)
                # for the case that an entire match was trimmed away in a previous round
                if cur_matches[k][1] == nothing 
                    continue
                end
                for l in (k+1):length(cur_matches)
                    if cur_matches[l][1] == nothing # for the case that an entire match was trimmed away
                        continue
                    end
                    # Step 3) Trim overlapping matches in each sequence
                    trimMatches(cur_matches[k], cur_matches[l])
                end
            end

            # Step 4) Calculate the coverage:
            #         Coverage within a match can be weighted based on the type of match/substitiution
            #         Nucleotides not in a match are counted as dissimilar.
            coverage = ntWeight(sequences[i], sequences[j], cur_matches)
            coverageA = coverage[1]
            coverageB = coverage[2]

            percentA = coverageA / length(sequences[i])
            percentB = coverageB / length(sequences[j])
            percentAB = (percentA + percentB) / 2   # average the coverage
            editPercent = 1 - percentAB 
            
            # Step 5) Set the edit distance for the pairwise genomes
            edit_distance_matrix[i, j] = editPercent

            # Step 6) Repeat until all pairs have a distance calculated for them
        end
    end

    # Convert to a matrix from array of arrays and fill the lower half
    row_len = size(edit_distance_matrix, 1)
    for i = 1:row_len, j=1:row_len
        if i>j
            edit_distance_matrix[i, j] = edit_distance_matrix[j, i]
        end
    end
    return edit_distance_matrix
end



# Defining a "Tree" composite type
module TreeMod
    mutable struct Tree
        name::String
        newick::String
        par::Any
        branch::Any
        chil::Dict
        connections::Dict
    
        Tree(name, newick) = new(name, newick)
        function Base.show(io::IO, t::Tree)
            println(t.name)
        end
    end
end

function MakeAnchorTree(D, debug::Bool = false)
    # Make Q1 Matrix
    function QMat(A)
        m = size(A,1)
        Q = copy(A)
        function q(i,j)
            if i==j
                return 0
            else
                n = size(D,1) - 1
                dij = D[i,j]
                sum_dik = sum(D[i,:])
                sum_djk = sum(D[:,j])
                return (m-2) * dij - sum_dik - sum_djk
            end
        end
    
        for i = 1:m
            for j = 1:m
                Q[i,j] = q(i,j)
            end
        end
        return Q
    end

    # Get the distance between a node k and nodes f and g
    function dist(k,f,g)
        dfk = D[f,k]
        dgk = D[g,k]
        dfg = D[f,g]
        return (1/2) * (dfk + dgk - dfg)
    end

    # Create the nodes, add to the tree, update their connections
    function make_nodes(f,g,index_trees)
        # Create the tree nodes for f and g, and then create a new neighbor node u that connects the two
        s = index_trees[f]
        t = index_trees[g]
        uname = "neighbor"* string(r)
        u = TreeMod.Tree(uname, uname) # name neighbor nodes sequentially

        # compute distance from f and g to u -- returns array of dist(f), dist(g)
        node_dists = deltas(f,g)
        #println("node distances: ", node_dists)
        
        # update connections list with branch lengths
        if isdefined(u, :connections)
            merge!(u.connections, Dict(s => node_dists[1], t => node_dists[2]))
        else
            u.connections = Dict(s => node_dists[1], t => node_dists[2])
        end

        if isdefined(s, :connections)
            merge!(s.connections, Dict(u => node_dists[1]))
        else
            s.connections = Dict(u => node_dists[1])
        end
        
        if isdefined(t, :connections)
            merge!(t.connections, Dict(u => node_dists[2]))
        else
            t.connections = Dict(u => node_dists[2])
        end
        # update the new node's Newick string
        newick = "($(s.newick):$(node_dists[1]), $(t.newick):$(node_dists[2]))$(u.newick)"
        u.newick = newick
        #println("new node's newick: ",newick)

        # update list of trees and new index list
        #println(index)
        #println(index_trees)
        index_trees = [index_trees[i] for i in index if i != f && i != g]

        pushfirst!(index_trees, u)
        push!(tree, u)

        return index_trees, tree
    end

    # Calculate distance nodes f and g to the connecting neighbor node
    function deltas(f,g)
        dfg = D[f,g]
        sum_dfk = sum(D[f,:])
        sum_dgk = sum(D[:,g])
        deltaf = ((1/2) * dfg + 1/(2(m-2)) * (sum_dfk - sum_dgk))
        deltag = dfg - deltaf
        return [deltaf,deltag]
    end

    index = [i for i in 1:n]
    index_ids = ["seq"*string(i) for i in index]

    global index_trees = [TreeMod.Tree(i, i) for i in index_ids]
    global tree
    tree = copy(index_trees)
               
    # Start out with n nodes that need to be connected and an n+1xn+1 matrix D
    m = n
    D = D
    r = 1

    # Until we get down to 2 nodes
    while m > 3
        if debug
            println("\nNEW ITERATION!!!\n====================")
            println("Number of nodes still to connect (m): ", m," ...at iteration: ",r)
            println("\nFinding sequences to connect\n====================")
            println("Converting absolute distances from D to proportional (?) distances in Q")
            display(D)
        end
                                
        Q = QMat(D)
        if debug display(Q) end

        # Find the smallest distance between two sequences (smallest Q value)
        farthest_seqs = argmin(Q)

        # Indices of f and g are the indices of D and Q where these values are
        f = farthest_seqs[1]
        g = farthest_seqs[2]

        if debug println("farthest apart nodes are f ($(index_ids[f])) at index $f and g ($(index_ids[g])) at index $g") end
        
        # make the nodes for f and g, connect them to a new internal node, and add distances and connections
        if debug println("\nConnecting the sequences\n====================") end
        index_trees, tree = make_nodes(f,g,index_trees)
        
        # if we have more nodes left to connect, continue on and prepare new distance matrix
        if debug println("\nUpdating distance matrix and indices\n====================") end
        # Initialize a new n-1 x n-1 distance matrix of zeros (because have joined two nodes)
        D_new = zeros(m-1,m-1)

        index_new = [i for i in index if i != g && i != f]
        index_ids_new = [index_ids[i] for i in index if i != g && i != f]

        # Make a new array without the rows and columns of nodes we have already connected
        D_old = D[setdiff(1:end, (f,g)), setdiff(1:end, (f,g))]

        # Distances not affected by the join are the same as in D
        # So place the values from the D_old matrix into the D_new matrix, after the first row
        # and column, which are reserved for the newly created node.      
        for i in 2:m-1
            for j in 2:m-1
              D_new[i,j] = D_old[i-1, j-1]
            end
        end

        # Calculate distances from the new node to the other nodes
        for i in 1:m-2            
            old_idx = index_new[i]
            thedist = dist(old_idx,f,g)
            D_new[1,i+1] = thedist
            D_new[i+1,1] = thedist
        end

        # Update the distance index:
        # First, remove the two sequences we have already placed from our index
        # and replace with a dummy index called 0 for the joined node
        index_clean = [i for i in 1:m-1]
        index = index_clean
        index_ids = pushfirst!(index_ids_new, "internal")

        # Make D_new our main distance matrix
        D = D_new

        # Update m
        m -= 1
        r += 1
    end

    # 3 nodes left... join the first two nodes s and t at a new node u
    s = index_trees[1]
    t = index_trees[2]
    uname = "neighbor"* string(r)
    node_dists = deltas(1,2)
    u = TreeMod.Tree(uname, uname) # name neighbor nodes sequentially
    
    # Add connections between s and u, t and u                                    
    if isdefined(s, :connections)
        merge!(s.connections, Dict(u => node_dists[1]))
    else
        s.connections = Dict(u => node_dists[1])
    end
        
    if isdefined(t, :connections)
        merge!(t.connections, Dict(u => node_dists[2]))
    else
        t.connections = Dict(u => node_dists[2])
    end

    if isdefined(u, :connections)
        merge!(u.connections, Dict(s => node_dists[1], t => node_dists[2]))
    else
        u.connections = Dict(s => node_dists[1], t => node_dists[2])
    end

    # Then connect the last node r with the newly created node u
    r = index_trees[3]
    r_u_dist = dist(3, 1, 2)

    if isdefined(r, :connections)
        merge!(r.connections, Dict(u => r_u_dist))
    else
        r.connections = Dict(u => r_u_dist)
    end

    merge!(u.connections, Dict(r => r_u_dist))
    
    # Add the new node to the tree
    push!(tree, u)
                                        
    return tree
end


function reconstruct_path(cameFrom, current)
    total_path = [current]
    total_sum = 0
    while current in collect(keys(cameFrom))
        current = cameFrom[current]
        push!(total_path, current)
    end
    # println("total path: $total_path")
    return total_path

end

function A_Star(start, goal, tree)
    # The set of nodes already evaluated
    closedSet = Array{Any, 1}()

    # The set of currently discovered nodes that are not evaluated yet.
    # Initially, only the start node is known.
    openSet = [start]

    # For each node, which node it can most efficiently be reached from.
    # If a node can be reached from many nodes, cameFrom will eventually contain the
    # most efficient previous step.
    cameFrom = Dict()

    # For each node, the cost of getting from the start node to that node.
    # Initialized at infinity
    gScore = Dict(i => Inf for i in tree) 

    # The cost of going from start to start is zero.
    gScore[start] = 0

    # For each node, the total cost of getting from the start node to the goal
    # by passing by that node. That value is partly known, partly heuristic.
    fScore = Dict()

    # For the first node, that value is completely heuristic.
    # But we don't have a heuristic function for this so we set it to 0
    fScore[start] = 0

    while !isempty(openSet)
        # println("============")
        
        scores = Dict(i => fScore[i] for i in openSet)
        current = findmin(scores)[2]
        if current == goal
            return [reconstruct_path(cameFrom, current), gScore[current]]
        end
                                                
        deleteat!(openSet, findfirst(isequal(current), openSet))
        push!(closedSet, current)

        # println("For neighbors -----")
        for neighbor in keys(current.connections)
            # println(neighbor)
            if neighbor in closedSet
                continue        # Ignore the neighbor which is already evaluated.
            end

            # The distance from start to a neighbor
            tentative_gScore = gScore[current] + current.connections[neighbor]

            if !(neighbor in openSet)   # Discover a new node
                push!(openSet, neighbor)
            elseif tentative_gScore >= gScore[neighbor]
                continue     
            end

            # This path is the best until now. Record it!
            cameFrom[neighbor] = current
            gScore[neighbor] = tentative_gScore
            fScore[neighbor] = gScore[neighbor] + 0
        end
    end
end

function getRootID(midpath)
    traveled = 0
    for i in 1:length(midpath)-1
        current = midpath[i]
        next = midpath[i+1]
        dist_to_next = current.connections[next]
        new_traveled = traveled + dist_to_next
        
        # If after going to the next node, we won't be at the midpoint yet, do it                              
        if new_traveled < midpoint
            current = next
            traveled = new_traveled
        
        # If after going to the next node, we will have passed the midpoint, stop and root                                           
        elseif new_traveled > midpoint
            # Create root node
            root = TreeMod.Tree("root", "root")
            push!(tree, root)
                                                    
            # Specify the parent => child relationship between root and the two nodes it is. between                                        
            root.chil = Dict(current => midpoint-traveled, next=> new_traveled-midpoint)
            current.par = Dict(root => midpoint-traveled)
            next.par = Dict(root => new_traveled-midpoint)
            
            # Root node is now between "curent" and "next" nodes, so remove eachother from connections
            delete!(current.connections, next)
            delete!(next.connections, current)
            
            # Return an array with the new root, and its two children
            return([current, next, root])
            break
        
        # In the unlikely event that the midpoint is right on the next node, rename and return that as the root
        elseif new_traveled == midpoint
            next.name = "root"
            root = current
            return(current)
            break
        end
    end
end

function update_children(not_an_orphan, childless)
    while !isempty(childless)
        just_got_a_parent = Array{Any, 1}()
        for node in childless
            #println("parent: $node")
            if isdefined(node, :connections)
                for c in keys(node.connections)
                    if !(c in not_an_orphan)
                        branch_l = node.connections[c]
                        c.par = Dict(node => branch_l)
                        if isdefined(node, :chil)
                            merge!(node.chil, Dict(c => branch_l))
                        else
                            node.chil = Dict(c => branch_l)
                        end
                        push!(not_an_orphan, c)
                        push!(just_got_a_parent, c)
                        push!(childless, c)
                    end
                end
            end
            deleteat!(childless, findfirst(isequal(node), childless))
        end
    end
end

function create_phyogenetic_tree(distanceMatrix, debug = false)
    global n                                        
    n = size(distanceMatrix, 1)                                  
    tree = MakeAnchorTree(distanceMatrix, false)

    node_dict = Dict()
    for node in tree
        merge!(node_dict, Dict(node.name => node))
    end

    # initiate an array of paths from each leaf
    leaves = [t for t in tree if length(keys(t.connections)) == 1]
    no_leaves = length(leaves)
    leaf_dists = Dict()
    leaf_paths = Dict()

    for i = 1:no_leaves, j = 1:no_leaves
        if i > j
            A = A_Star(leaves[i], leaves[j], tree)
            leaf_dists[(i,j)] = A[2]
            leaf_paths[(i,j)] = A[1]
        end
    end

    max_pairdist = findmax(leaf_dists)
    global midpoint                                                    
    midpoint = max_pairdist[1]/2

    p = max_pairdist[2][1]
    q = max_pairdist[2][2]

    midpath = leaf_paths[max_pairdist[2]]
    rooted = getRootID(midpath)

    if length(rooted) > 1
        root = rooted[3]
        childless = [rooted[1], rooted[2]]
    else
        root = rooted
        childless = Array{Any, 1}(root)
    end

    # Don't try to parent these connections -- they have already been seeen
    not_an_orphan = rooted
    update_children(not_an_orphan, childless)

    leaves = [i for i in tree if !isdefined(i, :chil)]
    for i in tree
        if isdefined(i, :par)
            p = i.par
        else
            p = i
        end
    end

    to_n = Array{Any, 1}()
    for i in leaves
        i.newick = i.name
        push!(to_n, collect(keys(i.par))[1])
    end

    while !isempty(to_n)
        p = pop!(to_n)
        children = collect(keys(p.chil))
        c1 = children[1]
        c2 = children[2]                                                                  
        p.newick = "($(c1.newick):$(p.chil[c1]),$(c2.newick):$(p.chil[c2]))" * (occursin("seq", p.name) ? "$(p.name)" : "")
        if p.name != "root"
            pushfirst!(to_n, collect(keys(p.par))[1])
        end
    end
    println(root.newick*";")
    return root.newick*";"
end

function write_file(buffer::String, file_name::String="outputfile.nwk")
    output_file = open(file_name, "w")
    write(output_file, buffer)
    close(output_file)
    println("Wrote ", file_name)
end

phyloTree = String("")
function main(sequences)
    pairwise_matches = get_pairwise_matches(sequences);
    distanceMatrix = findDistances(sequences, pairwise_matches);

    global phyloTree = create_phyogenetic_tree(distanceMatrix, false)
end

inp_file_name = ARGS[1]
println("Processing ", inp_file_name)
sequences_with_extra_newline = readlines(open(inp_file_name))


sequences = Array{String, 1}()
for seq in sequences_with_extra_newline
    if length(seq) > 0
        push!(sequences, seq)
    end
end

# Running it for first time
# main(sequences)
# # write_tree_file(phyloTree)
# write_file(phyloTree, inp_file_name*"_out.nwk")


# Profiling it when it is run the second time
start_time = time()
mem_alloc = @allocated main(sequences)
end_time = time()
elapsed_time = Float16(end_time - start_time)

stats = "Statistics:\n";
stats *= "ElapsedTime: "*string(elapsed_time)*"\n"
stats *= "MemAlloc: "*string(mem_alloc)*"\n"
write_file(stats, inp_file_name*"_stats.txt")