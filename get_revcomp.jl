function get_comp(c::Char)
    if c=='A' return 'T'
    elseif c=='T' return 'A'
    elseif c=='C' return 'G'
    elseif c=='G' return 'C'
    else
        println("$c is not A, T, C, or G.")
        return 'N'
    end
end

get_revcomp(s::String) = reverse(join([get_comp(c) for c in s]))