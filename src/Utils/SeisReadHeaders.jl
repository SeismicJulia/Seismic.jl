include("Header.jl")

function SeisReadHeaders(filename)

    filename_h = join([filename ".seish"])
    stream_h = open(filename_h)
    seek(stream_h, header_count["n1"])
    nhead = 27
    nx = int(filesize(stream_h)/124)
    h = Array(Header,nx)
    for j = 1 : nx
        h[j] = GrabHeader(stream_h,j)
    end
    close(stream_h)
    return h
    
end
