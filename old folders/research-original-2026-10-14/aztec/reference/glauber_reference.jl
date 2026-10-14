# Historical exploratory prototype. The active, tested implementation lives in
# `aztec/src/GlauberSquareGrid.jl`; this file is retained only to document the
# original boundary and local-update conventions.

using SharedArrays
using Distributed

@everywhere function maxf(i,j,n)
    if i==0 || i==n
        return j%2
    elseif j==0 || j==n
        return -(i%2)
    else
        return 0
    end
end

@everywhere function intonemax(i,j,k,L)
    if i==k || i==2*L-k || j==k || j==2*L-k
        return 2*k+maxf(i-k,j-k,2*L-2*k)
    else
        return 0
    end
end

@everywhere function intonemin(i,j,k,L)
    if i==k || i==2*L-k || j==k || j==2*L-k
        return -2*k+maxf(i-k,j-k,2*L-2*k)
    else
        return 0
    end
end


#Gives the maximum configuration on the square
@everywhere Amax(L)=reshape([intonemax(i,j,min(i,j,2*L-i,2*L-j),L) for i=0:2L for j=0:2L],2L+1,2L+1)

#Gives the mininum configuration on the square

@everywhere Amin(L)=reshape([intonemin(i,j,min(i,j,2*L-i,2*L-j),L) for i=0:2L for j=0:2L],2L+1,2L+1)


@everywhere function mcmove2(mat,a)
    pos=rand(2:(size(mat)[1]-1),2)
    if (mat[pos[1],pos[2]]- mat[pos[1]-1,pos[2]])*(mat[pos[1],pos[2]]- mat[pos[1]+1,pos[2]])==9 || (mat[pos[1],pos[2]]- mat[pos[1],pos[2]+1])*(mat[pos[1],pos[2]]- mat[pos[1],pos[2]-1])==9
            if (pos[1]+pos[2])%2==0
                    if rand()<0.5
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]-1
                    else
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]+3
                end
            elseif pos[1]%2==0
                    if rand()<a*a/(1.0+a*a)
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]-3
                    else
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]+1
                end
            elseif pos[1]%2==1
                    if rand()<1.0/(1.0+a*a)
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]-3
                    else
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]+1
                end
            end
               end
        return mat
end

@everywhere function dynamicscoal(L,a)
    A1=Amax(L)
    A2=Amin(L)
    N=0
    while A1[L,L]>A2[L,L]
        A1=mcmove2(A1,a)
        A2=mcmove2(A2,a)
    N=N+1
    end
    return N
end

function mccoalpar(filename,N,L,a)
    val=SharedArray{Int64}(N)
       @distributed (+) for i=1:N
        val[i]=dynamicscoal(L,a)
        end
    f=open(filename,"w")
        for i=1:N
        write(f,string(val[i]))
        write(f,"\n")
    end
end




function mccoal(filename,N,L,a)
    f=open(filename,"w")
    for i=1:N
        write(f,string(dynamicscoal(L,a)))
        write(f,"\n")
    end
end

@everywhere function dynamicscoalwhole(L,a)
        A1=Amax(L)
        A2=Amin(L)
    N=0
    while all(A1[L+1,L+1].>A2[L+1,L+1])
            A1=mcmove2(A1,a)
            A2=mcmove2(A2,a)
        end
        while any(A1[2:2L,2:2L].>A2[2:2L,2:2L])
            A1=mcmove2(A1,a)
            A2=mcmove2(A2,a)
        N=N+1
        end
        return N
end

function mccoalparwhole(filename,N,L,a)
    val=SharedArray{Int64}(N)
       @distributed (+) for i=1:N
        val[i]=dynamicscoalwhole(L,a)
        end
    f=open(filename,"w")
        for i=1:N
        write(f,string(val[i]))
        write(f,"\n")
    end
    close(f)
end





function spot(mat)
        A=[]
    L=size(mat)[1]-1
        for i=2:L
            for j=2:L
            if (mat[i,j]-mat[i,j+1])*(mat[i,j]-mat[i,j-1])==9 || (mat[i,j]-mat[i+1,j])*(mat[i,j]-mat[i-1,j])==9
                        A=append!(A,[i,j])
                end
            end
        end
        return A
    end


function mcmove(mat,pos,a)
        if (pos[1]+pos[2])%2==0
                    if rand()<0.5
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]-1
                    else
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]+3
                end
            elseif pos[1]%2==0
                    if rand()<a*a/(1.0+a*a)
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]-3
                    else
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]+1
                end
            elseif pos[1]%2==1
                    if rand()<1.0/(1.0+a*a)
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]-3
                    else
                        mat[pos[1],pos[2]]= mat[pos[1]-1,pos[2]]+1
                end
            end
        return mat
end

function dynamicsNsweepprint(filename,L,N,a)
        f=open(filename,"w")
    m=0
    A=Amax(L)
    while m<N
        write(f,"{")
        for i=1:(size(A)[1])
            write(f,"{")
                for j=1:(size(A)[1])
                    if j!= (size(A)[1])
                        write(f,string(A[i,j]))
                        write(f,",")
                    else
                        write(f,string(A[i,j]))
                        end
                end
            if i!=(size(A)[1])
                write(f,"},")
            else
                write(f,"}}")
                write(f,"\n")
                end
            end
        r=spot(A)
        x=trunc(Int,rand(1:size(r)[1]/2))
        y=[r[2x-1],r[2x]]
        A=mcmove(A,y,a)
        m=m+1
        end
        close(f)
end

#
#
#write(f,string(val[i]))
#
