function ind(nod,nnel,ndof)
    edof = nnel*ndof
    k=0
    index=[]
    for i=1:nnel
       start = (nod[i]-1)*ndof
       for j=1:ndof
         k=k+1
         append!(index,(start+j))
       end
    end
    return index
end
