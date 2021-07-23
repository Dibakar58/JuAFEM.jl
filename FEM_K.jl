function FEM_K(Elements,Nodes)
    (Nel,nnode)=size(Elements)
    (Node,np)=size(Nodes)
    dof=Node*np
    ndof=2
    K_Be=Dict()

    for el=1:Nel
        elem_connect=Elements[el,:]

        x1=Nodes[elem_connect[1],1]
        x2=Nodes[elem_connect[2],1]
        x3=Nodes[elem_connect[3],1]

        y1=Nodes[elem_connect[1],2]
        y2=Nodes[elem_connect[2],2]
        y3=Nodes[elem_connect[3],2]
        #create module
        area2=Area(x1,x2,x3,y1,y2,y3)

        B=zeros(Float64,3,6)


        dNdx=(1/area2)*[(y2-y3),(y3-y1),(y1-y2)]
        dNdy=(1/area2)*[(x3-x2),(x1-x3),(x2-x1)]
        ## use bmat_2d
        for i=1:3
            i1=(i-1)*2+1;
            i2=i1+1;
            B[1,i1]=dNdx[i]
            B[2,i2]=dNdy[i]
            B[3,i1]=dNdy[i]
            B[3,i2]=dNdx[i]
        end

        ### for index

        k=0
        index=zeros(Int,6)
        for i=1:3
            start = (elem_connect[i]-1)*2;
            for j=1:2
                k=k+1
                index[k]=start+j
            end
        end

        K_Be[el]=B
   end

    return K_Be

end
