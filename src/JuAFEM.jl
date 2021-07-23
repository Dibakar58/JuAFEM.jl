module AFEM

using LinearAlgebra
using DelimitedFiles

include("AREA.jl")
include("C.jl")
#include("DATA.jl")
include("FEM_K.jl")
include("ind.jl")

export Load
export fem
export nsfem
export afem
export solve

function Load(A,B)
    Elements=readdlm(A,'\t','\n')
    Elements=Array{Int64,2}(Elements)
    #Elements=Elements[:,2:4]

    Nodes=readdlm(B,'\t')
    Nodes=Array{Float64,2}(Nodes)
    #Nodes=Nodes[:,2:3]

    return Elements , Nodes

end

function fem(Elements,Nodes,nu,E)
    D=C(nu,E)
    (Nel,nnode)=size(Elements)
    (Node,np)=size(Nodes)
    dof=Node*np
    ndof=2
    K=zeros(Float64,dof,dof)

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
        BT=transpose(B)
        Ke=(BT*D*B)*(area2/2.0)


        K[index,index]=K[index,index]+Ke

    end

    return K

end

function nsfem(Elements,Nodes,nu,E)
    l=C(nu,E)
    K_Be=FEM_K(Elements,Nodes)
    node_adj=Dict()
    Nnode=size(Nodes)[1]
    Nel=size(Elements)[1]

    K=zeros(Nnode*2,Nnode*2)

    for i=1:Nnode
        ch=[]
        for j=1:Nel

            for check=1:3
                if i==Elements[j,check]
                    append!(ch,j)
                    continue

                end

            end

        end
        node_adj[i]=ch
    end
    area=[]
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
        append!(area,area2/2)
    end


    area_nod=zeros(Nnode,1)
    for i=1:Nnode
        ele_con=node_adj[i]
        for j=1:size(ele_con)[1]
            area_nod[i]=area_nod[i]+(area[ele_con[j]])/3
        end
    end


    for ino=1:Nnode
        ino_adj=node_adj[ino]
        area_ino=area_nod[ino]
        #println("ino_adj ==>>",ino_adj)
        ne=size(ino_adj,1)
        #println("Started =>>" , ino)
        nn=0
        nodB=[]
        for ie=1:ne
            ie_nods=Elements[ino_adj[ie],:]
            #println("ie_nods ==>>" , ie_nods)
            nn_ie=size(ie_nods,1)

            if ie==1
                #println("Inside if ie==1")
                append!(nodB,Elements[ino_adj[ie],:])
                nn=nn_ie
                for j=1:nn
                    if j==1
                        B_ino=1/3*area[ino_adj[ie]]*K_Be[ino_adj[ie]][1:3,2*j-1:2*j]
                    else
                        hj=1/3*area[ino_adj[ie]]*K_Be[ino_adj[ie]][1:3,2*j-1:2*j]
                        B_ino=hcat(B_ino,hj)
                    end
                end
                #B_ino=(1/(3*area[ino_adj[ie]]))*K_Be[ino_adj[ie]]
                #println("B_ino Done ==>>" , ie)
                #println(nn)
                #println(size(B_ino))
            else
                i0=0
                #println("Inside else")
                for iino=1:nn_ie
                    nod=Elements[ino_adj[ie],:][iino]
                    flag=0
                    #println(nn_ie)
                    #println(nn)
                    for j=1:nn
                        if nodB[j]==nod

                            B_ino[1:3,2*j-1:2*j]=B_ino[1:3,2*j-1:2*j] + 1/3*area[ino_adj[ie]]*(K_Be[ino_adj[ie]][1:3,2*iino-1:2*iino])
                            #println("B_ino done ==>>",iino)
                            flag=1
                            break
                        end
                    end

                    if flag==0
                        i0=i0+1
                        #println("Nod ==>>" , nod)

                        append!(nodB,nod)
                        #println("nodB ==>>" , nodB)
                        bi=1/3*area[ino_adj[ie]]*K_Be[ino_adj[ie]][1:3,2*iino-1:2*iino]
                        #println("bi ====>>>" , bi)
                        B_ino=hcat(B_ino,bi)
                    end
                end
                nn=nn+i0
            end

        end
        B_ino=B_ino/area_ino

        #println(B_ino)
        #println(nodB)
        BT_ino=transpose(B_ino)
        Ke=(BT_ino*l*B_ino)*(area_ino)

        index=ind(nodB,size(nodB,1),2)
        #println(index)
        K[index,index]=K[index,index]+Ke

    end

    return K

end



function afem(KFEM,KNSFEM,alpha)
    K=KNSFEM.*(1-alpha*alpha)+KFEM.*(alpha*alpha)
    return K
end

function solve(K,F)
    return K\F
end
end
