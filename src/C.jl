function C(nu,E)
    E1=E/(1-nu*nu)

    C=zeros(Float64,3,3)
    C[1,1]=E1
    C[1,2]=nu*E1
    C[2,2]=E1
    C[2,1]=nu*E1
    C[3,3]=(1-nu)*(0.5)*E1


    return C
end
