#Program to treat the data which are obtained from orientation.f95
#The information provided is the angle between the two ligands of a central atom according to its projection on the normal axis

function acos(x) { return atan2(sqrt(1-x*x), x) }

BEGIN{
    #Cell parameters
#    a=11.5888
#    b=13.3816
#    c=34.
    
    #Ratio to convert z atomic position into the associated the cells for tabz
#    div=0.25    #Space division according to x,y and z
    ratio=(int(param/div)+1)/param

    #Initialization of the table
    for(i=0;i<=int(param/div)+1;i++){
	tab_ref[i,1]=0
	tab_ref[i,2]=0
	norm[i]=0
    }

    rad2deg=45/atan2(1,1)

}

{
    if((getline f1 < axis) <= 0 || (getline u1 < x1) <= 0 || (getline v1 < y1) <= 0 || (getline w1 < z1) <= 0 || (getline u2 < x2) <= 0 || (getline v2 < y2) <= 0)#Recording the data of the 2 other files
	exit 
    else if(NR>2){
	split(f1,tab1," ")
	split(u1,tab2," ")
	split(v1,tab3," ")
	split(w1,tab4," ")
	split(u2,tab5," ")
	split(v2,tab6," ")
	for(i=2;i<=NF;i++){
	    if($i!="x"){
		tab_ref[int(tab1[i]*ratio),1]+=acos(((tab2[i]*tab5[i])+(tab3[i]*tab6[i])+(tab4[i]*$i))/sqrt((tab2[i]^2+tab3[i]^2+tab4[i]^2)*(tab5[i]^2+tab6[i]^2+$i^2))) #We sum the angles
		tab_ref[int(tab1[i]*ratio),2]+=sqrt(tab5[i]^2+tab6[i]^2+$i^2)
		norm[int(tab1[i]*ratio)]++
	    }
	}
    }
    else if(NR==2){
	nb_at=NF-1
    }
}


END{
    print "#proj_central_atom(Ang)   angle(L1-center-L2)(deg)   d(center-L2)(Ang)"
    for(i=0;i<=int(param/div)+1;i++){
	if(norm[i]!=0){#If there is at least one atom in this area...
	    print i/ratio" "tab_ref[i,1]/norm[i]*rad2deg" "tab_ref[i,2]/norm[i]#We print the average angle between the two ligands and the average length bond of ligand 2
	}
    }
}
