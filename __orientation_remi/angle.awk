#Program to treat the data which are obtained from orientation.f95
#The information provided is the angle between the two ligands of a central atom

function acos(x) { return atan2(sqrt(1-x*x), x) }

BEGIN{
    #Cell parameters
#    a=11.5888
#    b=13.3816
#    c=34.
    
    #Ratio to convert x/y/z atomic position into the associated the cells for tabxyz
#    div=0.2    #Space division according to x,y and z
    ratio_a=(int(a/div)+1)/a
    ratio_b=(int(b/div)+1)/b
    ratio_c=(int(c/div)+1)/c

    #Initialization of the table
    for(i=0;i<=int(a/div)+1;i++){
	for(j=0;j<=int(b/div)+1;j++){
	    for(k=0;k<=int(c/div)+1;k++){
		tab_xyz[i,j,k]=0
		norm[i,j,k]=0
	    }
	}
    }

    rad2deg=45/atan2(1,1)
}

{
    if((getline f1 < x0) <= 0 || (getline f2 < y0) <= 0 || (getline f3 < z0) <= 0 || (getline u1 < x1) <= 0 || (getline v1 < y1) <= 0 || (getline w1 < z1) <= 0 || (getline u2 < x2) <= 0 || (getline v2 < y2) <= 0)#Recording the data of the 2 other files
	exit 
    else if(NR>2){
	split(f1,tab1," ")
	split(f2,tab2," ")
	split(f3,tab3," ")
	split(u1,tab4," ")
	split(v1,tab5," ")
	split(w1,tab6," ")
	split(u2,tab7," ")
	split(v2,tab8," ")
	for(i=2;i<=NF;i++){
	    if($i!="x"){
		tab_xyz[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c)]+=acos(((tab4[i]*tab7[i])+(tab5[i]*tab8[i])+(tab6[i]*$i))/sqrt((tab4[i]^2+tab5[i]^2+tab6[i]^2)*(tab7[i]^2+tab8[i]^2+$i^2)))#We sum the angles
		norm[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c)]++#Normalization factor to calculate the average
	    }
	}
    }
    else if(NR==2){
	nb_at=NF-1
    }
}

END{
    print "#x0(Ang)   y0(Ang)   z0(Ang)   angle(L1-center-L2)(deg)"
    for(i=0;i<=int(a/div)+1;i++){
	for(j=0;j<=int(b/div)+1;j++){
	    for(k=0;k<=int(c/div)+1;k++){
		if(norm[i,j,k]!=0){#If there is at least one atom in this area...
		    print i/ratio_a" "j/ratio_b" "k/ratio_c" "tab_xyz[i,j,k]/norm[i,j,k]*rad2deg#We print the average angle between the two ligands
		}
	    }
	}
    }
}
