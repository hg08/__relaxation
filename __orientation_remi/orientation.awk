#Program to treat the data which are obtained from orientation.f95 (certainly to be used with gnuplot)
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
		tab_xyz[i,j,k,1]=0
		tab_xyz[i,j,k,2]=0
		tab_xyz[i,j,k,3]=0
		tab_xyz[i,j,k,4]=0
		tab_xyz[i,j,k,5]=0
		tab_xyz[i,j,k,6]=0
		norm[i,j,k]=0
	    }
	}
    }
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
		tab_xyz[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c),1]+=tab4[i]#We sum the vectors u1
		tab_xyz[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c),2]+=tab5[i]#We sum the vectors v1
		tab_xyz[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c),3]+=tab6[i]#We sum the vectors w1
		tab_xyz[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c),4]+=tab7[i]#We sum the vectors u2
		tab_xyz[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c),5]+=tab8[i]#We sum the vectors v2
		tab_xyz[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c),6]+=$i     #We sum the vectors w2
		norm[int(tab1[i]*ratio_a),int(tab2[i]*ratio_b),int(tab3[i]*ratio_c)]++#Normalization factor to calculate the average
	    }
	}
    }
    else if(NR==2){
	nb_at=NF-1
    }
}

END{
    #We are looking for the maximal value of norm
    max=0
    for(i=0;i<=int(a/div)+1;i++){
	for(j=0;j<=int(b/div)+1;j++){
	    for(k=0;k<=int(c/div)+1;k++){
		if(max<norm[i,j,k]){max=norm[i,j,k]}
	    }
	}
    }

    print "#x0(Ang)   y0(Ang)   z0(Ang)   v(x0-xL1)(Ang)   v(y0-yL1)(Ang)   v(z0-zL1)(Ang)   v(x0-xL2)(Ang)   v(y0-yL2)(Ang)   v(z0-zL2)(Ang)   occurence(highest occurence = 1)"
    for(i=0;i<=int(a/div)+1;i++){
	for(j=0;j<=int(b/div)+1;j++){
	    for(k=0;k<=int(c/div)+1;k++){
		if(norm[i,j,k]!=0){#If there is at least one atom in this area...
		    print i/ratio_a" "j/ratio_b" "k/ratio_c" "tab_xyz[i,j,k,1]/norm[i,j,k]" "tab_xyz[i,j,k,2]/norm[i,j,k]" "tab_xyz[i,j,k,3]/norm[i,j,k]" "tab_xyz[i,j,k,4]/norm[i,j,k]" "tab_xyz[i,j,k,5]/norm[i,j,k]" "tab_xyz[i,j,k,6]/norm[i,j,k]" "norm[i,j,k]/max#... we print the average vectors
		}
	    }
	}
    }
}
