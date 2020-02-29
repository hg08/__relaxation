#Program to treat the data which are obtained from orientation.f95
#The output is a data file where the angle between O-H and the normal axis is a function fo the distance of the interface
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
	tab_ref[i,3]=0
	tab_ref[i,4]=0
	tab_ref[i,5]=0
	tab_ref[i,6]=0
	norm[i]
    }

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
		tab_ref[int(tab1[i]*ratio),1]+=tab2[i]#We sum the vectors u1
		tab_ref[int(tab1[i]*ratio),2]+=tab3[i]#We sum the vectors v1
		tab_ref[int(tab1[i]*ratio),3]+=tab4[i]#We sum the vectors w1
		tab_ref[int(tab1[i]*ratio),4]+=tab5[i]#We sum the vectors u2
		tab_ref[int(tab1[i]*ratio),5]+=tab6[i]#We sum the vectors v2
		tab_ref[int(tab1[i]*ratio),6]+=$i     #We sum the vectors w2
		norm[int(tab1[i]*ratio)]   +=sqrt((tab3[i]+tab6[i]+$i)**2)#Sum of the norm of (w1+w2)
	    }
	}
    }
    else if(NR==2){
	nb_at=NF-1
    }
}


END{
    print "#Projection_axis(Ang)   x(ref-L1)(Ang)   y(ref-L1)(Ang)   z(ref-L1)(Ang)   x(ref-L2)(Ang)   y(ref-L2)(Ang)   z(ref-L2)(Ang)   Reliability(%)"
    #The reliability is equals to 100% if all the dipolar moments -- for a given z -- have the average orientation, so if all the water molecules with a given z are frozen with a specifc angle with the normal axis
    #The reliability is equals to 0% if for a given z is randomly oriented, therefore the average angle is meaningless.
    for(i=0;i<=int(param/div)+1;i++){
	if(norm[i]!=0){#If there is at least one atom in this area...
	    print i/ratio" "tab_ref[i,1]" "tab_ref[i,2]" "tab_ref[i,3]" "tab_ref[i,4]" "tab_ref[i,5]" "tab_ref[i,6]" "100*sqrt((tab_ref[i,3]+tab_ref[i,6])**2)/norm[i]#... we print the average vectors and the reliability
	}
    }
}
