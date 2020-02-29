#Program to treat the data from lifetime.f95

function statistics(inc,tab,start)
{
    if(inc>0){
	coef2=(inc+1)%2   ; div2=(inc+1-coef2)/2    ; coef2=coef2/2
	coef4=(inc+3)%4   ; div4=(inc+3-coef4)/4    ; coef4=coef4/4
	coef10=(inc+9)%10 ; div10=(inc+9-coef10)/10 ; coef10=coef10/10

	average=0 ; 
	for(i=1+start;i<=inc+start;i++){
	    average+=tab[i]
	}

	printf("%10.5f ",tab[start+1])                                                  #Minimum
	printf("%10.5f ",(1-coef10)*tab[div10+start]+coef10*tab[div10+1+start])         #First decile
	printf("%10.5f ",(1-coef4)*tab[div4+start]+coef4*tab[div4+1+start])             #First quartile
	printf("%10.5f ",(1-coef2)*tab[div2+start]+coef2*tab[div2+1+start])             #Median
	printf("%10.5f ",coef4*tab[inc-div4+start]+(1-coef4)*tab[inc-div4+1+start])     #Third quartile
	printf("%10.5f ",coef10*tab[inc-div10+start]+(1-coef10)*tab[inc-div10+1+start]) #Nineth decile
	printf("%10.5f ",tab[inc+start])                                                #Maximum
	printf("%10.5f\t",average/inc)
    }
    else{
	printf("%10s %10s %10s %10s %10s %10s %10s %10s\t","x","x","x","x","x","x","x","x")
    }
}

BEGIN{
    print "#Lifetime init (step) / Lifetime end (step) / Occurence of the lifetime / full system{ z-minimum (Ang) / 1st decile (Ang) / 1st quartile / z-median (Ang) / 3rd quartile (Ang) / 9th decile (Ang) / z-maximum (Ang) / z-average (Ang)} / Occurence on the first half (max over both half and over time=1) / First half{idem} / Occurence on the second half / Second half{ idem } / Symmetry applied {idem}" 
}

{
    if(NR>4){
	if(NR%2==0){#Odd lines
	    #Data are stored
	    inc2=NF-4
	    for(i=5;i<=NF;i++){
		tab[inc1+i-4]=$i
	    }

	    #Data are writen
	    printf("%i\t%i\t%10.5f ",$1,$2,$3)#General information
	    statistics(inc1+inc2,tab,0)#Data without any modification
	    printf("%10.5f ",occurence)
	    statistics(inc1,tab,0)#First side of the data
	    printf("%10.5f ",$4)
	    statistics(inc2,tab,inc1)#Second side of the data

	    #Symmetry
	    i=0;j=0;k=0
	    while(j!=inc1&&k!=inc2){
		i++
		if(tab[1+j]<2*center-tab[inc1+inc2-k]){
		    order[i]=tab[1+j]
		    j++
		}
		else{
		    order[i]=2*center-tab[inc1+inc2-k]
		    k++
		}
	    }

	    if(j==inc1){
		for(j=k+1;j<=inc2;j++){
		    order[i+j-k]=2*center-tab[inc1+inc2-j+1]
		}
	    }
	    else{
		for(k=j+1;k<=inc1;k++){
		    order[i+k-j]=tab[k]
		}
	    }


	    statistics(inc1+inc2,order,0)#Symmetry
	    printf("\n")

	}

	else{#Even lines
	    #Data are only stored
	    occurence=$4
	    inc1=NF-4
	    for(i=5;i<=NF;i++){
		tab[i-4]=$i
	    }
	}

    }
    else if(NR==3){center=$3}
}
