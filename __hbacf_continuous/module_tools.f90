      MODULE tools 
      !2020/2/15
      !===================================================
      ! The functions 
      !===================================================      
      !=======================================
      implicit none

      integer,parameter :: rk=4              
      !real(kind=rk)     :: direct_distance2, distance2, diff_axis 

      REAL FUNCTION direct_distance2(u1,v1,w1,u2,v2,w2)
          real(kind=rk), INTENT(IN) :: u1,v1,w1,u2,v2,w2)
          direct_distance2 = (u2-u1)**2 + (v2-v1)**2 + (w2-w1)**2
      END FUNCTION direct_distance2

      REAL FUNCTION distance2(u1,v1,w1,u2,v2,w2,a,b,c)
          real(kind=rk),INTENT(IN) :: u1,v1,w1,u2,v2,w2,a,b,c
          logical :: A1,A2,A3,B1,B2,B3,C1,C2,C3
          A1 = (abs(u1-u2) > a/2 .AND. u1 > u2)
          A2 = (abs(u1-u2) > a/2 .AND. u2 > u1)
          A3 = (abs(u1-u2) < a/2)
          B1 = (abs(v1-v2) > b/2 .AND. v1 > v2)
          B2 = (abs(v1-v2) > b/2 .AND. v2 > v1)
          B3 = (abs(v1-v2) < b/2)
          C1 = (abs(w1-w2) > c/2 .AND. w1 > w2)
          C2 = (abs(w1-w2) > c/2 .AND. w2 > w1)
          C3 = (abs(w1-w2) < c/2)
          if (A3 .and. B3 .and. C3) then
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B3 .and. C3) then
              u2 = u2-a
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B3 .and. C3) then
              u1 = u1-a
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B1 .and. C3) then
              v1 = v1 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B3 .and. C2) then
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B2 .and. C3) then
              v2 = v2 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B3 .and. C1) then
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B1 .and. C3) then
              u1 = u1 - a
              v1 = v1 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B2 .and. C3) then
              u1 = u1 - a
              v2 = v2 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B3 .and. C1) then
              u1 = u1 - a
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B3 .and. C2) then
              u1 = u1 - a
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B1 .and. C3) then
              u2 = u2 - a
              v1 = v1 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B2 .and. C3) then
              u2 = u2 - a
              v2 = v2 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B3 .and. C1) then
              u2 = u2 - a
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B3 .and. C2) then
              u2 = u2 - a
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B1 .and. C1) then
              v1 = v1 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B1 .and. C2) then
              v1 = v1 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B2 .and. C1) then
              v2 = v2 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B2 .and. C2) then
              v2 = v2 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B1 .and. C1) then
              u1 = u1 - a
              v1 = v1 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B1 .and. C2) then
              u1 = u1 - a
              v1 = v1 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B2 .and. C1) then
              u1 = u1 - a
              v2 = v2 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B2 .and. C2) then
              u1 = u1 - a
              v2 = v2 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B1 .and. C1) then
              u2 = u2 - a
              v1 = v1 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B1 .and. C2) then
              u2 = u2 - a
              v1 = v1 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B2 .and. C1) then
              u2 = u2 - a
              v2 = v2 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B2 .and. C2) then
              u2 = u2 - a
              v2 = v2 - b
              w2 = w2 - c
              distance2 = direct_distance2(u1,v1,w1,u2,v2,w2) 
      END FUNCTION distance2

      REAL FUNCTION diff_axis(u1,u2,a)
          logical :: A1,A2,A3
          real(kind=rk),INTENT(IN) :: u1, u2, a
          A1 = (abs(u1-u2) > a/2 .AND. u1 > u2)
          A2 = (abs(u1-u2) > a/2 .AND. u2 > u1)
          A3 = (abs(u1-u2) < a/2)
          if (A3) then
              diff_axis= u1 - u2
          elseif A1 then
              u1 = u1 - a
              diff_axis= u1-u2
          elseif A2 then
              u2 = u2 - a
              diff_axis= u1 - u2
      END FUNCTION diff_axis

      END MODULE 
