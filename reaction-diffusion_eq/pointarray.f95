      program pointarray
              implicit none
              integer,parameter :: n=5
              real, target :: T(n,n)
              real,pointer :: bot_boundary(:),&
                      top_boundary(:), center4(:,:)
              top_boundary => T(:,1)
              bot_boundary => T(:,n)
              center4 => T(n/2:n/2, n/2:n/2+1)
              
              call random_number(T)
              top_boundary=0.0
              bot_boundary=1.0
              center4=3.0

              print*,T
        end program pointarray
