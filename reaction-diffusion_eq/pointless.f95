      program pointless
              implicit none
              real, pointer :: p
              real,target :: a= 3.1, b=4.7

              p=>a
              print*,p
              p=4.0
              print*,a
              p=>b
              print*,b
              print*,associated (p)
              nullify(p)
              print*,associated(p)

       end program pointless
