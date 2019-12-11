program contactmap
   implicit none
   integer :: natoms=734
   integer :: nbeads=92
   double precision coor(3,92),coor2(3,1000), coor_C(3,92)
   integer i,j,k,count1,idmy1,id,jd
   integer atom_num(1000),res_seq(1000)
   integer xyz_atom_num(1000)
   integer, parameter :: dk = kind(0d0) 
   double precision :: cutoff=5.5d0,mindist=0.0d0,scaling=1.2
   character*10  cdmy1,cdmy2,cdmy3,cdmy4,cdmy5,cdmy6,cdmy7
   character*4 atom_name(1000),hydrogen(1000),ai,aj
   double precision dist(1000,1000),dist2(1000,1000),a,b,c
   double precision dist_CA(1000,1000),x,y,z,x1,y1,z1,Q
   integer conmap(1000,1000)
   character*100 pdbfile,xyzfile
   character*200 string,string1
   integer ios,diff,pdb_contacts,xyz_contacts
   integer :: maxlines=1000
    
!  store native contacts formed in native pdb as true in contact_n
!  store contacts formed after comparison with pdb as true in contact_u
   logical :: fractional,contact_n(92,92),contact_u(92,92) 
   CALL getarg(1,pdbfile)
   CALL getarg(2,xyzfile)


! list atomnames to be excluded from distance calc. from pdb here
    hydrogen(1)='   H'
    hydrogen(2)='  OH'
    hydrogen(3)='  NH'
    hydrogen(4)=' NH1'
    hydrogen(5)=' NH2'
!  read pdb and xyz files
   open (unit=7, file=trim(pdbfile),status='old')
   open (unit=8, file=trim(xyzfile),status='old')
!  write files
   open (unit=9, file='pdb_contact',status='unknown')
   open (unit=10, file='xyz_contact',status='unknown')
  ! open (unit=11, file='CA_dist_table',status='unknown')
  ! open (unit=12, file='out',status='unknown')
!  intialise
   do i=1,1000   
    do j=1,1000
      dist(i,j)=0.0D0
      dist2(i,j)=0.0d0
      atom_num(i)=0
    enddo
   enddo

!    
  do i=1,nbeads
   do j=1,nbeads
     contact_n(i,j)=.false.
     contact_u(i,j)=.false.
     dist_CA(i,j)= 0.0D0
   enddo
  enddo
!  generate contact map from all atom pdb 
   j=0  
   k=0
   write(*,*) 'Read pdbfile : ',pdbfile
   do i=1,maxlines
       read ( 7, '(a)', iostat = ios ) string
         if ( string(1:4) .eq. 'ATOM' ) then
          j=j+1
          cdmy1 = (trim(string(31:38)))
          read (cdmy1, *) coor2(1,j)
!          write(*,*) coor2(1,j),j
          cdmy2 = (trim(string(39:46)))
          read (cdmy2, *) coor2(2,j)
          cdmy3 = (trim(string(47:54)))
          read (cdmy3, *) coor2(3,j)
          cdmy4 = (trim(string(7:11))) 
          read (cdmy4, *) atom_num(j)
          cdmy5 = (trim(string(13:16)))
          read (cdmy5, *) atom_name(j)
          !write(*,*) atom_name(j)          
          cdmy6 = (trim(string(23:26)))
          read (cdmy6, *) res_seq(j)

!       store all alpha carbon coordinates.
          if (atom_name(j).eq.'CA') then
           k=k+1
           coor_C(1,k)=coor2(1,j)
           coor_C(2,k)=coor2(2,j)
           coor_C(3,k)=coor2(3,j)
          endif




!       residue number is not zero
          if (res_seq(j).eq.0) then
          write(*,*) "READ>>Residue number can't be zero. STOP.",res_seq(j)
          STOP
          endif
         endif

     enddo


   write(*,*) 'Number of atoms read from pdb = ',j
   write(*,*) 'Number of CA read from pdb = ',k
   natoms=j
   nbeads=k
   if (natoms.gt.1000) then
    write(*,*) 'Too many atoms. STOP. Set less than 1000. Change array sizes.'
    STOP
   endif

! make calpha dist matrix
   do i=1,nbeads
    do j=1,nbeads
      x=(coor_C(1,i)-coor_C(1,j))
      y=(coor_C(2,i)-coor_C(2,j))
      z=(coor_C(3,i)-coor_C(3,j))
      dist_CA(i,j) = sqrt( (x**2) + (y**2) + (z**2))
!      write(12,*) coor_C(1,i)
      if (isnan(dist_CA(i,j))) then
      write(*,*) 'NaN distance encountered!, STOP',coor_C(1,i),coor_C(1,j),&
    &coor_C(2,i),coor_C(2,j),coor_C(3,i),coor_C(3,j),i,j,dist_CA(i,j)
      STOP
      endif
!      write(11,*) i,j,dist_CA(i,j)
    enddo
   enddo
   
   j=0
   k=0
   count1=0

!    calculate CA-CA distances for all residues
!     count1=0
!     do i=1,734
!      do j=1,734
!      id=res_seq(i)
!      jd=res_seq(j)
!      ai=trim(atom_name(i))                                                                                  
!      aj=trim(atom_name(j)) 
!      if (ai.eq.' CA'.and.aj.eq.'CA') then
!         dist_CA(id,jd)=dist2(i,j)
!         write(11,*) ai,aj,id,jd,dist_CA(id,jd),dist2(i,j)
!         count1=count1+1
!         endif
!      enddo
!     enddo
!      if (count1.ne.nbeads) then 
!        write(*,*) 'CA atoms different from no. of beads', count1,nbeads
!        stop
!      endif



    k=0
!  atom-atom distances, cut-offs, etc.
   do i=1,natoms
    do j = 1,natoms
       x1=(coor2(1,i)-coor2(1,j))   
       y1=(coor2(2,i)-coor2(2,j))   
       z1=(coor2(3,i)-coor2(3,j))   
      dist2(i,j)=sqrt((x1**2) +(y1**2)+(z1**2))
      id=res_seq(i)
      jd=res_seq(j)
      ai=trim(atom_name(i))
      aj=trim(atom_name(j)) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1     
          if ((id.eq.0)) then
           write(*,*) "ID>>Residue sequence can't be zero. STOP."             
           write(*,*) atom_num(i),atom_name(i),res_seq(i),i,j
           STOP
          endif
          if ((jd.eq.0)) then
           write(*,*) "JD>>:Residue sequence can't be zero. STOP."
           write(*,*) atom_num(i),atom_name(j),res_seq(j),i,j
           STOP
          endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       apply cutoff criterion
        if ((dist2(i,j).lt.cutoff).and.(jd>id+3)) then
!          if (ai.eq.'CA'.and.aj.eq.'CA'.and.dist_CA(id,jd).ne.(dist2(i,j)))then  
!          write (*,*) 'Something wrong with reading CAs. STOP'
!         write (12,*) 'coor_C = ',coor_C(1,id),coor_C(2,id),coor_C(3,id),coor_C(1,jd),coor_C(2,jd),coor_C(3,jd)
!          write (12,*) 'coor2  = ',coor2(1,i),coor2(2,i),coor2(3,i),coor2(1,j),coor2(2,j),coor2(3,j)
!!          stop
!          endif
          contact_n(id,jd)=.true.
          k=k+1
!          write(9,*) contact_n(id,jd),i,j,id,jd,ai,aj,dist_CA(id,jd),dist2(i,j)        
        endif


   enddo
  enddo 

  write(*,*) 'Number of contacts with duplicates = ', k
  k=0
! count contacts without any duplicates and write to file.
! all contacts are stored in logical contact_n(i,j)
   count1=0
   do i=1,nbeads
     do j=1,nbeads
     if (contact_n(i,j).and.(j.gt.i+3)) then
       write(9,*) i,j,dist_CA(i,j)
        count1=count1+1
     endif     
     enddo
   enddo
   pdb_contacts=count1
   write(*,*) 'Number of contacts without duplicates = ', pdb_contacts 

!   read grofile
   write(*,*) 'Reading xyz file :', xyzfile
    read(8,*) cdmy1
    read(8,*) cdmy1
   do i=1,nbeads
      read(8,*)cdmy1,coor(1,i),coor(2,i),coor(3,i)
!      read(8,*)coor(1,i),coor(2,i),coor(3,i)
!     write(*,*) cdmy1,cdmy1,xyz_atom_num(i)
   enddo
   i=0
   j=0
   k=0
   count1=0
!    do i = 1,nbeads
!!  nm to ang!
!       coor(1,i)=10*coor(1,i)
!       coor(2,i)=10*coor(2,i)
!       coor(3,i)=10*coor(3,i)
!      enddo



!  xyz file here.
   count1=0
   do i=1,nbeads
    do j=1,nbeads
       x=(coor(1,i)-coor(1,j))   
       y=(coor(2,i)-coor(2,j))   
       z=(coor(3,i)-coor(3,j))          
      
     dist(i,j)=sqrt((x**2) + (y**2) + (z**2))
        if ((dist(i,j).le.scaling*dist_CA(i,j).and.(contact_n(i,j)))) then
            count1=count1+1
            contact_u(i,j)=.true.
            write(10,*)i,j,dist(i,j),dist_CA(i,j)
        endif
    enddo
   enddo
   xyz_contacts=count1
   Q=(real(xyz_contacts,dk)/(pdb_contacts))
   write(*,*)'Number of contacts mapped in gro file to pdb file =',xyz_contacts
   write(*,*)'Q=',Q,pdb_contacts,xyz_contacts
 
   close(7)
   close(8)
   close(10)
   close(9)
end

