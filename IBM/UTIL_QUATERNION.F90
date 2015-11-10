subroutine quaternion_form(axis,angle,quat)
use operation
implicit none

real*8,intent(in) :: axis(3),angle
real*8,intent(out) :: quat(4)

real*8 :: a(3),b

a=axis/mo(axis,3)
a=a*sin(0.5d0*angle)
b=cos(0.5*angle)
quat(1)=b
quat(2:4)=a

end subroutine quaternion_form
!--------------------------------------------------------
!
!--------------------------------------------------------
subroutine quaternion_rotation(vec,quat)
use operation
implicit none

real*8,intent(in) :: quat(4)
real*8,intent(inout) :: vec(3)

real*8 :: invQuat(4),qVec(4),temp(4)

invQuat=(/quat(1),-quat(2),-quat(3),-quat(4)/)
qVec(1)=0.0d0
qVec(2:4)=vec

temp=qm(quat,qVec)
temp=qm(temp,invQuat)
vec=temp(2:4)


end subroutine quaternion_rotation