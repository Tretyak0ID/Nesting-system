!Module for fundamental constants
module const_mod
implicit none
real(kind=8), parameter :: pi = acos(-1._8)

real(kind=8), parameter :: Earth_grav = 9.80616_8   !Earth gravity acceleration m/s^2
real(kind=8), parameter :: Earth_radii = 6371.22_8 * 1000.0_8 !Earth radius
real(kind=8), parameter :: Earth_sidereal_T = 23*60*60+56*60+4!Earth sidereal period
real(kind=8), parameter :: Earth_omega = 2.0_8*pi/Earth_sidereal_T!Earth angular speed

real(kind=8), parameter :: Day24h_sec = 24._8*3600._8 !day length in seconds

real(kind=8), parameter :: rgaz = 0.2870597E+03_8 !dry air gas-constant
real(kind=8), parameter :: Cp   = 3.5_8*rgaz
real(kind=8), parameter :: Cv   = 2.5_8*rgaz
real(kind=8), parameter :: kappa = 2.0_8 / 7.0_8

real(kind=8), parameter :: pcori = 2.0_8 * 7.292e-5_8

contains

end module const_mod
