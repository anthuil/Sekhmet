program Sekhmet
  ! Release 1 (20251022), from Version 6.6
  ! Antoine Thuillier (antoine.thuillier@uliege.be)
  !
  ! This code models the evolution of the orbital parameters of a substellar body while its host star ascends the Red-Giant Branch.
  !
  ! Instructions :
  !    Enter the input parameters in the 'input.par' file
  !    Run Sekhmet
  !    Two output files are produced :
  !       - "Sekhmet.log" that records important informations during the runs
  !       - "Sekhmet_output_[...]" that contains the resulting data. The output file's name reflects the initial parameters set
  !         for the run when different from the reference scenario (detailed in "input.par").
  !    Run Amon-Ra (with Python3) to extract graphs from the data
  !
  ! Compilation (fast)
  ! gfortran -o Sekhmet.x Sekhmet_Release-1.f90 -fopenmp -x f95-cpp-input -O3 -ffast-math -march=native -funroll-loops -fno-protect-parens -ftree-vectorize -finline-limit=1500 -finline-functions-called-once -finit-local-zero -flto -ffree-form -fdefault-real-8 -fdefault-double-8 -finit-local-zero -I .
  ! Compilation (debug)
  ! gfortran -o Sekhmet.x Sekhmet_Release-1.f90 -fopenmp -x f95-cpp-input -Wsurprising -Wall -fbounds-check -Wunderflow -fbacktrace -ffpe-trap=zero -finit-local-zero -g -ffree-form -fdefault-real-8 -fdefault-double-8 -finit-local-zero -I .
  ! Execution
  ! ./Sekhmet.x


  ! DECLARATION OF CONSTANTS AND VARIABLES
  implicit none
  integer :: error_check
  integer :: formalism
  integer, parameter :: Logfile = 14
  integer, parameter :: Special_print_file = 15
  integer :: Special_print = 0 ! Prints more data in a separate file

  ! Converters from cgs to SI units
  real, parameter :: unit_conv_dist = 100.
  real, parameter :: unit_conv_mass = 1000.
  real, parameter :: unit_conv_L = 10000000.

  real, parameter :: &
    pi        = 3.1415926535d0, &
    pim2      = 2.d0 * pi, &
    pim4      = 4.d0 * pi, &
    sec       = 3.1557807d7, & ! Number of seconds in a year.
    c         = 2.99792458d8 * unit_conv_dist, & ! [m/s]
    g         = 6.67430d-11 * unit_conv_dist**3 / unit_conv_mass, &
    boltz     = 1.380649d-23 * unit_conv_mass * unit_conv_dist**2, & ! [J K-1]
    Nav       = 6.0221367d23, &
    rk        = Nav * boltz, &      ! 8.3145112d7
    M_Sun     = 1.9884d30 * unit_conv_mass, & ! [kg]
    R_Sun      = 6.957d8 * unit_conv_dist, & ! [m]
    L_Sun      = 3.8275d26 * unit_conv_L, &
    au        = 1.49597870700d11 * unit_conv_dist, & ! [m]
    M_Earth    = 5.9722e24 * unit_conv_mass, & ! [kg]
    M_Jupiter = 1.8986e27 * unit_conv_mass, & ! [kg]
    R_Jupiter  = 6.9911e7 * unit_conv_dist, & ! [m]
    M_Sun_year = M_Sun / sec ! [kg/s]

  real, parameter :: Cd = 0.9 ! Dimensionless drag coefficient for a sphere.
  real, parameter :: cF = 0.9 ! Frequency parameter

  real :: L_star      ! Star's luminosity [L_Sun]
  real :: R_eff_star  ! Star's photospheric radius [R_Sun]
  real :: T_eff       ! Star's effective temperature [K]
  real :: dot_M_star  ! Star's mass loss [M_Sun.yr-1]
  real :: M_star      ! Star's mass [M_Sun]
  real :: dt_star     ! Timestep of STAREVOL [yr]
  real :: age_star    ! Star's age [yr]
  real :: sub_step    ! Coefficient used for interpolation.
  real :: M_env       ! Mass of the envelope of the star [kg].
  real :: R_env       ! Radius of the base of the envelope of the star [m].
  real :: R_star      ! Radius of the star [m].
  real :: k2_star     ! Love number for the second-order harmonic potential of the star.
  real :: I_star      ! Moment of inertia [kg m2].
  real :: dot_I_star  ! Variation of the moment of inertia [kg m2 s-1].
  real :: J_star      ! Angular momentum [kg m^2 s-1].
  real :: dot_J_star  ! Stellar torque [kg m^2 s-2].
  real :: omega_star  ! Surface rotation of the star [rad s-1].
  real :: dot_omega_star ! Variation of the self-rotation speed of the star [rad s-2]

  integer, allocatable :: model_array(:)   ! Stellar model number
  real, allocatable :: L_star_array(:)     ! Star's luminosity [L_Sun]
  real, allocatable :: R_eff_star_array(:) ! Star's photospheric radius [R_Sun]
  real, allocatable :: T_eff_array(:)      ! Star's effective temperature [K]
  real, allocatable :: M_star_array(:)     ! Star's mass [M_Sun]
  real, allocatable :: dt_star_array(:)    ! timestep of STAREVOL [yr].
  real, allocatable :: age_star_array(:)   ! Star's age [yr]
  real, allocatable :: Menv_base_array(:)  ! "Mass coordinate of the base of the convective envelope, in [M_Sun]. The envelope mass is simply M-Menv_base" according to STAREVOL documentation.
  real, allocatable :: Renv_base_array(:)  ! Radius of the base of the convective envelope, in [R_Sun].
  real, allocatable :: k2_star_array(:)    ! "Love number for the second-order harmonic potential" of the star, provided by Starevol.
  real, allocatable :: M_env_array(:)      ! Star's envelope mass [kg]
  real, allocatable :: I_star_array(:)     ! Moment of inertia of the star [kg m2]

  ! Planet's parameters
  real :: a         ! Planet's semi-major axis [m]
  real :: q         ! Mass ratio : M_planet / M_star.
  real :: v_orb     ! Orbital speed [m/s]
  real :: P_orb     ! Orbital period [s]
  real :: rho       ! Density of the medium [kg/m3] (or [g/m3] in cgs)
  real :: f_drag    ! Frictional drag
  real :: f_grav    ! Gravitational drag
  real :: M_planet  ! Planet's mass [kg]
  real :: R_planet  ! Planet's radius [m]
  real :: mu        ! Mean molecular weigth of the environment.
  real :: Q_p       ! Tidal dissipation parameter of the planet.
  real :: T_wind    ! Temperature of the medium [K].
  real :: c_s       ! Sound speed [m/s] (or [cm/s] in cgs).
  real :: tau       ! Eddy turnover timescale in the stellar enveloppe.
  real :: dot_a     ! Planet's semi-major axis variation [m/s].
  real :: dot_a_1   ! Planet's semi-major axis variation from mass conservation, drag and stellar tides (using V14) [m/s].
  real :: dot_a_2   ! Planet's semi-major axis variation from tidal interactions (using V14) [m/s].
  real :: J_planet          ! Angular momentum of the planet [kg m2 s-1].
  real :: omega_planet      ! Rotation speed of the planet [rad/s].
  real :: dot_omega_planet  ! Variation of the rotation speed of the planet [rad s-2].
  real :: I_planet          ! Momentum of inertia of the planet [kg m2].
  real :: I_planet_coef     ! I/(MR^2) | >0 and <0.4

  real :: e           ! Planet's eccentricity [ø].
  real :: dot_e       ! Planet's eccentricity variation [ø/s].
  real :: dot_e_1     ! Planet's eccentricity variation from angular momentum (using V14) [ø/s].
  real :: dot_e_2     ! Planet's eccentricity variation from stellar tides (using V14) [ø/s].
  real :: dot_e_eq_tides_star   ! Planet's eccentricity variation from stellar tides (using F23) [ø/s].
  real :: dot_e_eq_tidal_planet ! Planet's eccentricity variation from planetary tides (using F23) [ø/s].
  real :: n           ! Planet's mean motion (= 2*pi / P_orb).
  real :: Omega_e     !
  real :: fe          !
  real :: alpha_p     ! Square of the radius of gyration (Matsumura 2010) [ø]
  real :: k2dt_planet ! Product of k2 (Love number for the second-order harmonic potential of the planet) and delta_t (time-lag of the planet).
  real :: rotation_planet !

  ! Other parameters
  real :: dt        ! Time interval between two iterations of the model [s].
  real :: time0     ! Age of the star when the computation should start [s].
  real :: dt_max    ! Number of years to model (between time0 and the last computation) [s].
  real :: time_end  ! Time of the last computation to do [s].
  real :: time      ! [s]
  real :: eps       ! Maximum relative variation in "a" and "e" for the timestep.

  real :: wind_speed      ! Speed of the stellar wind [m/s].
  real :: f_orb           ! Orbital frequency of the planet.
  real :: J_syst          ! Total angular momentum of the system [kg m2 s-1].
  real :: J_orb ! Angular momentum contained in the the orbital motion [kg m2 s-1].
  real :: dot_J_orb ! Variation of the angular momentum contained in the the orbital motion [kg m2 s-2].

  real :: Pe          ! Periapsis of the planet [m].
  real :: Ap          ! Apoapsis of the planet [m].
  real :: f1, f2, f3  ! Frequency components (see Villaver 2014 for details).
  real :: g2, g3, g4, g5 ! Eccentricity function (see Villaver 2014 for details).
  character(128) :: Stellar_model
  character(141) :: headline
  character(100) :: Output_file_name
  real :: dum
  integer :: iter, nlines, istart, io, i, k, len_str
  integer :: force_continue ! 0 : stop if dt too small | 1 : force a bigger value for dt if too small.
  integer :: verbose ! 0 : no written details | 1 : small details | 2 : full details
  character(10) :: OFN_Mlrp, OFN_Star_mass, OFN_a, OFN_Mp, OFN_Rp, OFN_e, OFN_formalism, OFN_Omega_p, OFN_k2dt_p, OFN_Coef_I_p
  character(10) :: OFN_T_wind, OFN_wind_speed, OFN_mu, OFN_Omega_s, OFN_t0, OFN_dt_max, OFN_dt_max_2
  character(10) :: OFN_eps, OFN_force_continue, OFN_Qp, OFN_alpha_p
    ! Strings for output file name.

  call read_parameters(Stellar_model, a, e, R_planet, M_planet, Q_p, omega_planet, alpha_p, dt_max, T_wind, time0, wind_speed, mu, &
    eps, force_continue, verbose, formalism, k2dt_planet, omega_star, I_planet_coef)

  ! ##############################################################################################################################
  ! ##############################################################################################################################

  dt = 100. * sec ! Initial timestep (will be adjusted by the code).
  time = time0
  time_end = time0 + dt_max
  c_s = sqrt(((5/3.) * Rk * T_wind) / mu)
    ! 5/3 : adiabatic coef
  I_planet = M_planet * R_planet**2 * I_planet_coef

  ! ##############################################################################################################################
  ! Operations

  if (Logfile /= 6) open(Logfile, file='./Sekhmet.log', status='unknown', access='sequential', &
    action='write', form='formatted', iostat=error_check)
  if (error_check .ne. 0) stop "Error while opening log file"

  if (Special_print == 1) then
    open(Special_print_file, file='./Special_print.csv', status='unknown', access='sequential', &
      action='write', form='formatted', iostat=error_check)
      write(Special_print_file, *) "[ø], [s], [m/s], [m/s], [m/s], [m/s], [m/s]"
      write(Special_print_file, *) "iter, time, dot_a, dot_a_drag, dot_a_eq_tides_star, dot_a_eq_tides_planet, dot_a_ang_mom"
  endif

  ! Reading the stellar model.
  write(Logfile,*) "Opening STAREVOL model file..."
  open(unit=12, file=trim(Stellar_model), status='old', form='formatted', iostat=error_check)
  if (error_check .ne. 0) then
    print *,'error =', error_check,' file = ',trim(Stellar_model)
    stop "Error while opening STAREVOL file"
  endif
  write(Logfile,*) "... Done"

  read(12, FMT='(141A)') headline ! Reading the first line (headline).
  !write(*,FMT='(141A)') headline

  ! Determining the number of lines in the input file (stellar model).
  nlines = 0
  do
    read(12, *, iostat=io)
    if (io /= 0) exit
    nlines = nlines + 1
  enddo
  write(Logfile,*) nlines, " lines detected in Starevol file."

  ! Looking for the closest iteration of STAREVOL to time0. //////////////////////////////////////////////////////////////////////
  write(Logfile,*) "Looking for stellar properties at time0..."

  allocate(model_array(nlines), L_star_array(nlines), R_eff_star_array(nlines), T_eff_array(nlines), M_star_array(nlines), &
    dt_star_array(nlines), age_star_array(nlines), I_star_array(nlines), k2_star_array(nlines), Menv_base_array(nlines), &
    Renv_base_array(nlines), M_env_array(nlines))

  if (verbose >= 2) then
    write(Logfile, '(2x, "model", 5x, "L_star", 9x, "R_star", 6x, "T_star", 6x, "M_star", 8x, "dt_star", 8x, "age_star", 8x, &
      & "M_env_base", 8x, "R_env_base", 8x, "I_star", 9x, "k2_star")')
  end if

  ! Load the data in a file.
  rewind (12)
  read(12, '(1x)')
  do i = 1, nlines
    read(12, *) model_array(i), L_star_array(i), R_eff_star_array(i), dum, T_eff_array(i), dum, M_star_array(i), &
      dt_star_array(i), age_star_array(i), I_star_array(i), k2_star_array(i), Menv_base_array(i), Renv_base_array(i)
    if (age_star_array(i) < time0 / sec) istart = i
    if (verbose >= 2) then
      write(Logfile, FMT='(I6, F15.8, F14.8, F10.2, F14.8, ES15.4, ES16.6, ES18.8, ES18.8, F15.8, F14.7)') &
        model_array(i), L_star_array(i), R_eff_star_array(i), T_eff_array(i), M_star_array(i), dt_star_array(i), &
        age_star_array(i), Menv_base_array(i), Renv_base_array(i), I_star_array(i), k2_star_array(i)
    end if
  enddo
  close(12)

  if (verbose >= 2) then
    write(Logfile, '(2x, "model", 5x, "L_star", 9x, "R_star", 6x, "T_star", 6x, "M_star", 8x, "dt_star", 8x, "age_star", 8x, &
      & "M_env_base", 8x, "R_env_base", 8x, "I_star", 9x, "k2_star")')
  end if

  ! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ! Initial interpolation on the stellar model to find the stellar parameters at time0.

  ! ------------------------------------------------------------------------------------------------------------------------------
  ! Conversion to SI units.
  L_star_array     = L_star_array * L_Sun     ! [W] from L_Sun
  R_eff_star_array = R_eff_star_array * R_Sun ! [m] from R_Sun
  M_star_array     = M_star_array * M_Sun     ! [kg] from M_Sun
  dt_star_array    = dt_star_array * sec      ! [s] from year
  age_star_array   = age_star_array * sec     ! [s] from year
  Menv_base_array  = Menv_base_array * M_Sun  ! [kg] from M_Sun
  Renv_base_array  = Renv_base_array * R_Sun  ! [m] from R_Sun
  M_env_array      = M_star_array - Menv_base_array  ! [kg]
  I_star_array     = I_star_array * M_Sun * R_Sun**2 ! [kg m2] from M_Sun R_sun**2

  sub_step = (time0 - age_star_array(istart)) / (age_star_array(istart+1) - age_star_array(istart)) ! [s / s]

  L_star     = L_star_array(istart)     + (sub_step * (L_star_array(istart+1)     - L_star_array(istart)))     ! [W]
  R_eff_star = R_eff_star_array(istart) + (sub_step * (R_eff_star_array(istart+1) - R_eff_star_array(istart))) ! [m]
  T_eff      = T_eff_array(istart)      + (sub_step * (T_eff_array(istart+1)      - T_eff_array(istart)))      ! [K]
  M_star     = M_star_array(istart)     + (sub_step * (M_star_array(istart+1)     - M_star_array(istart)))     ! [kg]
  dt_star    = dt_star_array(istart)    + (sub_step * (dt_star_array(istart+1)    - dt_star_array(istart)))    ! [s]
  age_star   = age_star_array(istart)   + (sub_step * (age_star_array(istart+1)   - age_star_array(istart)))   ! [s]
  R_env      = Renv_base_array(istart)  + (sub_step * (Renv_base_array(istart+1)  - Renv_base_array(istart)))  ! [m]
  M_env      = M_env_array(istart)      + (sub_step * (M_env_array(istart+1)      - M_env_array(istart)))      ! [kg]
  I_star     = I_star_array(istart)     + (sub_step * (I_star_array(istart+1)     - I_star_array(istart)))     ! [kg m2]

  J_star = I_star * omega_star
  J_orb = M_star * M_planet * sqrt(G * a * (1. - e**2) / M_star)

  ! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ! Creating an output file to write the results in.

  Output_file_name = "./Sekhmet_output"


  ! Looking for the mass of the star.
  if (index(Stellar_model, "1.0") /= 0) then
    OFN_Star_mass = "_Ms1.0"
  elseif (index(Stellar_model, "1.2") /= 0) then
    OFN_Star_mass = "_Ms1.2"
  elseif (index(Stellar_model, "1.5") /= 0) then
    if(verbose >= 1) then ! If verbose >= 1 print the full details in the name.
      OFN_Star_mass = "_Ms1.5"
    else
      OFN_Star_mass = ""
    endif
  elseif (index(Stellar_model, "1.7") /= 0) then
    OFN_Star_mass = "_Ms1.7"
  elseif (index(Stellar_model, "2.0") /= 0) then
    OFN_Star_mass = "_Ms2.0"
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_Star_mass)

  ! Looking for the mass loss rate prescription (Mlrp).
  if (index(Stellar_model, "reimers0.2") /= 0) then
    OFN_Mlrp = "_r0.2"
  elseif (index(Stellar_model, "reimers0.5") /= 0) then
    OFN_Mlrp = "_r0.5"
  elseif (index(Stellar_model, "schroeder") /= 0) then
    if(verbose >= 1) then ! If verbose >= 1 print the full details in the name.
      OFN_Mlrp = "_schr"
    else
      OFN_Mlrp = ""
    endif
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_Mlrp)

  ! Looking for the initial semi-major axis of the planet.
  if (a / au /= 2) then
    write (OFN_a, '(F4.2)') a / au
    OFN_a = "_a" // OFN_a
  elseif (verbose >= 1) then
    write (OFN_a, '(F4.2)') a / au
    OFN_a = "_a" // OFN_a
  else
    OFN_a = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_a)

  ! Looking for the mass of the planet.
  if (M_planet / M_Jupiter > 1.02 .or. M_planet / M_Jupiter < 0.98) then
    write (OFN_Mp, '(F4.2)') M_planet / M_Jupiter
    OFN_Mp = "_Mp" // OFN_Mp
  elseif (verbose >= 1) then
    write (OFN_Mp, '(F4.2)') M_planet / M_Jupiter
    OFN_Mp = "_Mp" // OFN_Mp
  else
    OFN_Mp = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_Mp)
  !write(*,*) "M_planet / M_Jupiter ", M_planet / M_Jupiter


  ! Looking for the radius of the planet.
  if (R_planet / R_Jupiter > 1.02 .or. R_planet / R_Jupiter < 0.98) then
    write (OFN_Rp, '(F4.2)') R_planet / R_Jupiter
    OFN_Rp = "_Rp" // OFN_Rp
  elseif (verbose >= 1) then
    write (OFN_Rp, '(F4.2)') R_planet / R_Jupiter
    OFN_Rp = "_Rp" // OFN_Rp
  else
    OFN_Rp = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_Rp)

  ! Looking for the initial eccentricity of the planet.
  if (e /= 0) then
    write (OFN_e, '(F4.2)') e
    OFN_e = "_e" // OFN_e
  elseif (verbose >= 1) then
    write (OFN_e, '(F4.2)') e
    OFN_e = "_e" // OFN_e
  else
    OFN_e = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_e)

  ! Looking for the formalism
  if (formalism == 1) then
    OFN_formalism = "_F23"
  elseif (formalism == 2) then
    OFN_formalism = "_V14"
  elseif (formalism == 3) then
    if(verbose >= 1) then ! If verbose >= 1 print the full details in the name.
      OFN_formalism = "_Hme"
    else
      OFN_formalism = ""
    endif
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_formalism)

  ! Looking for the initial eccentricity of the planet.
  if (Omega_planet /= 7e-5) then
    write (OFN_Omega_p, '(A3, F3.1, A1, I2)') "_Op", 10. * Omega_planet / (10.**int(log10(Omega_planet))), &
      "e", int(log10(Omega_planet)) - 1
  elseif (verbose >= 1) then
    write (OFN_Omega_p, '(A3, F3.1, A1, I2)') "_Op", 10. * Omega_planet / (10.**int(log10(Omega_planet))), &
      "e", int(log10(Omega_planet)) -1
  else
    OFN_Omega_p = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_Omega_p)
  !write (*, *) "_Op", Omega_planet / (10.**int(log10(Omega_planet))), "e", int(log10(Omega_planet))


  ! Looking for the k2dt of the planet.
  if (k2dt_planet /= 1e-2) then
    write (OFN_k2dt_p, '(A4, F3.1, A1, I2)') "_ktp", 10. * k2dt_planet / (10.**int(log10(k2dt_planet))), &
      "e", int(log10(k2dt_planet)) -1
  elseif (verbose >= 2) then
    write (OFN_k2dt_p, '(A4, F3.1, A1, I2)') "_ktp", 10. * k2dt_planet / (10.**int(log10(k2dt_planet))), &
      "e", int(log10(k2dt_planet)) -1
  else
    OFN_k2dt_p = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_k2dt_p)
  !write (*, *) "_ktp", k2dt_planet / (10.**int(log10(k2dt_planet))), "e", int(log10(k2dt_planet))

  ! Looking for the coefficient for the moment of inertia of the planet.
  if (I_planet_coef /= 0.3) then
    write (OFN_Coef_I_p, '(F4.2)') I_planet_coef
    OFN_Coef_I_p = "_Ip" // OFN_Coef_I_p
  elseif (verbose >= 1) then
    write (OFN_Coef_I_p, '(F4.2)') I_planet_coef
    OFN_Coef_I_p = "_Ip" // OFN_Coef_I_p
  else
    OFN_Coef_I_p = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_Coef_I_p)

  ! Looking for the temperature of the medium (wind).
  if (T_wind /= 2500.) then
    write (OFN_T_wind, '(I4)') int(T_wind)
    OFN_T_wind = "_Tw" // OFN_T_wind
  elseif (verbose >= 2) then
    write (OFN_T_wind, '(I4)') int(T_wind)
    OFN_T_wind = "_Tw" // OFN_T_wind
  else
    OFN_T_wind = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_T_wind)

  ! Looking for the speed of the wind.
  if (wind_speed / unit_conv_dist > 5001. .or. wind_speed / unit_conv_dist < 4999.) then
    write (OFN_wind_speed, '(I4)') int(wind_speed / unit_conv_dist)
    OFN_wind_speed = "_Ws" // OFN_wind_speed
  elseif (verbose >= 2) then
    write (OFN_wind_speed, '(I4)') int(wind_speed / unit_conv_dist)
    OFN_wind_speed = "_Ws" // OFN_wind_speed
  else
    OFN_wind_speed = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_wind_speed)

  ! Looking for mu.
  if (mu /= 1.38) then
    write (OFN_mu, '(F4.2)') mu
    OFN_mu = "_mu" // OFN_mu
  elseif (verbose >= 2) then
    write (OFN_mu, '(F4.2)') mu
    OFN_mu = "_mu" // OFN_mu
  else
    OFN_mu = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_mu)

  ! Looking for the initial rotation rate of the star.
  if (omega_star /= 1e-7) then
    write (OFN_Omega_s, '(A3, F3.1, A1, I2)') "_Os", 10. * omega_star / (10.**int(log10(omega_star))), &
      "e", int(log10(omega_star)) - 1
  elseif (verbose >= 1) then
    write (OFN_Omega_s, '(A3, F3.1, A1, I2)') "_Os", omega_star / (10.**int(log10(omega_star))), &
      "e", int(log10(omega_star))
  else
    OFN_Omega_s = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_Omega_s)

  ! Looking for the initial time.
  if (time0 / sec /= 2.56e9) then
    write (OFN_t0, '(A3, F5.3, A1, I1)') "_ti", (time0 / sec) / (10.**int(log10(time0 / sec))), "e", int(log10(time0 / sec))
  elseif (verbose >= 1) then
    write (OFN_t0, '(A3, F5.3, A1, I1)') "_ti", (time0 / sec) / (10.**int(log10(time0 / sec))), "e", int(log10(time0 / sec))
  else
    OFN_t0 = ""
  endif
  !write(*,*) "time0 ", time0
  Output_file_name = trim(Output_file_name) // trim(OFN_t0)

  ! Looking for the duration.
  if (dt_max / sec > 1.1e8 .or. dt_max / sec < 9.9e7) then
    write (OFN_dt_max, '(A3, F3.1, A1, I1)') "_dt", (dt_max / sec) / (10**int(log10(dt_max / sec))), "e", int(log10(dt_max / sec))
  elseif (verbose >= 1) then
    write (OFN_dt_max, '(A3, F3.1, A1, I1)') "_dt", (dt_max / sec) / (10**int(log10(dt_max / sec))), "e", int(log10(dt_max / sec))
  else
    OFN_dt_max = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_dt_max)

  ! Looking for the uncertainty coefficient.
  if (eps /= 0.01) then
    write (OFN_eps, '(A4, F3.1, A1, I2)') "_eps", 10. * eps / (10.**int(log10(eps))), "e", int(log10(eps)) - 1
  elseif (verbose >= 1) then
    write (OFN_eps, '(A4, F3.1, A1, I2)') "_eps", 10. * eps / (10.**int(log10(eps))), "e", int(log10(eps)) - 1
  else
    OFN_eps = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_eps)

  ! Looking for the end condition.
  if (force_continue == 0) then
    OFN_force_continue = "_fc0"
  elseif (force_continue == 1) then
    if(verbose >= 1) then ! If verbose >= 1 print the full details in the name.
      OFN_force_continue = "_fc1"
    else
      OFN_force_continue = ""
    endif
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_force_continue)

  ! Looking for the Q factor of the planet.
  if (formalism == 2 .and. Q_p /= 1.e5) then
    write (OFN_Qp, '(A3, F3.1, A1, I1)') "_Qp", Q_p / (10.**int(log10(Q_p))), "e", int(log10(Q_p))
  elseif (verbose >= 1) then
    write (OFN_Qp, '(A3, F3.1, A1, I1)') "_Qp", Q_p / (10.**int(log10(Q_p))), "e", int(log10(Q_p))
  else
    OFN_Qp = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_Qp)

  ! Looking for the alpha parameter of the planet.
  if (formalism == 2 .and. alpha_p /= 0.26) then
    write (OFN_alpha_p, '(F4.2)') alpha_p
    OFN_alpha_p = "_alp" // OFN_alpha_p
  elseif (verbose >= 1) then
    write (OFN_alpha_p, '(F4.2)') alpha_p
    OFN_alpha_p = "_alp" // OFN_alpha_p
  else
    OFN_alpha_p = ""
  endif
  Output_file_name = trim(Output_file_name) // trim(OFN_alpha_p)

  Output_file_name = trim(Output_file_name) // ".csv"
  !write(*,*) trim(Output_file_name)

  if (verbose >= 1) write(Logfile,*) Output_file_name

  open(13, file=Output_file_name, status='unknown', access='sequential', &
    action='write', form='formatted', iostat=error_check)
  if (error_check .ne. 0) stop "Error while opening file"

  write(13, FMT='(A38, F8.4, A22, F8.4, A22, F6.3, A21, F4.2, A28, ES8.2, A11, ES8.2, A15, F4.2, A19, ES8.2)') &
    "Planet's parameters | R_planet [Rj] : ", R_planet / R_Jupiter, "   |  M_planet [Mj] : ", M_planet / M_Jupiter, &
    "   | a initial [au] : ", a / au, "   | e initial [ø] : ", e, "   | omega_planet [rad/s] : ", omega_planet, "   | Q_p : ", &
    Q_p, "   | alpha_p : ", alpha_p, "   | k2dt_planet : ", k2dt_planet
  write(13, FMT='(A41, A50, A18, F6.0, A24, F6.0, A10, F4.2, A18, ES8.2)') "Environment parameters | Stellar model : ", &
    trim(Stellar_model), "   | T_wind [K] : ", T_wind, "   | wind_speed [m/s] : ", wind_speed / unit_conv_dist, "   | mu : ", &
    mu , "   | omega_star : ", omega_star
  write(13, FMT='(A67, ES9.3, A25, ES9.3, A34, ES8.2, A17, I1)') &
    "Programmatical parameters | Sekhmet version : 5.1   | time0 [yr] : ", time0 / sec, "   | max_duration [yr] : ", &
    dt_max / sec, "   | max. rel. variat. a and e : ", eps, "   | formalism : ", formalism

  if (verbose < 1) then ! If no details are asked.
    write(13,'("time", 17x, "M_star", 6x, "R_eff_star", 6x, "L_star", 6x, "dot_M_star", 7x, "M_env", 8x, "R_env", 12x, &
      & "a", 11x, "dot_a", 10x, "e", 9x, "omega_planet", 3x, "omega_star")')
    write(13,'("[yr]", 17x, "[M_Sun]", 6x, "[R_Sun]" ,8x, "[L_Sun]", 5x, "[M_Sun/yr]", 6x, "[M_Sun]", 7x, "[R_Sun]", 8x, "[au]", &
      & 9x, "[au/yr]", 8x, "[ø]", 10x, "[rad/s]", 8x, "[rad/s]")')
  else if (verbose < 2) then ! If it is asked to print more details.
    write(13,'("time", 17x, "M_star", 6x, "R_eff_star", 6x, "L_star", 6x, "dot_M_star", 7x, "M_env", 8x, "R_env", 12x, &
      & "a", 11x, "dot_a", 10x, "e", 9x, "omega_planet", 3x, "omega_star", 2x, "dot_omega_star", 4x, "J_orb", 7x, "dot_J_orb", &
      & 6x, "J_star", 6x, "dot_J_star", 6x, "I_star")')
    write(13,'("[yr]", 17x, "[M_Sun]", 6x, "[R_Sun]" ,8x, "[L_Sun]", 5x, "[M_Sun/yr]", 6x, "[M_Sun]", 7x, "[R_Sun]", 8x, "[au]", &
      & 9x, "[au/yr]", 8x, "[ø]", 10x, "[rad/s]", 8x, "[rad/s]")')
  else ! If it is asked to print even more details.
    write(13,'("time", 17x, "M_star", 6x, "R_eff_star", 6x, "L_star", 6x, "dot_M_star", 7x, "M_env", 8x, "R_env", 12x, &
      & "a", 11x, "dot_a", 10x, "e", 9x, "omega_planet", 3x, "omega_star", 2x, "dot_omega_star", 4x, "J_orb", 7x, "dot_J_orb", &
      & 6x, "J_star", 6x, "dot_J_star", 7x, "I_star", 7x, "J_planet", 7x, "J_syst", 5x, "dot_omega_planet")')
    write(13,'("[yr]", 17x, "[M_Sun]", 6x, "[R_Sun]" ,8x, "[L_Sun]", 5x, "[M_Sun/yr]", 6x, "[M_Sun]", 7x, "[R_Sun]", 8x, "[au]", &
      & 9x, "[au/yr]", 8x, "[ø]", 10x, "[rad/s]", 8x, "[rad/s]", 5x, "[rad/s2]", 5x, "[kg m2 s-1]", 3x, "[kg m2 s-2]", &
      & 3x, "[kg m2 s-1]", 3x, "[kg m2 s-2]", 5x, "[kg m2]", 5x, "[kg m2 s-1]", 3x, "[kg m2 s-1]", 5x, "[rad/s2]")')
  end if

! ################################################################################################################################
!  MAIN LOOP
! ################################################################################################################################

  iter = 0
  i = 1

  do while (time < time_end .and. i < nlines)

    iter = iter + 1

    do k = i, nlines
      if (age_star_array(k) > time) then
        i = k
        exit
      endif
    enddo

    ! Checking if age is correctly set.
    if (age_star_array(i-1) > time .or. age_star_array(i)< time) then
      print *, i, age_star_array(i-1) - time, age_star_array(i) - time
      stop '***** Error finding age *****'
    endif

    ! Interpolation of the stellar model.
    sub_step = (time - age_star_array(i-1)) / (age_star_array(i) - age_star_array(i-1))

    L_star     = L_star_array(i-1)     + (sub_step * (L_star_array(i)     - L_star_array(i-1)))
    R_eff_star = R_eff_star_array(i-1) + (sub_step * (R_eff_star_array(i) - R_eff_star_array(i-1)))
    R_star     = R_eff_star_array(i-1) + (sub_step * (R_eff_star_array(i) - R_eff_star_array(i-1))) ! Using R_eff as the stellar radius.
    T_eff      = T_eff_array(i-1)      + (sub_step * (T_eff_array(i)      - T_eff_array(i-1)))
    M_star     = M_star_array(i-1)     + (sub_step * (M_star_array(i)     - M_star_array(i-1)))
    dt_star    = dt_star_array(i-1)    + (sub_step * (dt_star_array(i)    - dt_star_array(i-1)))
    age_star   = age_star_array(i-1)   + (sub_step * (age_star_array(i)   - age_star_array(i-1)))
    k2_star    = k2_star_array(i-1)    + (sub_step * (k2_star_array(i)    - k2_star_array(i-1)))
    R_env      = Renv_base_array(i-1)  + (sub_step * (Renv_base_array(i)  - Renv_base_array(i-1)))
    M_env      = M_env_array(i-1)      + (sub_step * (M_env_array(i)      - M_env_array(i-1)))
    I_star     = I_star_array(i-1)     + (sub_step * (I_star_array(i)     - I_star_array(i-1)))

    dot_M_star = (M_star_array(i) - M_star_array(i-1)) / (age_star_array(i) - age_star_array(i-1))
    dot_I_star = (I_star_array(i) - I_star_array(i-1)) / (age_star_array(i) - age_star_array(i-1))

    if (M_star .lt. 0.) then ! If the mass of the star reaches 0 for any reason, we stop the computations.
      write(Logfile,*) "Star's mass reached 0, end of the computations."
      !write(Logfile, FMT='(I10, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5)') &
      !  time/sec, M_star, R_eff_star, L_star, dot_M_star, M_env, R_env, a, dot_a
      write(Logfile, FMT='(ES14.5, ES14.5, ES14.5, ES14.5)') M_star_array(i-1), sub_step, M_star_array(i), &
        (M_star_array(i) - M_star_array(i-1))
      exit
    endif

    if (a .gt. 10. * au) then ! If "a" reaches 10 AU, we consider the planet as ejected and we stop the computations.
      write(Logfile,*) "The planet has been ejected, end of the computations."
      !write(14, FMT='(I10, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5)') &
      !  time/sec, M_star, R_eff_star, L_star, dot_M_star, M_env, R_env, a, dot_a
      !write(14, FMT='(ES14.5, ES14.5, ES14.5, ES14.5)') M_star_array(i-1), sub_step, M_star_array(i), (M_star_array(i) - M_star_array(i-1))
      write(13, *) "# Survival = 3   contact time : ", -1

      exit
    endif

    !write(*, FMT='(A, i5, A, ES16.10, A,es16.10)') "Star age (iter) ", iter, "  (time) : ", time / sec , " years, dt = ", dt / sec
    write(*, FMT='(A, i5, A, ES16.10, A, ES16.10, A, F11.8)') "Iter ", iter, "  Star age (time) : ", time / sec , &
      " years, a = ", a / au, " AU, e = ", e

    write(Logfile,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    write(Logfile,*) "age_star = ", age_star / sec, " years"
    if(verbose >=2) then
      write(Logfile,*) "sub_step = ", sub_step
      write(Logfile,*) "M_star = ", M_star, " kg"
      write(Logfile,*) "R_eff_Star = ", R_eff_Star, " m"
      write(Logfile,*) "M_env = ", M_env
      write(Logfile,*) "R_env = ", R_env
      write(Logfile,*) "L_star = ", L_star
    end if

    ! ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ! Compute equation terms
    ! ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (formalism == 1) then ! Fellay's formalism.
      call Compute_F23 (dot_a, dot_e, G, pi, wind_speed, Cd, c_s, dot_M_star, M_star, R_star, M_env, R_env, L_star, &
        omega_star, a, e, M_planet, R_planet, k2dt_planet, verbose, dot_J_star)

      J_star = J_star + dot_J_star * dt
      omega_star = J_star / I_star

      if (verbose > 0) then
        dot_omega_star = (dot_J_star * I_star - J_star * dot_I_star) / I_star**2 ! Only for comparison with other formalisms.
        dot_J_orb = (M_star * M_planet * sqrt(G * a * (1. - e**2) / M_star)) - J_orb ! Only for comparison with other formalisms.
        J_orb = M_star * M_planet * sqrt(G * a * (1. - e**2) / M_star) ! Only for comparison with other formalisms.
        if (verbose > 1) then
          omega_planet = 2. * pi / sqrt(4. * pi**2 * a**3 / G * M_star) ! The planet is assumed synchronized with the star in this formalism, so its self-rotation period is the same as its orbital period.
          J_planet = I_planet * omega_planet
          J_syst = J_orb + J_star + J_planet
        end if
      end if

    elseif (formalism == 2) then ! Villaver's formalism.

      call Compute_V14(dot_a, dot_e, dot_omega_planet, G, M_star, a, dot_M_star, pi, wind_speed, Cd, R_planet, M_planet, c_s, &
        verbose, M_env, R_star, R_env, L_star, pim2, cF, e, Q_p, omega_planet)

      omega_planet = omega_planet + dot_omega_planet

    elseif (formalism == 3) then ! Home-made formalism.

      call Compute_Homemade (G, pi, wind_speed, Cd, c_s, k2dt_planet, M_planet, I_planet, M_star, a, dot_M_star, R_planet, e, &
        M_env, R_star, R_env, L_star, I_star, omega_star, J_orb, omega_planet, verbose, dot_a, dot_e, dot_J_star, &
        dot_J_orb, dot_omega_planet, dot_omega_star, dot_I_star)


      J_star = J_star + dot_J_star * dt
      J_orb  = J_orb  + dot_J_orb  * dt
      omega_planet = omega_planet + dot_omega_planet * dt
      omega_star   = omega_star   + dot_omega_star   * dt

      if (verbose >= 2) then
        J_planet = I_planet * omega_planet
        J_syst = J_orb + J_star + J_planet
      end if
    end if

    a = a + dot_a * dt
    e = e + dot_e * dt

    ! ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Ap = a * (1. + e) ! Apoapsis
    Pe = a * (1. - e) ! Periapsis

    if (e < 1.d-10) then ! If e is too small, set it to 0.
      e = 0.
    end if

    if (e < 0.) then
      write(Logfile,*) "Eccentricity set to 0 due to negative value."
      stop 'Eccentricity cannot become < 0 !'
      e = 0.
    end if

    ! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ! Writing results in the output file.

    if (R_eff_star > 0.) then ! Writing in the result file.
      if (verbose < 1) then ! If no details are asked.
        write(13, FMT='(ES16.10, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5)') &
          time / sec, M_star / M_Sun, R_eff_star / R_Sun, L_star / L_Sun, dot_M_star / M_Sun_year, M_env / M_Sun, &
          R_env / R_Sun, a / au, dot_a / au * sec, e, omega_planet, omega_star
      else if(verbose < 2) then ! If it is asked to print details.
        write(13, FMT='(ES16.10, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, &
          & ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5)') &
          time / sec, M_star / M_Sun, R_eff_star / R_Sun, L_star / L_Sun, dot_M_star / M_Sun_year, M_env / M_Sun, &
          R_env / R_Sun, a / au, dot_a / au * sec, e, omega_planet, omega_star, dot_omega_star, J_orb, dot_J_orb, J_star, &
          dot_J_star, I_star
      else
        write(13, FMT='(ES16.10, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, &
          & ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5, ES14.5)') &
          time / sec, M_star / M_Sun, R_eff_star / R_Sun, L_star / L_Sun, dot_M_star / M_Sun_year, M_env / M_Sun, &
          R_env / R_Sun, a / au, dot_a / au * sec, e, omega_planet, omega_star, dot_omega_star, J_orb, dot_J_orb, J_star, &
          dot_J_star, I_star, J_planet, J_syst, dot_omega_planet
      end if
    end if


    if (Pe <= R_eff_star) then ! If the planet reaches the star's outer layers, we stop the computations and record the time of impact.
      write(Logfile,*) "The planet collided with the star, end of the computations."
      write(*,*) "The planet collided with the star, end of the computations."
      write(13, *) "# Survival = 0   contact time : ", time / sec
      exit
    endif

    ! Determining the next timestep to keep delta_a/a and delta_e/e < eps and larger than 0.1*eps.

    call get_timestep(dt, dot_a, a, dot_e, e, dot_J_star, J_star, time, age_star_array(i-1), age_star_array(i), eps, error_check)
    if (error_check == 1) then
      write(Logfile,*) "***** Error with get_timestep function, end of the computations. *****"
      exit
    end if

    time = time + dt

    if (dt <= sec .and. force_continue == 0) then ! If the timestep is very small and the continuation is not forced.
      write(Logfile,*) "Stopping due to dt being too small."
      write(*,*) "Stopping due to dt being too small.", dt
      write(13, *) "# Survival = 1   contact time : ", -1

      exit ! We stop there.
    end if

  end do

  if (time > time_end) then
    write(Logfile,*) "End of time range reached, end of the computation."
    write(*,*) "End of time range reached, end of the computation."
    write(13, *) "# Survival = 1   contact time : ", -1
  end if

  close(13) ! Output file
  close(12) ! Closing STAREVOL model file.
  close(Logfile) ! Log file
  if (Special_print == 1) close(Special_print_file)

contains


  ! ##############################################################################################################################
  ! ##############################################################################################################################

   subroutine read_parameters(Stellar_model, a, e, R_planet, M_planet, Q_p, omega_planet, alpha_p, dt_max, &
        T_wind, time0, wind_speed, mu, eps, force_continue, verbose, formalism, k2dt_planet, omega_star, I_planet_coef)

    implicit none
    real, intent(out) :: a, e, R_planet, M_planet, Q_p, omega_planet, alpha_p, dt_max, T_wind, time0, wind_speed, mu, eps, &
      k2dt_planet, omega_star, I_planet_coef
    integer, intent(out) :: force_continue, verbose, formalism
    character(128), intent(out) :: Stellar_model
    logical :: iexist

    namelist /parameters/ Stellar_model, a, e, R_planet, M_planet, Q_p, omega_planet, alpha_p, dt_max, T_wind, time0, &
      wind_speed, mu, eps, force_continue, verbose, formalism, k2dt_planet, omega_star, I_planet_coef

    inquire (file='input.par', exist=iexist)
    if (.not. iexist) then
      print *,'You must provide the input : input.dat'
      stop 'Input file input.par not found'
    else
      open(unit=33, file='input.par')
      read(33, parameters)
    endif

    ! Convert inout parameters in SI units
    a          = a * au ! [m] (from au).
    R_planet   = R_planet * R_Jupiter ! [m] (from km).
    M_planet   = M_planet * M_Jupiter ! [kg]
    dt_max     = dt_max * sec ! [s] (from years).
    time0      = time0 * sec ! [s] (from years).
    wind_speed = wind_speed * unit_conv_dist

   end subroutine read_parameters

  ! ##############################################################################################################################
  ! ##############################################################################################################################

  subroutine get_timestep (dt, dot_a, a, dot_e, e, dot_J_star, J_star, time, age_star_ant, age_star_sup, eps, error_check)

    ! ----------------------------------------------------------------------------------------------------------------------------

    implicit none
    real, intent(inout) :: dt
    real, intent(in) :: a, dot_a, dot_e, e, dot_J_star, J_star
    real, intent(in) :: time, age_star_ant, age_star_sup, eps
    integer, intent(out) :: error_check

    real :: dt_a, dt_e, dt_J

    !write(Logfile,*) "Entering get_timestep subroutine."

    ! If the timestep is larger than the remaining time until next Starevol data point, we lower it to this point.
    if (time > age_star_sup) then
      write(Logfile,*) "####################################### time bug #######################################"
      write(Logfile,*) "Time is ahead of the next Starevol data point. Terminating program."
      write(Logfile,*) "########################################################################################"

      write(Logfile,*) "time : ", time / sec, " years"
      write(Logfile,*) "age_star_sup : ", age_star_sup / sec, " years,   age_star_ant : ", age_star_ant / sec, " years"
      error_check = 1
    end if

    dt_e = dt
    dt_a = dt
    dt_J = dt

    ! Constraint based on variation of a
    call calc_dt(dt_a, dot_a, a, time, age_star_sup, age_star_ant, eps)

    ! Constraint based on variation of e
    if (e > 1.e-3) then ! If there is an eccentricity.
      call calc_dt(dt_e, dot_e, e, time, age_star_sup, age_star_ant, eps)
      dt = min(dt_a, dt_e)
    else
      dt = dt_a
    endif

    if (formalism == 1 .or. formalism == 3) then ! If F23's formalism is used we check the evolution of omega_star.
      call calc_dt(dt_J, dot_J_star, J_star, time, age_star_sup, age_star_ant, eps)
      dt = min(dt, dt_J)
    end if

    if (force_continue == 1 .and. dt < sec) then ! If dt is too small and we want to force it to a bigger value.
      ! Checking if adding 2e-10 is more than the time to the next Starevol datapoint.
      if (verbose == 1) then
        write(Logfile,*) "Forcing dt. Value : ", dt
        write(Logfile,*) "   Time : ", time
        write(Logfile,*) "   age_star_sup : ", age_star_sup
        write(Logfile,*) "   delta t : ", age_star_sup - time
      end if

      if (time + sec > age_star_sup) then ! If adding a year to "time" makes id greater than the next Starevol point.
        dt = (age_star_sup - time) ! We fix the time to Starevol's point
      else
        dt = sec ! Otherwise we add set dt to 1 year.
      end if
    end if

    !write(Logfile,*) "New time step : ", dt, " s   (equal to ", dt / sec, " years)"
    !write(Logfile,*) "Exiting get_timestep subroutine."

  end subroutine get_timestep

  ! ##############################################################################################################################
  ! ##############################################################################################################################


  subroutine calc_dt(dt, dot_x, x, time, age_star_sup, age_star_ant, eps)
    ! Subroutine that computes the length of a time step.
    ! ----------------------------------------------------------------------------------------------------------------------------

    real, intent(in) :: dot_x, x, time, age_star_sup, age_star_ant, eps
    real, intent(inout) :: dt

    ! If the variation of "x" is greater than 1% of its value, we reduce the time step (to increase precision).
    do while (abs((dt * dot_x) / x) > eps)
      dt = dt * 0.8

      if (dt + time > age_star_sup) then
        write(Logfile,*) "time step reduced to avoid reaching age_star_sup."
        dt = age_star_sup - time
        exit
      end if

      write (Logfile,*) "timestep lowered to : ", dt / sec, " years    (abs(delta_x / x) = ", &
        abs((dt * dot_x) / x), ")   delta_x = ", dot_x * dt

      if (isnan(x) .eqv. .True.) then
        write(*,*) "x = ", x, "dot_x = ", dot_x, "dt = ", dt
        write(*,*) "*************************** ERROR *****************************"
        stop "*************** 'x' is a NaN ***************"
      end if
    end do

    ! If the variation of "x" is lower than 0.1% of its value, we increase the time step (to save computation time).
    do while (abs((dt * dot_x) / x) < (0.1 * eps))
      dt = dt * 1.2
      write (Logfile,*) "timestep increased to : ", dt / sec, " years    (abs(delta_x / x) = ", &
        abs((dt * dot_x) / x), ")   x = ", x

      ! While increasing the time step, if it reaches 50% of the difference between age_star_sup and age_star_ant,
      ! we cap it at this value (50% of (age_star_sup - age_star_ant))
      if (dt > 0.5 * (age_star_sup - age_star_ant)) then
        dt = 0.5 * (age_star_sup - age_star_ant)
        write(Logfile,*) "time step capped to ", dt / sec, " years."
        exit
      end if
    end do

  end subroutine calc_dt

  ! ##############################################################################################################################
  ! ##############################################################################################################################

  subroutine Compute_V14 (dot_a, dot_e, dot_omega_planet, G, M_star, a, dot_M_star, pi, wind_speed, Cd, R_planet, M_planet, c_s, &
      verbose, M_env, R_star, R_env, L_star, pim2, cF, e, Q_p, omega_planet) ! Villaver 2014's formalism
    real, intent(in) :: G, M_star, a, dot_M_star, pi, wind_speed, Cd, R_planet, M_planet, c_s
    real, intent(in) :: M_env, R_star, R_env, L_star, pim2, cF, e, Q_p, omega_planet
    integer, intent(in) :: verbose
    real, intent(inout) :: dot_a, dot_e, dot_omega_planet
    real :: dot_a_eq_tides_star

    v_orb = sqrt(G * (M_star + M_planet) / a) ! Approximation from sqrt(G * M_star * ((2 / r) - (1 / a)))

    rho = dot_M_star / (pim4 * a**2 * wind_speed) ! Formula 9 from Villaver 2009
    f_drag = 0.5 * Cd * rho * v_orb**2 * pi * R_planet**2 ! Formula 4 from Villaver 2014
    f_grav = (pim4 * (M_planet * G)**2 * rho * 0.5) / c_s**2 ! Formula 5 from Villaver 2014 | with mach number := 0.5

    tau = ((M_env * (R_star - R_env)**2) / (3. * L_star))**(1./3.) ! Formula 8 from Villaver 2014

    P_orb = sqrt((4. * pi**2 * a**3) / (G * M_star)) ! Kepler 3
    n = pim2 / P_orb

    f1 = (9./2.) * amin1(1., (pim2 / (1. * n * cF * tau))**2)
    f2 = (9./2.) * amin1(1., (pim2 / (2. * n * cF * tau))**2)
    f3 = (9./2.) * amin1(1., (pim2 / (3. * n * cF * tau))**2)

    g2 = 1. + (15./2.) * e**2 + (45./8.) * e**4 + (5./16.) * e**6 ! Formula 13 from Villaver 2014.
    g3 = 1. + (15./4.) * e**2 + (15./8.) * e**4 + (5./64.) * e**6 ! Formula 14 from Villaver 2014.
    g4 = 1. + (3./2.)  * e**2 + (1./8.)  * e**4 ! Formula 15 from Villaver 2014.
    g5 = 1. +  3.      * e**2 + (3./8.)  * e**4 ! Formula 16 from Villaver 2014.

    q = M_planet / M_star

    dot_a_eq_tides_star = (1. / (9. * tau)) * (M_env / M_star) * q * (1. + q) * &
      (R_Star / a)**8 * (2. * f2 + e**2 * ((7./8.) * f1 - 10. * f2 + (441./8.) * f3)) ! Formula 6 from Villaver 2014.

    ! --------------------------------------------------------------------------------------------------------------------------
    ! Computation of the eccentricity variation.

    dot_e_1 = e * ((-1. / (36. * tau)) * (M_env / M_star) * q * (1. + q) * (R_star / a)**8 * &
      ((5./4.) * f1 - (2. * f2) + (147. / 4.) * f3)) ! Formula 7 from Villaver 2014.
    dot_e_2 = e * ((-81. * n)/(2. * Q_p) / q * (R_planet / a)**5 * &
      (g3 / sqrt((1 - e**2)**13) - (11. * g4 * (1. - e**2)**(-5) * omega_planet) / (18. * n))) ! Formula 10 from Villaver 2014.

    dot_e = dot_e_1 + dot_e_2

    ! --------------------------------------------------------------------------------------------------------------------------
    ! Computation of the change in separation.

    dot_a_1 = a * (-dot_M_star / (M_star + M_planet) - (2. / (M_planet * v_orb) * (f_drag + f_grav)) - dot_a_eq_tides_star) ! Formula 3 from Villaver 2014.
    dot_a_2 = a * (((2. * e * dot_e) / (1. - e**2)) - ((9. / (Q_p * q)) * (R_planet / a)**5) * &
      (g2 * sqrt(1. - e**2)**(-13) * n - g5 * (1. - e**2)**(-5) * omega_planet)) ! Formula 11 from Villaver 2014.

    dot_a = dot_a_1 + dot_a_2

    ! --------------------------------------------------------------------------------------------------------------------------
    ! Computation of the variation of the rotation of the planet.

    dot_omega_planet = ((9. * n**2) / (2. * alpha_p * Q_p * q))   *   (R_planet / a)**3   *   &
      (g2 * (1. - e**2)**(-6) - g5 * sqrt(1. - e**2)**(-9) * (omega_planet / n)) ! Formula 12 from Villaver 2014

  end subroutine Compute_V14

  ! ##############################################################################################################################
  ! ##############################################################################################################################

  subroutine Compute_F23 (dot_a, dot_e, G, pi, wind_speed, Cd, c_s, dot_M_star, M_star, R_star, M_env, R_env, L_star, &
      omega_star, a, e, M_planet, R_planet, k2dt_planet, verbose, dot_J_star)

    real, intent(in) :: G, pi, wind_speed, Cd, c_s
    real, intent(in) :: dot_M_star, M_star, R_star, M_env, R_env, L_star, omega_star
    real, intent(in) :: a, e, M_planet, R_planet, k2dt_planet ! k2dt : k2 * delta_t planet
    integer, intent(in) :: verbose
    real, intent(inout) :: dot_a, dot_e, dot_J_star
    real :: dot_a_drag, dot_e_eq_tides_star, dot_a_eq_tides_star, dot_a_eq_tides_planet, dot_e_eq_tides_planet, dot_a_ang_mom
    real :: dot_J_eq_star, dot_J_Mloss_star, dot_J_sys
    real :: v_orb, rho, f_drag, f_grav, f_orb, P_orb, tau, q, omega_e, fe, n, f1, f2, rotation_planet
    real :: dot_J_star_tides, mu

    dot_a_ang_mom = -a * dot_M_star / (M_star + M_planet) ! Adaptation from formula 3 of Villaver 2014 (first term).

    v_orb = sqrt(G * M_star / a) ! Approximation from (G * M_star * ((2 / r) - (1 / a)))**(0.5)
    rho = dot_M_star / (4 * pi * a**2 * wind_speed) ! Formula 9 from Villaver 2009.
    f_drag = 0.5 * Cd * rho * v_orb**2 * pi * R_planet**2 ! Formula 4 from Villaver 2014.
    f_grav = (pim4 * (G * M_planet)**2 * rho * 0.5) / c_s**2 ! Formula 5 from Villaver 2014 | with mach number := 0.5

    dot_a_drag = -a * (2. / (M_planet * v_orb)) * (f_drag + f_grav) ! Adaptation from formula 3 of Villaver 2014 (second term).

    P_orb = ((4. * pi**2 * a**3) / (G * M_star))**0.5 ! Kepler 3

    tau = (M_env * (R_star - R_env)**2 / (3. * L_star))**(1./3.) ! Formula 9 of Fellay 2023.

    if ( tau < P_orb / 2. ) then ! Formula 10 of Fellay 2023 and paragraph above.
      f_orb = 1.
    else
      f_orb = (P_orb / (2. * tau))**2
    end if

    q = M_planet / M_star
    f1 = (1. + (31./2.) * e**2 + (255./8.) * e**4 + (185./16.) * e**6 + (25./64.) * e**8) / sqrt((1. - e**2)**15) ! Formula 6 of Fellay 2023.
    f2 = (1. + (15./2.) * e**2 + (45./8.) * e**4 + (5./16.) * e**6) / (1. - e**2)**6 ! Formula 5 of Fellay 2023.
    n = pim2 / P_orb

    dot_a_eq_tides_star = a * (f_orb / tau) * (M_env / M_star) * q * (1. + q) * (R_star / a)**8 * (f2 * (omega_star / n) - f1) ! Formula 3 of Fellay 2023.

    rotation_planet = (1. + (15./2.) * e**2 + (45./8.) * e**4 + (5./16.) * e**6) / &
      ((1. + 3. * e**2 + (3./8.)*e**4) * sqrt((1. - e**2)**3)) ! Formula 13 of Fellay 2023.

    dot_a_eq_tides_planet = a * 6. * k2dt_planet * (1. / q) * (R_planet / a)**5 * n**2 * (f2 * rotation_planet - f1) ! Formula 11 of Fellay 2023.

    omega_e = (1. + (3./2.) * e**2 + (1./8.) * e**4) / (1. - e**2)**5 ! Formula 7 of Fellay 2023.
    fe = (1. + (15./4.) * e**2 + (15./8.) * e**4 + (5./64.) * e**6) / sqrt((1. - e**2)**13) ! Formula 8 of Fellay 2023.

    dot_e_eq_tides_star = e * (11./4.) * (f_orb / tau) * (M_env / M_star) * q * (1. + q) * (R_star / a)**8 * &
      (omega_e * (omega_star / n) - (18./11.) * fe) ! Formula 4 of Fellay 2023.
    dot_e_eq_tidal_planet = e * (33./2.) * k2dt_planet * (1. / q) * (R_planet / a)**5 * n**2 * &
      (omega_e * rotation_planet - (18./11.) * fe) ! Formula 12 of Fellay 2023.

    dot_a = dot_a_ang_mom + dot_a_drag + dot_a_eq_tides_star + dot_a_eq_tides_planet ! Sum of the components computed above.
    dot_e = dot_e_eq_tides_star + dot_e_eq_tidal_planet ! Sum of the components computed above.


    dot_J_Mloss_star = -(2./3.) * dot_M_star * omega_star * R_star**2 ! Formula 23 of Siess 2013.
    !write (Logfile,*) "dot_J_Mloss_star : ", dot_J_Mloss_star
    !write (Logfile,*) "dot_M_star : ", dot_M_star
    !write (Logfile,*) "omega_star : ", omega_star
    !write (Logfile,*) "R_star : ", R_star

    mu = M_star * M_planet / (M_star + M_planet)

    dot_J_star_tides = - sqrt(G * mu * a * (1 - e**2)) * &
      (0.5 * (dot_a_eq_tides_star + dot_a_eq_tides_planet) / a - &
      (e / (1 - e**2)) * (dot_e_eq_tides_star + dot_e_eq_tidal_planet)) ! Formula 18 of Fellay 2023
    !dot_J_star_tides = 0.
    !write (Logfile,*) "dot_J_star_tides : ", dot_J_star_tides
    !write (Logfile,*) "  - sqrt(G * mu * a * (1 - e**2)) : ", - sqrt(G * mu * a * (1 - e**2))
    !write (Logfile,*) "  0.5 * (dot_a_eq_tides_star + dot_a_eq_tides_planet) / a : ", 0.5 * &
    !  (dot_a_eq_tides_star + dot_a_eq_tides_planet) / a


    dot_J_star = dot_J_Mloss_star + dot_J_star_tides
    !write (Logfile,*) "dot_J_star : ", dot_J_star
    !write (Logfile,*) "  dot_J_star_tides : ", dot_J_star_tides
    !write (Logfile,*) "  dot_J_Mloss_star : ", dot_J_Mloss_star

    if (Special_print == 1) then
      write(Special_print_file, FMT='(I6, ES15.6, ES15.6, ES15.6, ES15.6, ES15.6, ES15.6)') &
        iter, time/sec, dot_a/au, dot_a_drag/au, dot_a_eq_tides_star/au, dot_a_eq_tides_planet/au, &
          (dot_a_ang_mom + dot_a_eq_tides_star + dot_a_eq_tides_planet)/au
    endif

  end subroutine Compute_F23

  ! ##############################################################################################################################
  ! ##############################################################################################################################

  subroutine Compute_Homemade (G, pi, wind_speed, Cd, c_s, k2dt_planet, M_planet, I_planet, M_star, a, dot_M_star, R_planet, e, &
    M_env, R_star, R_env, L_star, I_star, omega_star, J_orb, omega_planet, verbose, dot_a, dot_e, dot_J_star, dot_J_orb, &
    dot_omega_planet, dot_omega_star, dot_I_star)

    real, intent(in) :: G, pi, wind_speed, Cd, c_s
    real, intent(in) :: M_star, dot_M_star, M_env, R_star, R_env, L_star, omega_star, I_star, J_orb, dot_I_star
    real, intent(in) :: a, e, R_planet, M_planet, omega_planet, I_planet, k2dt_planet
    real, intent(inout) :: dot_a, dot_e, dot_J_star, dot_J_orb, dot_omega_planet, dot_omega_star
    integer, intent(in) :: verbose
    real :: dot_a_drag, dot_e_eq_tides_star, dot_a_eq_tides_star, dot_a_eq_tides_planet, dot_e_eq_tides_planet
    real :: dot_J_eq_star, dot_J_Mloss_star, dot_J_sys, dot_a_tides_ang_mom
    real :: v_orb, rho, f_drag, f_grav, f_orb, P_orb, tau, q, omega_e, fe, n, f1, f2, rotation_planet
    real :: dot_J_planet      ! Variation of the angular momentum of the planet [kg m2 s-2].

    v_orb = sqrt(G * M_star / a) ! Approximation from (G * M_star * ((2 / r) - (1 / a)))**(0.5)
    rho = dot_M_star / (4. * pi * a**2 * wind_speed) ! Formula 9 from Villaver 2009.
    f_drag = 0.5 * Cd * rho * v_orb**2 * pi * R_planet**2 ! Formula 4 from Villaver 2014.
    f_grav = (4. * pi * (G * M_planet)**2 * rho * 0.5) / c_s**2 ! Formula 5 from Villaver 2014 | with mach number := 0.5
    dot_a_drag = -a * (2. / (M_planet * v_orb)) * (f_drag + f_grav) ! Adaptation from formula 3 of Villaver 2014 (second term).

    P_orb = ((4. * pi**2 * a**3) / (G * M_star))**0.5 ! Kepler 3

    tau = (M_env * (R_star - R_env)**2 / (3. * L_star))**(1./3.) ! Formula 9 of Fellay 2023.

    if (tau < P_orb / 2.) then ! Formula 10 of Fellay 2023 and paragraph above.
      f_orb = 1.
    else
      f_orb = (P_orb / (2. * tau))**2
    end if

    omega_e = (1. + (3./2.) * e**2 + (1./8.) * e**4) / (1. - e**2)**5 ! Formula 7 of Fellay 2023.
    fe = (1. + (15./4.) * e**2 + (15./8.) * e**4 + (5./64.) * e**6) / sqrt((1. - e**2)**13) ! Formula 8 of Fellay 2023.

    n = 2. * pi / P_orb
    q = M_planet / M_star
    f1 = (1. + (31./2.) * e**2 + (255./8.) * e**4 + (185./16.) * e**6 + (25./64.) * e**8) / sqrt((1. - e**2)**15) ! Formula 6 of Fellay 2023.
    f2 = (1. + (15./2.) * e**2 + (45./8.) * e**4 + (5./16.) * e**6) / (1. - e**2)**6 ! Formula 5 of Fellay 2023.

    dot_e_eq_tides_star = e * (11./4.) * (f_orb / tau) * (M_env / M_star) * q * (1. + q) * (R_star / a)**8 * &
      (omega_e * (omega_star / n) - (18./11.) * fe) ! Formula 4 of Fellay 2023.

    dot_a_eq_tides_star = a * (f_orb / tau) * (M_env / M_star) * q * (1. + q) * (R_star / a)**8 * (f2 * (omega_star / n) - f1) ! Formula 3 of Fellay 2023.

    dot_J_eq_star = J_orb * (((1./2.) * dot_a_eq_tides_star / a) - ((e / (1 - e**2)) * dot_e_eq_tides_star))

    dot_J_Mloss_star = -(2./3.) * dot_M_star * omega_star * R_star**2

    dot_J_star = dot_J_eq_star

    rotation_planet = omega_planet / n

    dot_a_eq_tides_planet = a * 6. * k2dt_planet * (1. / q) * (R_planet / a)**5 * n**2 * (f2 * rotation_planet - f1) ! Formula 11 of Fellay 2023.
    dot_e_eq_tides_planet = e * (33./2.) * k2dt_planet * (1. / q) * (R_planet / a)**5 * n**2 * &
      (omega_e * rotation_planet - (18./11.) * fe) ! Formula 12 of Fellay 2023.

    dot_J_planet = J_orb * (((1./2.) * dot_a_eq_tides_planet / a) - ((e / (1 - e**2)) * dot_e_eq_tides_planet))

    dot_J_sys = dot_M_star * (M_planet / M_star) * J_orb / (M_star + M_planet)

    dot_e = dot_e_eq_tides_star + dot_e_eq_tides_planet ! Sum of the components computed above.

    dot_a_tides_ang_mom = a * (2. * ((dot_J_star + dot_J_planet + dot_J_sys) / J_orb) - (2. * dot_M_star / M_star) + &
      dot_M_star / (M_star + M_planet) + (2. * e * dot_e) / (1. - e**2)) ! Stellar and planetary tides components as well as mass loss component included.

    dot_a = dot_a_tides_ang_mom + dot_a_drag

    dot_J_orb = J_orb * (((1./2.) * dot_a / a) + (dot_M_star / M_star) - (1./2. * dot_M_star / (M_star + M_planet)) - &
      ((e / (1 - e**2)) * dot_e_eq_tides_star))

    dot_omega_planet = dot_J_planet / I_planet
    dot_omega_star   = (dot_J_star * I_star - J_star * dot_I_star) / I_star**2

    if (Special_print == 1) then
      write(Special_print_file, FMT='(I6, ES15.6, ES15.6, ES15.6, ES15.6, ES15.6, ES15.6)') &
        iter, time/sec, dot_a/au, dot_a_drag/au, dot_a_eq_tides_star/au, dot_a_eq_tides_planet/au, dot_a_tides_ang_mom/au
    endif

  end subroutine Compute_Homemade

end program Sekhmet
