enable_language(Fortran)
add_executable(1DDG Driver.f90 TimeIntegration.f90 DGMeshClass.f90 DGElementClass.f90 NodalDGStorage.f90 DGtoolbox.f90 FluxRoutines.f90 Diverses.f90)
