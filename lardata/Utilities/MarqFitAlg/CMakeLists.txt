#set(ENV{SPACK_CC} /packages/intel/20/compilers_and_libraries_2020/linux/bin/intel64/icc)
#set(ENV{SPACK_COMPILER_IMPLICIT_RPATHS} /packages/intel/20/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin)
#set(ENV{SPACK_COMPILER_SPEC} intel@19.1.1.217)
#set(ENV{SPACK_CPPFLAGS} '-Wno-deprecated -Wno-deprecated-declarations -xCORE-AVX512')
#set(ENV{SPACK_CXX} /packages/intel/20/compilers_and_libraries_2020/linux/bin/intel64/icpc)
#set(ENV{SPACK_ENV_PATH} /home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env:/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/case-insensitive:/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel)
#set(ENV{SPACK_F77} /packages/intel/20/compilers_and_libraries_2020/linux/bin/intel64/ifort)
#set(ENV{SPACK_FC} /packages/intel/20/compilers_and_libraries_2020/linux/bin/intel64/ifort)
#set(CC /home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel/icc)
#set(CXX /home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel/icc)
#set(FC /home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel/ifort)
#set(F77 /home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel/ifort)

add_custom_command(TARGET gshf_MarqFitAlg 
                   COMMAND ${CMAKE_COMMAND} -E env SPACK_CC=/packages/intel/20/compilers_and_libraries_2020/linux/bin/intel64/icc SPACK_COMPILER_IMPLICIT_RPATHS=/packages/intel/20/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin SPACK_COMPILER_SPEC=intel@19.1.1.217 SPACK_CPPFLAGS='-Wno-deprecated -Wno-deprecated-declarations -xCORE-AVX512' SPACK_CXX=/packages/intel/20/compilers_and_libraries_2020/linux/bin/intel64/icpc SPACK_ENV_PATH=/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env:/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/case-insensitive:/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel SPACK_F77=/packages/intel/20/compilers_and_libraries_2020/linux/bin/intel64/ifort SPACK_FC=/packages/intel/20/compilers_and_libraries_2020/linux/bin/intel64/ifort CC=/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel/icc CXX=/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel/icpc F77=/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel/ifort FC=/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/intel/ifort
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E env SPACK_CC=/home/users/sberkman/spack_dec17/gcc/8.2.0/linux-ubuntu18.04-x86_64-gcc-9.3.0-5nayizasz6ynib4fcem6ngcydi2ytf6w/bin/gcc SPACK_COMPILER_IMPLICIT_RPATHS=/home/users/sberkman/spack_dec17/gcc/8.2.0/linux-ubuntu18.04-x86_64-gcc-9.3.0-5nayizasz6ynib4fcem6ngcydi2ytf6w/lib/gcc/x86_64-pc-linux-gnu/8.2.0:/home/users/sberkman/spack_dec17/gcc/8.2.0/linux-ubuntu18.04-x86_64-gcc-9.3.0-5nayizasz6ynib4fcem6ngcydi2ytf6w/lib64 SPACK_COMPILER_SPEC=gcc@8.2.0 SPACK_CPPFLAGS='-Wno-deprecated -Wno-deprecated-declarations' SPACK_CXX=/home/users/sberkman/spack_dec17/gcc/8.2.0/linux-ubuntu18.04-x86_64-gcc-9.3.0-5nayizasz6ynib4fcem6ngcydi2ytf6w/bin/g++ SPACK_ENV_PATH=/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env:/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/case-insensitive:/home/users/sberkman/spack_dec17/spack/rollout2/NULL/lib/spack/env/gcc SPACK_F77=/home/users/sberkman/spack_dec17/gcc/8.2.0/linux-ubuntu18.04-x86_64-gcc-9.3.0-5nayizasz6ynib4fcem6ngcydi2ytf6w/bin/gfortran SPACK_FC=/home/users/sberkman/spack_dec17/gcc/8.2.0/linux-ubuntu18.04-x86_64-gcc-9.3.0-5nayizasz6ynib4fcem6ngcydi2ytf6w/bin/gfortran
)

include(FindOpenMP)
if(OPENMP_FOUND)
  # even if OpenMP is found on a SLF6 machine, it cannot be used.
  # see the discussion in https://cdcvs.fnal.gov/redmine/issues/23001
  if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6.32.*" )
    message(STATUS "OpenMP is not supported on SLF6")
  else()
    add_definitions(-DWITH_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif(OPENMP_FOUND)

add_library (gshf_MarqFitAlg STATIC MarqFitAlg.cxx)
set_property(TARGET gshf_MarqFitAlg PROPERTY POSITION_INDEPENDENT_CODE ON)

install_headers()
install_fhicl()
install_source()
