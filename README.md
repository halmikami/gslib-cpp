⏺ # gslib-cpp                                                                                                                             
                                                                                                                                  
  Modern C++17 implementation of [GSLIB](http://www.statios.com/GSLIB/)                                                                   
  geostatistical algorithms with Python bindings.                                                                                         
                                                                                                                                          
  ## What is this?                                                                                                                        
                                                                                                                                          
  A drop-in replacement for GSLIB's Fortran programs, written in modern C++                                                               
  with OpenMP parallelization and Python bindings via pybind11.
                                                                                                                                          
  **No Fortran compiler needed.** No f2py. No gfortran version issues.                                                                    
   
  ## Implemented Algorithms                                                                                                               
                                                        
  | Function | GSLIB Equivalent | Description |                                                                                           
  |----------|-----------------|-------------|
  | `gamv()` | GAMV | Experimental variogram for irregular data |                                                                         
  | `kt3d()` | KT3D | Ordinary / Simple Kriging (3D) |                                                                                    
  | `cova3()` | COVA3 | Covariance model evaluation |                                                                                     
  | `declus()` | DECLUS | Cell declustering |                                                                                             
                                                                                                                                          
  ## Installation                                       
                                                                                                                                          
  ```bash                                               
  pip install gslib-cpp

  Or build from source:                                                                                                                   
   
  git clone https://github.com/halmikami/gslib-cpp.git                                                                                    
  cd gslib-cpp                                          
  pip install .
                                                                                                                                          
  Requirements
                                                                                                                                          
  - C++17 compiler (GCC 7+, Clang 5+, MSVC 2019+)       
  - CMake 3.14+
  - Python 3.9+                                                                                                                           
  - OpenMP (optional, for parallelization)
                                                                                                                                          
  Quick Start                                           

  import gslib_cpp
  import numpy as np

  # Load drillhole data                                                                                                                   
  x = np.array([...])  # easting
  y = np.array([...])  # northing                                                                                                         
  z = np.array([...])  # elevation                      
  values = np.array([...])  # grade (e.g., Cu %)                                                                                          
                                                                                                                                          
  # Experimental variogram                                                                                                                
  result = gslib_cpp.gamv(                                                                                                                
      x, y, z, values,                                  
      n_lags=15,
      xlag=20.0,       # lag distance (m)                                                                                                 
      xltol=10.0,      # lag tolerance                                                                                                    
      azm=0.0,         # azimuth (0=N)                                                                                                    
      atol=90.0,        # angular tolerance                                                                                               
      bandwh=1e10,      # horizontal bandwidth                                                                                            
  )                                                                                                                                       
  print(result.lags, result.semivariance, result.pair_counts)                                                                             
                                                                                                                                          
  # Ordinary Kriging                                                                                                                      
  estimates = gslib_cpp.kt3d(                                                                                                             
      x, y, z, values,                                                                                                                    
      bhid=drillhole_ids,                               
      xout=grid_x, yout=grid_y, zout=grid_z,                                                                                              
      search_radius=200.0,
      ndmax=16, ndmin=4, noct=4,                                                                                                          
      nugget=0.1,                                                                                                                         
      model_types=[1],     # 1=spherical                                                                                                  
      model_cc=[0.9],      # sill contribution                                                                                            
      model_aa=[150.0],    # range                                                                                                        
  )                                                                                                                                       
  print(estimates.values, estimates.variances)                                                                                            
                                                                                                                                          
  Why not just use PyGSLIB?
                                                                                                                                          
  ┌─────────────────┬────────────────────────────────────┬───────────────────────────────┐                                                
  │                 │              PyGSLIB               │           gslib-cpp           │
  ├─────────────────┼────────────────────────────────────┼───────────────────────────────┤                                                
  │ Language        │ Fortran 77 via f2py                │ C++17 via pybind11            │
  ├─────────────────┼────────────────────────────────────┼───────────────────────────────┤
  │ Build           │ Requires gfortran, Python 3.8 only │ Any C++ compiler, Python 3.9+ │                                                
  ├─────────────────┼────────────────────────────────────┼───────────────────────────────┤                                                
  │ Parallelization │ Single-threaded                    │ OpenMP multi-threaded         │                                                
  ├─────────────────┼────────────────────────────────────┼───────────────────────────────┤                                                
  │ Spatial index   │ Brute force O(n²)                  │ KDTree O(n log n)             │
  ├─────────────────┼────────────────────────────────────┼───────────────────────────────┤                                                
  │ Maintenance     │ Last updated 2021                  │ Active                        │
  ├─────────────────┼────────────────────────────────────┼───────────────────────────────┤                                                
  │ API             │ Dict-based (GSLIB style)           │ Pythonic (keyword args)       │
  └─────────────────┴────────────────────────────────────┴───────────────────────────────┘                                                
   
  Compatibility                                                                                                                           
                                                        
  Results are numerically equivalent to GSLIB Fortran within floating-point                                                               
  precision (tested with np.testing.assert_allclose(rtol=1e-6)).
                                                                                                                                          
  References                                            
                                                                                                                                          
  - Deutsch, C.V. & Journel, A.G. (1997). GSLIB: Geostatistical Software                                                                  
  Library and User's Guide. Oxford University Press.
  - Journel, A.G. & Huijbregts, C.J. (1978). Mining Geostatistics.                                                                        
  Academic Press.                                       
                                                                                                                                          
  License                                               
                                                                                                                                          
  MIT. Original GSLIB algorithms are public domain (Stanford University).                                                                 
   
  Contributing                                                                                                                            
                                                        
  PRs welcome. Please include tests that compare output against GSLIB Fortran.    
