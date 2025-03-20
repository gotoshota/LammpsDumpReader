# LammpsDumpReader Fortran Module Guide

## Build & Test Commands
```bash
# Compile the module
gfortran -c lammpsio.f90

# Compile and run the test program
gfortran -o test_lammpstrj lammpsio.f90 ../test/main.f90
./test_lammpstrj

# Run with specific input file
sed -i 's|input_filename = ".*"|input_filename = "../data/yourfile.lammpstrj"|' ../test/main.f90
```

## Code Style Guidelines
- **Module Structure**: Fortran modules with clear type definitions and subroutines
- **Naming**: Use descriptive lowercase names with underscores (snake_case)
- **Documentation**: Use Doxygen-style comments with `!>` prefix for documentation
- **Error Handling**: Use print statements and `stop` for fatal errors
- **Types**: Define structured types with clear component descriptions
- **Memory Management**: Explicitly allocate/deallocate arrays as needed
- **Parameters**: Use `intent(in)`, `intent(out)`, or `intent(inout)` for all parameters
- **Functions**: Include detailed parameter descriptions in comments
- **Internationalization**: Comments may be in Japanese or English