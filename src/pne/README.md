# PineMP - Pin Aware Macro Placer

PineMP is a new macroplacer design to explore the relationship between Floorplanning steps, using Pin Asssignment techniques to better optimize the placement.

## Commands

### pine_mp
Performs macro placement using the PineMP algorithm.

```tcl
pine_mp [-num_threads num_threads]
```

**Options:**
- `-num_threads`: Number of threads to use for placement (default: 1)

**Example:**
```tcl
pine_mp -num_threads 4
```

### pine_mp_debug
Enables or disables debug output for PineMP.

```tcl
pine_mp_debug [-debug]
```

**Options:**
- `-debug`: Enable debug output

**Example:**
```tcl
pine_mp_debug -debug
```

## Implementation Status

This is a blank project template following OpenROAD development guidelines. The core placement algorithm needs to be implemented in the `PineMP::place()` method.

## Development

The project follows the standard OpenROAD tool structure:
- `src/`: Source files and implementation
- `include/pne/`: Public header files
- `test/`: Test cases and regression tests

## Dependencies

- OpenROAD database (odb)
- OpenROAD utilities (utl)
- OpenSTA timing analysis
- TCL interface