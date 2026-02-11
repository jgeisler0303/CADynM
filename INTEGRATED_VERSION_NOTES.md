# Integrated Version - Multi-Backend Support

## Overview
The `integrated_version` branch combines the functionality of both the `main` branch (MATLAB Symbolic Toolbox) and the `mamas` branch (MAMaS/Maxima interface) into a single unified codebase.

## Key Changes

### 1. **Removed Type Annotations from Symbolic Properties**
   - Properties that previously had type specifications like `msym` or `sym` now use generic types or no type
   - Examples:
     - `gravity (:,1) msym = []` → `gravity (:,1) = []`
     - `q (:,1) msym = []` → `q (:,1) = []`
   - Files affected: `MultiBodySystem.m`, `Body.m`, `ElasticBody.m`, `Parameters.m`

### 2. **Added `symbolicBackend` Property to MultiBodySystem**
   - New property: `symbolicBackend (1,1) string = "msym"`
   - Selectable values: `"sym"` for MATLAB Symbolic Toolbox, `"msym"` for MAMaS/Maxima (default)
   - Set via constructor 4th argument: `MultiBodySystem(name, dof, input, symbolicBackend)`

### 3. **Symbolic Abstraction Layer**
   New methods in `MultiBodySystem` provide centralized symbolic creation:

   #### Core Factory Methods
   - **`createSymbolic(varargin)`** - Main factory for creating symbolic variables
     - Handles numeric values: `createSymbolic(42)`
     - Handles variable names: `createSymbolic('x')` or `createSymbolic('x', 'real')`
     - Delegates to backend-specific helpers

   - **`symbolicEye(n, m)`** - Create identity matrix
   - **`symbolicZeros(rows, cols)`** - Create zero matrix

   #### Backend-Specific Helpers (private)
   - **`createSymbolicSym(varargin)`** - MATLAB Symbolic Toolbox implementation
   - **`createSymbolicMsym(varargin)`** - MAMaS/Maxima implementation

### 4. **Updated Dependencies to Use Factory Methods**
   All direct `msym(...)` and `sym(...)` calls replaced with `obj.createSymbolic(...)`:
   - In `MultiBodySystem.m`: Constructor, derivatives, internal symbolic creations
   - In `Body.m`: Coordinate transformations, rotations, external transformations
   - In `ElasticBody.m`: Elastic element initializations
   - In `ElasticTaylor.m`: Taylor expansion element handling
   - In `Parameters.m`: Parameter symbolic variable creation

### 5. **System Reference in Parameters**
   - `Parameters` class now accepts a system reference: `Parameters(struct, system)`
   - Uses system's `createSymbolic` method if available, defaults to `msym`
   - Allows proper symbolic backend selection when parameters are created

### 6. **Method Adjustments**
   - `Body.rotationMatrix()` changed from static to instance method to access system's symbolic backend
   - Updated all callers: `Body.rotationMatrix(...)` → `obj.rotationMatrix(...)`

## Usage

### Creating a System with MATLAB Symbolic Toolbox
```matlab
system = MultiBodySystem('MySystem', 'q1', 'u1', 'sym');
```

### Creating a System with MAMaS/Maxima (default)
```matlab
system = MultiBodySystem('MySystem', 'q1', 'u1');
% or explicitly:
system = MultiBodySystem('MySystem', 'q1', 'u1', 'msym');
```

### Checking Current Backend
```matlab
backend = system.symbolicBackend;
fprintf('Using backend: %s\n', backend);
```

## Maintenance Strategy

### Minimal Conditional Code
The design minimizes parts that require conditional execution:

1. **All symbolic creation goes through factory methods** - No backend-specific code scattered throughout
2. **Type annotations removed** - No need for conditional type casting
3. **Argument specifications simplified** - Methods accept generic input types
4. **Centralized backend selection** - Only `createSymbolic*` methods contain backend-specific logic

### Adding New Features
When adding features that create symbolic variables:
1. Use `obj.createSymbolic(...)` instead of `msym(...)` or `sym(...)`
2. Use `obj.symbolicEye(...)` instead of `msym(eye(...))` or `sym(eye(...))`
3. Use `obj.symbolicZeros(...)` instead of hardcoding zeros and converting

### Conditional Patterns (if needed)
For operations that differ significantly between backends, add helper methods:
```matlab
% In MultiBodySystem.m
function result = backendSpecificOperation(obj, ...)
    if strcmp(obj.symbolicBackend, "sym")
        % MATLAB Symbolic version
        result = ...;
    else
        % MAMaS/Maxima version
        result = ...;
    end
end
```

## Files Modified

- `toolbox/MultiBodySystem.m` - Core integration, factory methods
- `toolbox/Body.m` - Removed type specs, use system's factory
- `toolbox/ElasticBody.m` - Removed type specs, use system's factory
- `toolbox/ElasticTaylor.m` - Updated to use system's factory
- `toolbox/Parameters.m` - Accept system reference, use backend-aware creation

## Testing

To verify both backends work correctly:
1. Create system with `'sym'` backend and run existing tests
2. Create system with `'msym'` backend and run existing tests
3. Compare results to ensure both produce equivalent models

## Future Enhancements

- Consider lazy initialization of symbolic backend (allow switching after construction with warning)
- Add utility method to convert models between backends
- Document performance characteristics of each backend for different model sizes
