# Advanced Usage

## Customizing Simulation Settings

### Controlling Simulation Steps

Options to skip specific phases:
```bash
# Skip preprocessing
./exprorer_msmd protocol.yaml --skip-preprocess

# Skip simulation
./exprorer_msmd protocol.yaml --skip-simulation

# Skip postprocessing
./exprorer_msmd protocol.yaml --skip-postprocess
```

### Independent Trial Settings

```yaml
general:
  iter_index: 0,1,2  # Execute 3 independent trials
  # Notation:
  # "1-3" => 1,2,3
  # "5-9:2" => 5,7,9
  # "1-3,5-9:2" => 1,2,3,5,7,9
```

### Adjusting Simulation Parameters

```yaml
exprorer_msmd:
  general:
    dt: 0.002          # Time step (ps)
    temperature: 300   # Temperature (K)
    pressure: 1.0      # Pressure (bar)
    pbc: xyz          # Periodic boundary conditions

  sequence:
    - name: pr        # Production run
      type: production
      nsteps: 20000000  # Number of steps (40 ns)
      nstxtcout: 5000   # Output frequency (10 ps)
```

## Customizing Analysis Settings

### PMAP Calculation Settings

```yaml
map:
  type: pmap
  snapshot: 2001-4001:1  # Snapshot range to use
  maps:
    - suffix: nVH        # Heavy atoms only
      selector: (!@VIS)&(!@H*)
    - suffix: nV         # All atoms
      selector: (!@VIS)
  map_size: 80           # Map size (Å)
  normalization: total   # Normalization method, 'snapshot' and 'GFE' are also available
```

### Inverse MSMD Related Settings

```yaml
probe_profile:
  resenv:
    map: nVH            # Map to use
    threshold: 0.001    # Probability threshold
  profile:
    types:
      - name: anion # Map name
        atoms:
          - ["ASP", " CB "]
          - ["GLU", " CB "]