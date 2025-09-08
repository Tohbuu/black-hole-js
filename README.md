# 3D Black Hole Simulation

A combined JavaScript (Three.js) visualization and C++ headless simulation exploring black hole effects: gravity, photon capture, simple time-dilation, and CSV exports for offline analysis.

## Quick links (open these in the repo)
- Web demo / main UI: [blackhole.html](blackhole.html) — see [`BlackHole`](blackhole.html), [`MatterManager`](blackhole.html), [`Simulation`](blackhole.html)
- Test/demo scene: [tests/blackhole1.html](tests/blackhole1.html) — contains [`Photon`](tests/blackhole1.html) and [`LightSource`](tests/blackhole1.html)
- C++ simulation & utilities: [copy.cpp](copy.cpp) — implements [`BlackHole`](copy.cpp), [`Particle`](copy.cpp), [`BlackHoleSimulation::runSimulation`](copy.cpp), [`visualizeSpacetimeCurvature`](copy.cpp)
- Alternative C++ helpers: [black-hole.cpp](black-hole.cpp) — contains [`Vector3`](black-hole.cpp), [`getGravitationalAcceleration`](black-hole.cpp), [`getTimeDilationFactor`](black-hole.cpp)
- Sample outputs: [blackhole_simulation.csv](blackhole_simulation.csv), [spacetime_curvature.csv](spacetime_curvature.csv)

## Features
- Interactive 3D scene with orbit controls, accretion disk, corona and particle tracers ([blackhole.html](blackhole.html)).
- Matter creation (click or UI), trails, capture effects and simple Newtonian gravity approximation.
- Headless C++ simulation for higher-precision offline runs and CSV export ([copy.cpp](copy.cpp)).

## Run / Development

Web UI (recommended)
- Serve the repository folder (avoids CORS issues):
  - Python 3: python -m http.server 8000
  - Open: http://localhost:8000/blackhole.html
- Or open [blackhole.html](blackhole.html) directly (some features may require a server).

C++ simulation (command-line)
- Build:
  - g++ -std=c++11 copy.cpp -o bh_sim -lm
- Run:
  - ./bh_sim
- Outputs:
  - [blackhole_simulation.csv](blackhole_simulation.csv) — particle traces from [`BlackHoleSimulation::runSimulation`](copy.cpp)
  - [spacetime_curvature.csv](spacetime_curvature.csv) — curvature / time-dilation grid from [`visualizeSpacetimeCurvature`](copy.cpp)

## Notes & pointers
- JS uses simplified Newtonian-style gravity for visuals (`getGravitationalAcceleration` in [blackhole.html](blackhole.html) / [black-hole.cpp](black-hole.cpp)).
- Units: many values scaled by `METERS_PER_UNIT` in [blackhole.html](blackhole.html).
- The C++ code ([copy.cpp](copy.cpp)) can be extended to output JSON/binary or to increase numerical accuracy.

## Contributing
- Replace simplified accelerations in [blackhole.html](blackhole.html) with geodesic integration for more realism.
- Add automated tests / visual scenarios in [tests/blackhole1.html](tests/blackhole1.html).
- Extend C++ output options and formats in [copy.cpp](copy.cpp).

License: Add your preferred license file if

