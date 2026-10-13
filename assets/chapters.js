// The table of contents. The sidebar, the home page and the previous/next links are built from it.
window.CHAPTERS = [
  { part: 'Getting started', id: '00-introduction', title: 'Introduction' },
  { part: 'Getting started', id: '01-sandbox', title: 'The Sandbox' },

  { part: 'Foundations', id: '02-particles', title: 'Particles and Newton' },
  { part: 'Foundations', id: '03-time-integration', title: 'Time Integration' },
  { part: 'Foundations', id: '04-springs-damping', title: 'Springs, Damping and Energy' },

  { part: 'Constraints', id: '05-constraints', title: 'Constraints and Multipliers' },
  { part: 'Constraints', id: '06-pbd', title: 'Position-Based Dynamics' },
  { part: 'Constraints', id: '07-xpbd', title: 'XPBD: Compliance' },
  { part: 'Constraints', id: '08-substepping', title: 'Small Steps' },
  { part: 'Constraints', id: '09-linear-solvers', title: 'Iterative and Direct Solvers' },

  { part: 'Rods', id: '10-rotations', title: 'Rotations and Rigid Bodies' },
  { part: 'Rods', id: '11-beams', title: 'Bending, Twisting, Buckling' },
  { part: 'Rods', id: '12-cosserat-rods', title: 'Cosserat Rods' },
  { part: 'Rods', id: '13-tendons', title: 'Cable-Driven Bending' },

  { part: 'Contact', id: '14-collision-detection', title: 'Collision Detection' },
  { part: 'Contact', id: '15-contact-friction', title: 'Contact and Friction' },

  { part: 'Soft bodies', id: '16-cloth-shells', title: 'Cloth and Thin Shells' },
  { part: 'Soft bodies', id: '17-pressure', title: 'Pressure and Inflated Bodies' },
  { part: 'Soft bodies', id: '18-coupling', title: 'Coupling Two Solvers' },
  { part: 'Soft bodies', id: '19-active-materials', title: 'Muscles and Travelling Waves' },

  { part: 'Continuum and MPM', id: '20-deformation', title: 'The Deformation Gradient' },
  { part: 'Continuum and MPM', id: '21-hyperelasticity', title: 'Stress and Hyperelasticity' },
  { part: 'Continuum and MPM', id: '22-pic-flip', title: 'Particles Meet Grids: PIC and FLIP' },
  { part: 'Continuum and MPM', id: '23-apic', title: 'APIC' },
  { part: 'Continuum and MPM', id: '24-mpm', title: 'The Material Point Method' },
  { part: 'Continuum and MPM', id: '25-plasticity', title: 'Plasticity: Snow and Sand' },
  { part: 'Continuum and MPM', id: '26-implicit-mpm', title: 'Implicit MPM' },
  { part: 'Continuum and MPM', id: '27-pbmpm', title: 'Position-Based MPM' },

  { part: 'Reading the papers', id: '28-reading-apic', title: 'Reading: APIC (2015)' },
  { part: 'Reading the papers', id: '29-reading-mpm-course', title: 'Reading: The MPM Course (2016)' },
  { part: 'Reading the papers', id: '30-reading-pbmpm', title: 'Reading: PB-MPM (2024)' },

  { part: 'Making it real', id: '31-stability', title: 'Hunting Explosions' },
  { part: 'Making it real', id: '32-verification', title: 'Verifying a Simulator' },
];
