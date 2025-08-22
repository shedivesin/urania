// Convert a date (in milliseconds) to J2000 days or centuries.
function j2000(date) {
  return (date - 946728000000) / 86400000;
}

function j2000c(date) {
  return j2000(date) / 36525;
}


// Convert Keplerian orbital elements to Cartesian coordinates.
// https://ssd.jpl.nasa.gov/planets/approx_pos.html
function kepler(a, e, I, L, w1, N) {
  // Compute the argument of perihelion and the mean anomaly.
  const w = w1 - N;
  const M = L - w1;

  // Solve Kepler's equation.
  let E = M + e * Math.sin(M);
  for(;;) {
    const dM = M - (E - e * Math.sin(E));
    const dE = dM / (1 - e * Math.cos(E));
    E += dE;

    if(Math.abs(dE) <= 2e-8) {
      break;
    }
  }

  // Compute coordinates in orbital plane.
  const u = a * (Math.cos(E) - e);
  const v = a * Math.sqrt(1 - e * e) * Math.sin(E);

  // Convert to (and return) coordinates in ecliptic plane.
  const sin_I = Math.sin(I);
  const cos_I = Math.cos(I);
  const sin_N = Math.sin(N);
  const cos_N = Math.cos(N);
  const sin_w = Math.sin(w);
  const cos_w = Math.cos(w);
  return [
    u * (cos_w * cos_N - sin_w * sin_N * cos_I) - v * (sin_w * cos_N + cos_w * sin_N * cos_I),
    u * (cos_w * sin_N + sin_w * cos_N * cos_I) - v * (sin_w * sin_N - cos_w * cos_N * cos_I),
    u * sin_w * sin_I + v * cos_w * sin_I,
  ];
}


// Add two 3D vectors together.
function add(a, b) {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}


// atan2 to degrees, 0..360
function atan2(y, x) {
  return Math.atan2(-y, -x) * (180 / Math.PI) + 180;
}


// Convert Cartesian coordinates to longitude.
// FIXME: We could compute latitude and distance here, too, if we wanted.
function longitude(v) {
  return Math.atan2(-v[1], -v[0]) * (180 / Math.PI) + 180;
}


// Compute the positions of the planets for a given date.
// J. L. Simon et al. "Numerical expressions for precession formulae and
// mean elements for the Moon and the planets." Astronomy and Astrophysics
// 282 (1994): 663-683.
function planets(date) {
  const t = j2000c(date);

  const sun = kepler(
    1.0000010178,
    0.0167086342 - t *    0.00004203654,
    0,
    4.8950631131 + t *  628.33196540650,
    4.9381883009 + t *    0.03001023491,
    0,
  );

  const D = 5.1984665887 + t * 7771.37714559371;
  const F = 1.6279050815 + t * 8433.46615691637;
  const l = 2.3555557435 + t * 8328.69142571909;
  const M = 6.2400601269 + t *  628.30195517140;

  return {
    sun: longitude(sun),
    moon: longitude(
      kepler(
        0.0025628558 + t *    0.00000000003 + 22730e-9 * Math.cos(2 * D)
                                            -  4249e-9 * Math.cos(2 * D - l)
                                            -  1575e-9 * Math.cos(l)
                                            +  1458e-9 * Math.cos(2 * D - M)
                                            +  1210e-9 * Math.cos(2 * D + l),
        0.055545526  - t *    0.000000016   + 14216e-6 * Math.cos(2 * D - l)
                                            +  8551e-6 * Math.cos(2 * D - 2 * l)
                                            -  1383e-6 * Math.cos(l)
                                            +  1356e-6 * Math.cos(2 * D + l)
                                            -  1147e-6 * Math.cos(4 * D - 3 * l)
                                            -   914e-6 * Math.cos(4 * D - 2 * l)
                                            +   869e-6 * Math.cos(2 * D - M - l)
                                            -   627e-6 * Math.cos(2 * D)
                                            -   394e-6 * Math.cos(4 * D - 4 * l)
                                            +   282e-6 * Math.cos(2 * D - M - 2 * l)
                                            -   279e-6 * Math.cos(D - l)
                                            -   236e-6 * Math.cos(2 * l)
                                            +   231e-6 * Math.cos(4 * D)
                                            +   229e-6 * Math.cos(6 * D - 4 * l)
                                            -   201e-6 * Math.cos(2 * l - 2 * F), 
        0.0900012160 - t *    0.00000000039 + 23575e-7 * Math.cos(2 * D - 2 * F)
                                            -  1946e-7 * Math.cos(2 * D)
                                            +  1818e-7 * Math.cos(2 * F)
                                            +  1247e-7 * Math.cos(2 * l - 2 * F)
                                            +   968e-7 * Math.cos(2 * D - M - 2 * F),
        3.8103442782 + t * 8399.70911096274 - 16158e-6 * Math.sin(2 * D)
                                            +  5805e-6 * Math.sin(2 * D - l)
                                            -  3212e-6 * Math.sin(M)
                                            +  1921e-6 * Math.sin(l)
                                            -  1057e-6 * Math.sin(2 * D - M),
        1.4547885347 + t *   71.01768524366 - 26960e-5 * Math.sin(2 * D - l)
                                            - 16828e-5 * Math.sin(2 * D - 2 * l)
                                            -  4747e-5 * Math.sin(l)
                                            +  4550e-5 * Math.sin(4 * D - 3 * l)
                                            +  3639e-5 * Math.sin(4 * D - 2 * l)
                                            +  2578e-5 * Math.sin(2 * D + l)
                                            +  1689e-5 * Math.sin(4 * D - 4 * l)
                                            -  1657e-5 * Math.sin(2 * D - M - l)
                                            -  1227e-5 * Math.sin(6 * D - 4 * l)
                                            -  1152e-5 * Math.sin(2 * D)
                                            -  1006e-5 * Math.sin(2 * D - 3 * l)
                                            -   913e-5 * Math.sin(2 * l)
                                            -   842e-5 * Math.sin(6 * D - 5 * l)
                                            +   788e-5 * Math.sin(M)
                                            -   664e-5 * Math.sin(6 * D - 3 * l),
        2.1824391966 - t *   33.75704595363 -  2614e-5 * Math.sin(2 * D - 2 * F)
                                            -   262e-5 * Math.sin(M)
                                            -   214e-5 * Math.sin(2 * D)
                                            +   205e-5 * Math.sin(2 * F)
                                            -   140e-5 * Math.sin(2 * l - 2 * F),
      ),
    ),
    mercury: longitude(
      add(
        sun,
        kepler(
          0.3870983098,
          0.2056317526 + t *    0.00002040653,
          0.1222600741 + t *    0.00003179069,
          4.4026088425 + t * 2608.81470576909,
          1.3518643031 + t *    0.02716431730,
          0.8435332140 + t *    0.02070155118,
        ),
      ),
    ),
    venus: longitude(
      add(
        sun,
        kepler(
          0.7233298200,
          0.0067719164 - t *    0.00004776521,
          0.0592480270 + t *    0.00001751758,
          3.1761466970 + t * 1021.35294171618,
          2.2962197935 + t *    0.02447216844,
          1.3383170775 + t *    0.01572618080,
        ),
      ),
    ),
    mars: longitude(
      add(
        sun,
        kepler(
          1.5236793419 + t *    0.00000000003,
          0.0934006477 + t *    0.00009048438,
          0.0322838173 - t *    0.00001049081,
          6.2034761129 + t *  334.08562607873,
          5.8653575674 + t *    0.03213095395,
          0.8649518975 + t *    0.01347427507,
        ),
      ),
    ),
    jupiter: longitude(
      add(
        sun,
        kepler(
          5.2026032092 + t *    0.00000019132,
          0.0484979255 + t *    0.00016322542,
          0.0227462998 - t *    0.00009593223,
          0.5995471051 + t *   52.99348050854,
          0.2501267457 + t *    0.02814579341,
          1.7534346836 + t *    0.01781941774,
        ),
      ),
    ),
    saturn: longitude(
      add(
        sun,
        kepler(
          9.5549091915 - t *    0.00000213896,
          0.0555481426 - t *    0.00034664062,
          0.0434391294 - t *    0.00006520932,
          0.8740162840 + t *   21.35429658204,
          1.6241551868 + t *    0.03427410072,
          1.9838372649 + t *    0.01579288747,
        ),
      ),
    ),
  };
}

function angles(date, lat, lon) {
  const t = j2000(date);

  // https://aa.usno.navy.mil/faq/GAST
  const lst = 4.8949613 + 6.3003880989482 * t + lon * (Math.PI / 180);
  const sin_lst = Math.sin(lst);
  const cos_lst = Math.cos(lst);
  const ε = 0.409093 - 0.000000007 * t;
  const sin_ε = Math.sin(ε);
  const cos_ε = Math.cos(ε);

  return {
    // https://radixpro.com/a4a-start/the-ascendant/
    ascendant: atan2(
      cos_lst,
      -(sin_ε * Math.tan(lat * (Math.PI / 180)) + cos_ε * sin_lst),
    ),
    // https://radixpro.com/a4a-start/medium-coeli/
    midheaven: atan2(sin_lst, cos_lst * cos_ε),
  };
}
